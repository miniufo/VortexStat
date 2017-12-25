//
package pack;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.Typhoon;
import org.json.JSONObject;
import org.python.core.PyFunction;
import org.python.core.PyMethod;
import org.python.core.PyObject;
import org.python.core.PyString;
import org.python.util.PythonInterpreter;


//
public final class DownoladData2{
	//
	static boolean prsData=true;
	
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		//return year==1985&&ty.getName().equalsIgnoreCase("Elena");
		return year==2012&&(Integer.parseInt(ty.getID())==1208);
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="E:/Data/VortexStat/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	private static final DateTimeFormatter fmt=DateTimeFormatter.ofPattern("yyyy-MM-dd").withLocale(Locale.ENGLISH);
	
	
	//
	public static void main(String[] args){
		List<KeyValue> data=all.stream().filter(cond).flatMap(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			
			String id=ty.getID();
			
			MDate tstr=new MDate(ty.getTime(0));
			MDate tend=new MDate(ty.getTime(ty.getTCount()-1));
			
			String date=tstr.toFormattedDate(fmt)+"/to/"+tend.toFormattedDate(fmt);
			String target=path+ds+"/"+year+"/"+id+"/"+ty.getID();
			
			if(prsData) return Stream.of(
				new KeyValue(id+"_sfc",request2DERA5(   date,target+"sfcERA5.nc")),
				new KeyValue(id+"_pl ",request3DPrsERA5(date,target+ "plERA5.nc"))
			);
			else return Stream.of(
				new KeyValue(id+"_pt ",request3DPTERA5( date,target+ "pt.grb"))
			);
			
		}).collect(Collectors.toList());
		
		JythonDownload(data);
	}
	
	static void JythonDownload(List<KeyValue> data){
		PythonInterpreter interp=new PythonInterpreter();
		interp.execfile("D:/Java/ecmwf-api-client-python/ecmwfapi/api.py");
		interp.exec("parseFunc = json.loads");
		interp.exec("server = ECMWFDataServer()");
		interp.exec("retrieve = server.retrieve");
		
		PyFunction parser=interp.get("parseFunc"   ,PyFunction.class);
		PyMethod retrieve=interp.get("retrieve",PyMethod.class);
		
		data.forEach(kv->{
			System.out.println("\n\n\nstart downloading "+kv.getKey());
			
			JSONObject json=kv.getVal();
			
			Path p=Paths.get(json.get("target").toString()).getParent();
			
			if(!Files.exists(p)){
				try{Files.createDirectories(p);}
				catch(IOException e){ e.printStackTrace(); System.exit(0);}
			}
			
			PyObject pyJson=parser.__call__(new PyString(json.toString()));
			retrieve.__call__(pyJson);
		});
		
		interp.close();
	}
	
	static JSONObject request3DPrsInterim(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("dataset" , "interim");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("repres" , "ll");	// sh/ll/gg
		request.put("levtype" , "pl");	// ml/sfc/pl/pv/pt/dp
		request.put("levelist" , "1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50");
		request.put("param"   , "z/t/q/u/v/w");	//
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/06:00:00/12:00:00/18:00:00");
		request.put("type"    , "an");	// an/fc
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "84/81/-15/222"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.3/0.3");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request3DPrsERA5(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("class  " , "ea");
		request.put("dataset" , "era5");
		request.put("expver"  , "1");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("type"    , "an");	// an/fc/4v
		request.put("repres"  , "ll");	// sh/ll/gg
		request.put("levtype" , "pl");	// ml/sfc/pl/pv/pt/dp
		request.put("levelist", "1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50");
		request.put("param"   , "z/t/q/u/v/w");	//
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/02:00:00/04:00:00/06:00:00/08:00:00/10:00:00/12:00:00/14:00:00/16:00:00/18:00:00/20:00:00/22:00:00");
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "54/81/-6/153"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.3/0.3");
		request.put("format"  , "netcdf");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request3DPTInterim(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("dataset" , "interim");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("repres" , "ll");	// sh/ll/gg
		request.put("levtype" , "pt");	// ml/sfc/pl/pv/pt/dp
		request.put("levelist" , "265/275/285/300/315/330/350/370/395/430");
		request.put("param"   , "u/v/d/vo/mont/pres/q");	//
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/06:00:00/12:00:00/18:00:00");
		request.put("type"    , "an");	// an/fc
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "84/81/-15/228"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.75/0.75");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request3DPTERA5(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("class  " , "ea");
		request.put("dataset" , "era5");
		request.put("expver"  , "1");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("type"    , "an");	// an/fc/4v
		request.put("repres" , "ll");	// sh/ll/gg
		request.put("levtype" , "pt");	// ml/sfc/pl/pv/pt/dp
		request.put("levelist" , "265/275/285/300/315/330/350/370/395/430");
		request.put("param"   , "u/v/d/vo/mont/pres/q");	//
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/02:00:00/04:00:00/06:00:00/08:00:00/10:00:00/12:00:00/14:00:00/16:00:00/18:00:00/20:00:00/22:00:00");
		request.put("type"    , "an");	// an/fc
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "84/81/-15/222"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.3/0.3");
		request.put("format"  , "netcdf");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request2DInterim(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("dataset" , "interim");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("repres" , "ll");	// sh/ll/gg
		request.put("levtype" , "sfc");	// ml/sfc/pl
		request.put("param"   , "z/msl/10u/10v/sp/sst");	//sshf/slhf/2t/
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/06:00:00/12:00:00/18:00:00");
		request.put("type"    , "an");	// an/fc
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "84/81/-15/222"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.3/0.3");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request2DERA5(String date,String fname){
		JSONObject  request = new JSONObject();
		
		request.put("class  " , "ea");
		request.put("dataset" , "era5");
		request.put("expver"  , "1");
		request.put("date"    , date);
		request.put("stream"  , "oper");
		request.put("type"    , "an");	// an/fc
		request.put("repres"  , "ll");	// sh/ll/gg
		request.put("levtype" , "sfc");	// ml/sfc/pl
		request.put("param"   , "z/msl/10u/10v/sp/sst");	//sshf/slhf/2t/
		request.put("step"    , "0");
		request.put("time"    , "00:00:00/02:00:00/04:00:00/06:00:00/08:00:00/10:00:00/12:00:00/14:00:00/16:00:00/18:00:00/20:00:00/22:00:00");
		if(ds==DataSets.NHC) request.put("area"    , "63/231/-9/357"); // north/west/south/east
		else request.put("area"    , "54/81/-6/153"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.3/0.3");
		request.put("format"  , "netcdf");
		request.put("target"  , fname);
		
		return request;
	}
	
	
	private static final class KeyValue{
		private String     key=null;
		private JSONObject val=null;
		
		//
		public KeyValue(String key,JSONObject val){
			this.key=key;
			this.val=val;
		}
		
		
		/*** getor and setor ***/
		public String getKey(){ return key;}
		
		public JSONObject getVal(){ return val;}
	}
}
