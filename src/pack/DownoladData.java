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

import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.Typhoon;

import org.ecmwf.DataServer;
import org.json.JSONObject;


//
public final class DownoladData{
	//
	static boolean prsData=true;
	
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		//return year==1985&&ty.getName().equalsIgnoreCase("Elena");
		return year==2015&&(Integer.parseInt(ty.getID())==1517);
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="E:/Data/VortexStat/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	private static final DateTimeFormatter fmt=DateTimeFormatter.ofPattern("yyyy-MM-dd").withLocale(Locale.ENGLISH);
	
	
	public static void main(String[] args){
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			
			MDate tstr=new MDate(ty.getTime(0));
			MDate tend=new MDate(ty.getTime(ty.getTCount()-1));
			
			System.out.println("\n\n\nstart downloading "+ty.getID());
			
			if(prsData){
				download2DVars(
					tstr.toFormattedDate(fmt)+"/to/"+tend.toFormattedDate(fmt),
					path+ds+"/"+year+"/"+ty.getID()+"/"+ty.getID()+"sfc.grb"
				);
				
				//download3DVarsPrs(
				//	tstr.toFormattedDate(fmt)+"/to/"+tend.toFormattedDate(fmt),
				//	path+ds+"/"+year+"/"+ty.getID()+"/"+ty.getID()+"pl.grb"
				//);
				
			}else{
				download3DVarsPT(
					tstr.toFormattedDate(fmt)+"/to/"+tend.toFormattedDate(fmt),
					path+ds+"/Theta/"+year+"/"+ty.getID()+"/"+ty.getID()+"pt.grb"
				);
			}
		});
	}
	
	static void download2DVars(String time,String fname){
		Path p=Paths.get(fname).getParent();
		
		if(!Files.exists(p)){
			try{Files.createDirectories(p);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		DataServer server = new DataServer();
		
		try{
			server.retrieve(request2D(time,fname));
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void download3DVarsPrs(String time,String fname){
		Path p=Paths.get(fname).getParent();
		
		if(!Files.exists(p)){
			try{Files.createDirectories(p);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		DataServer server = new DataServer();
		
		try{
			server.retrieve(request3DPrs(time,fname));
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void download3DVarsPT(String time,String fname){
		Path p=Paths.get(fname).getParent();
		
		if(!Files.exists(p)){
			try{Files.createDirectories(p);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		DataServer server = new DataServer();
		
		try{
			server.retrieve(request3DPT(time,fname));
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}
	
	static JSONObject request3DPrs(String date,String fname){
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
		request.put("grid"    , "0.75/0.75");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request3DPT(String date,String fname){
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
		else request.put("area"    , "84/81/-15/222"); // north/west/south/east  (1517 TC: 84/81/-15/228)
		request.put("resol"   , "av");
		request.put("grid"    , "0.75/0.75");
		request.put("target"  , fname);
		
		return request;
	}
	
	static JSONObject request2D(String date,String fname){
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
		request.put("grid"    , "0.75/0.75");
		request.put("target"  , fname);
		
		return request;
	}
}
