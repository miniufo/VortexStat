//
package pack;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Predicate;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.Typhoon;


//
public class BatchFlow{
	//
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		return year==2014&&(Integer.parseInt(ty.getID())==1409);
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="E:/Data/VortexStat/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	
	//
	public static void main(String[] args){
		//batchForGribCtl(); System.exit(0);
		
		ConcurrentUtil.initDefaultExecutor(1);
		
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			//System.out.println(ty.getID());
			
			//MainProcess.extractGribData(ty,path+ds+"/"+year+"/"+ty.getID()+"/",29);System.out.println("\n\n\n");
			//MainProcess.verticalInterpolation(ty,path+ds+"/"+year+"/"+ty.getID()+"/",true);
			MainProcess.modelingWithoutBLH(ty,path+ds+"/"+year+"/"+ty.getID()+"/",true);
		});
		
		ConcurrentUtil.shutdown();
	}
	
	public static void batchForGribCtl(){
		StringBuilder sb=new StringBuilder();
		
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			String id=ty.getID();
			
			sb.append("grib2ctl.pl "+path+ds+"/"+year+"/"+id+"/"+id+"pl.grb  > "+path+ds+"/"+year+"/"+id+"/"+id+"pl.ctl \n");
			sb.append("grib2ctl.pl "+path+ds+"/"+year+"/"+id+"/"+id+"sfc.grb > "+path+ds+"/"+year+"/"+id+"/"+id+"sfc.ctl \n");
			sb.append("gribmap  -i "+path+ds+"/"+year+"/"+id+"/"+id+"pl.ctl \n");
			sb.append("gribmap  -i "+path+ds+"/"+year+"/"+id+"/"+id+"sfc.ctl \n");
			
			try{
				Files.deleteIfExists(Paths.get(path+ds+"/"+year+"/"+id+"/"+id+"sfc.ctl"));
				
			}catch(IOException e){ e.printStackTrace();System.exit(0);}
		});
		
		try(FileWriter fw=new FileWriter(path+ds+"/genGribCtl.bat")){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
}
