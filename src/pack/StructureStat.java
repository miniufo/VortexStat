package pack;
//
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;

import miniufo.application.statisticsModel.StatisticsBasicAnalysisMethods;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.util.TicToc;


//
public class StructureStat{
	//
	static final int[] count=new int[]{23,26,23,21,29,23,23,24,22,22,14,21,20};
	
	static final Predicate<Record> cond=id->{ return true;};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="/lustre/home/qianyk/Data/";
	static final String prefix="_spl_2nd";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons(path+"/Typhoons/JMA/JMA.txt","",ds);
	
	//
	public static void main(String[] args){
		//generateDataGS(); System.exit(0);
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"/VortexStat/JMA/StructureStatAll"+prefix+".ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Range r=new Range("",dd);
		
		String[] vars=new String[]{
			"ut","vr","w","t","q","a","z","the","gdwm","sigma","bsinb","csinb",
			"lh","fsfc","edyhadv","edyangmht","edyangmvt",
			"dhff","edyhaff","adhff","edyadvff","edyconff","frictff","tiltff"
		};
		
		Variable[] vmeans=new Variable[vars.length];
		Variable[] vstds =new Variable[vars.length];
		
		for(int i=0,I=vars.length;i<I;i++){
			TicToc.tic("computing "+vars[i]);
			Variable v=df.getVariables(r,vars[i])[0];
			
			 vstds[i]=StatisticsBasicAnalysisMethods.cTStandardDeviation(v);
			vmeans[i]=v.averageAlong(Dimension.T);
			TicToc.toc(TimeUnit.MINUTES);
		}
		
		DataWrite dw=null;
		dw=DataIOFactory.getDataWrite(dd,path+"/VortexStat/JMA/StructMeanAll"+prefix+".dat");
		dw.writeData(dd,vmeans); dw.closeFile();
		
		dw=DataIOFactory.getDataWrite(dd,path+"/VortexStat/JMA/StructSTDAll"+prefix+".dat");
		dw.writeData(dd,vstds); dw.closeFile();
	}
	
	static void generateDataGS(){
		int cc=0;
		
		StringBuilder sb=new StringBuilder();
		
		sb.append("'set gxout fwrite'\n");
		sb.append("'set fwrite "+path+"VortexStat/JMA/StructureStatAll"+prefix+".dat'\n\n");
		
		for(int yy=2000;yy<2013;yy++)
		for(int i=1;i<=count[yy-2000];i++){
			String id=String.format(String.valueOf(yy).substring(2)+"%02d",i);
			
			List<Typhoon> ls=AccessBestTrack.getTyphoons(all,"id="+id);
			
			if(ls.size()!=1) throw new IllegalArgumentException("find more than one TC");
			
			Typhoon ty=ls.get(0);
			
			sb.append("'open "+path+"VortexStat/JMA/"+yy+"/"+id+"/"+id+"_model"+prefix+".ctl'\n");
			
			for(int l=0,L=ty.getTCount();l<L;l++){
				Record r=ty.getRecord(l);
				
				if(cond.test(r)){
					sb.append("'set t '"+(l+1)+"\n");
					sb.append("'d ut'\n");
					sb.append("'d vr'\n");
					sb.append("'d w'\n");
					sb.append("'d t'\n");
					sb.append("'d q'\n");
					sb.append("'d a'\n");
					sb.append("'d z'\n");
					sb.append("'d the'\n");
					sb.append("'d gdwM'\n");
					sb.append("'d sigma'\n");
					sb.append("'d BsinB'\n");
					sb.append("'d CsinB'\n");
					sb.append("'d lh'\n");
					sb.append("'d fsfc'\n");
					sb.append("'d edyhadv'\n");
					sb.append("'d edyangmht'\n");
					sb.append("'d edyangmvt'\n");
					sb.append("'d dhff'\n");
					sb.append("'d edyhaff'\n");
					sb.append("'d adhff'\n");
					sb.append("'d edyadvff'\n");
					sb.append("'d edyconff'\n");
					sb.append("'d frictff'\n");
					sb.append("'d tiltff'\n");
					
					cc++;
				}
			}
			
			sb.append("'close 1'\n\n");
		}
		
		sb.append("'disable fwrite'\n");
		sb.append("'reinit'\n");
		
		System.out.println(" there are "+cc+" TC samples");
		
		try(FileWriter fw=new FileWriter(path+"VortexStat/JMA/StructureStatData"+prefix+".gs")){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
}
