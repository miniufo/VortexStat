package pack;
//
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;
import miniufo.application.statisticsModel.StatisticsApplication;
import miniufo.basic.ArrayUtil;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;


//
public class EFCRVProfile{
	//
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		
		return year==2013||year==2014||year==2015;
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="F:/Data/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	
	//
	public static void main(String[] args){
		Variable[][] r1=processOneVariable("aamHFCR");
		Variable[][] r2=processOneVariable("aamVFCR");
		Variable[][] r3=processOneVariable("htHFC"  );
		Variable[][] r4=processOneVariable("htVFC"  );
		
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor(path+"VortexStat/JMA/2013/1301/1301_model.ctl");
		
		CtlDataWriteStream cdws=null;
		cdws=new CtlDataWriteStream(path+"VortexStat/JMA/Stat/zProfile.dat");
		cdws.writeData(dd,ArrayUtil.concatAll(Variable.class,r1[0],r2[0],r3[0],r4[0])); cdws.closeFile();
		cdws=new CtlDataWriteStream(path+"VortexStat/JMA/Stat/yProfile.dat");
		cdws.writeData(dd,ArrayUtil.concatAll(Variable.class,r1[1],r2[1],r3[1],r4[1])); cdws.closeFile();
	}
	
	static Variable[][] processOneVariable(String vname){
		List<Sample> samples=mapTyphoonsToSamples(all,vname);
		
		Variable[] rTD=sumUpAllSamples(samples,vname,TYPE.TD);
		Variable[] rTS=sumUpAllSamples(samples,vname,TYPE.TS);
		Variable[] rTY=sumUpAllSamples(samples,vname,TYPE.TY);
		Variable[] rEC=sumUpAllSamples(samples,vname,TYPE.EC);
		
		return new Variable[][]{
			{rTD[0],rTD[1],rTS[0],rTS[1],rTY[0],rTY[1],rEC[0],rEC[1]},
			{rTD[2],rTD[3],rTS[2],rTS[3],rTY[2],rTY[3],rEC[2],rEC[3]}
		};
	}
	
	static Variable[] sumUpAllSamples(List<Sample> samples,String vname,TYPE type){
		int z=samples.get(0).zProfile.length;
		int y=samples.get(0).yProfile.length;
		
		AtomicInteger counter=new AtomicInteger(0);
		
		float undef=-9999;
		
		double[] zP=new double[z];
		double[] zS=new double[z];
		double[] yP=new double[y];
		double[] yS=new double[y];
		
		samples.stream().filter(s->{return s.type==type;}).forEachOrdered(s->{
			for(int k=0;k<z;k++){ zP[k]+=s.zProfile[k]; zS[k]+=s.zStdErr[k];}
			for(int j=0;j<y;j++){ yP[j]+=s.yProfile[j]; yS[j]+=s.yStdErr[j];}
			counter.incrementAndGet();
		});
		
		if(counter.get()!=0){
			for(int k=0;k<z;k++){ zP[k]/=counter.get(); zS[k]/=counter.get();}
			for(int j=0;j<y;j++){ yP[j]/=counter.get(); yS[j]/=counter.get();}
			
		}else{
			for(int k=0;k<z;k++){ zP[k]=undef; zS[k]=undef;}
			for(int j=0;j<y;j++){ yP[j]=undef; yS[j]=undef;}
		}
		
		Variable zProf=new Variable(vname+"zP"+type,true,new Range(1,z,1,1));
		Variable zStdE=new Variable(vname+"zS"+type,true,new Range(1,z,1,1));
		Variable yProf=new Variable(vname+"yP"+type,true,new Range(1,1,y,1));
		Variable yStdE=new Variable(vname+"yS"+type,true,new Range(1,1,y,1));
		
		zProf.setUndef(undef); zStdE.setUndef(undef);
		yProf.setUndef(undef); yStdE.setUndef(undef);
		
		float[][][] zpdata=zProf.getData()[0];
		float[][][] zsdata=zStdE.getData()[0];
		float[][][] ypdata=yProf.getData()[0];
		float[][][] ysdata=yStdE.getData()[0];
		
		for(int k=0;k<z;k++) zpdata[k][0][0]=(float)zP[k];
		for(int k=0;k<z;k++) zsdata[k][0][0]=(float)zS[k];
		for(int j=0;j<y;j++) ypdata[0][j][0]=(float)yP[j];
		for(int j=0;j<y;j++) ysdata[0][j][0]=(float)yS[j];
		
		return new Variable[]{zProf,zStdE,yProf,yStdE};
	}
	
	static List<Sample> mapTyphoonsToSamples(List<Typhoon> all,String vname){
		int totalLength=all.stream().filter(cond).mapToInt(Typhoon::getTCount).sum();
		
		List<Sample> samples=new ArrayList<>(totalLength);
		
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			
			DiagnosisFactory df=DiagnosisFactory.parseFile(path+"VortexStat/JMA/"+year+"/"+ty.getID()+"/"+ty.getID()+"_model.ctl");
			DataDescriptor dd=df.getDataDescriptor();
			df.setPrinting(false);
			
			int t=ty.getTCount(),z=dd.getZCount(),y=dd.getYCount();
			
			Variable v=df.getVariables(new Range("",dd),false,vname)[0].multiplyEq(86400);
			
			if(ty.getTime(0)!=dd.getTimes()[0]||ty.getTime(t-1)!=dd.getTimes()[t-1])
			throw new IllegalArgumentException("time ranges are not the same");
			
			for(int l=0;l<t;l++){
				Sample s=new Sample();
				
				float[][] buf=new float[z][y];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=v.getData()[k][j][0][l];
				
				s.name=ty.getName();
				s.type=Typhoon.getType(ty.getRecord(l).getDataValue(Typhoon.Type));
				s.tstep=l+1;
				s.ID=ty.getID();
				s.r =ty.getRecord(l);
				
				float[][] re=profileOneSampleData(buf,v.getUndef());
				
				s.zProfile=re[0]; s.zStdErr=re[1];
				s.yProfile=re[2]; s.yStdErr=re[3];
				
				samples.add(s);
			}
		});
		
		return samples;
	}
	
	static float[][] profileOneSampleData(float[][] data,float undef){
		int z=data.length,y=data[0].length;
		
		float[] zProfile=new float[z];
		float[] zStdErr =new float[z];
		float[] yProfile=new float[y];
		float[] yStdErr =new float[y];
		
		for(int k=0;k<z;k++){
			float[] buf=new float[y];
			
			for(int j=0;j<y;j++) buf[j]=data[k][j];
			
			zProfile[k]=StatisticsApplication.cArithmeticMean(buf,undef);
			zStdErr [k]=StatisticsApplication.cStandardError (buf,undef);
			
		}
		
		for(int j=0;j<y;j++){
			float[] buf=new float[z];
			
			for(int k=0;k<z;k++) buf[k]=data[k][j];
			
			yProfile[j]=StatisticsApplication.cArithmeticMean(buf,undef);
			yStdErr [j]=StatisticsApplication.cStandardError (buf,undef);
		}
		
		return new float[][]{zProfile,zStdErr,yProfile,yStdErr};
	}
	
	
	static final class Sample{
		//
		int tstep=0;	// start from 1
		
		float[] zProfile=null;
		float[] zStdErr =null;
		float[] yProfile=null;
		float[] yStdErr =null;
		
		TYPE type=TYPE.OTHERS;
		
		String name=null;
		String ID=null;
		
		Record r=null;
		
		public String toString(){ return ID+"("+tstep+")";}
	}
}
