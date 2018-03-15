package pack;
//
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

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


//
public class Discriminant{
	//
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		
		return year==2013||year==2014||year==2015;
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="E:/Data/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	//
	public static void main(String[] args){
		int totalLength=all.stream().filter(cond).mapToInt(Typhoon::getTCount).sum();
		
		System.out.println("total length is "+totalLength);
		
		List<Sample> samples=new ArrayList<>(totalLength);
		
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			int count=ty.getTCount();
			
			System.out.println("For TC "+ty.getID());
			
			String fname=path+"VortexStat/JMA/"+year+"/"+ty.getID()+"/"+ty.getID()+"_model.ctl";
			
			DiagnosisFactory df=DiagnosisFactory.parseFile(fname);
			DataDescriptor dd=df.getDataDescriptor();
			
			if(ty.getTime(0)!=dd.getTimes()[0]||ty.getTime(count-1)!=dd.getTimes()[count-1])
			throw new IllegalArgumentException("time ranges are not the same");
			
			samples.addAll(Arrays.asList(mapToSamples(ty,df)));
		});
		
		Variable[] all=mapToVariables(samples,s->{return true;});
		Variable[] TD =mapToVariables(samples,s->{return s.r.getDataValue(Typhoon.Type)==0;});
		Variable[] TS =mapToVariables(samples,s->{return s.r.getDataValue(Typhoon.Type)==1;});
		Variable[] TY =mapToVariables(samples,s->{return s.r.getDataValue(Typhoon.Type)==2;});
		Variable[] EC =mapToVariables(samples,s->{return s.r.getDataValue(Typhoon.Type)==3;});
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("E:/Data/VortexStat/JMA/Stat/DiscriminantAll.dat");
		cdws.writeData(ArrayUtil.concatAll(Variable.class,all,TD,TS,TY,EC)); cdws.closeFile();
	}
	
	static Variable[] mapToVariables(List<Sample> samples,Predicate<Sample> condSmpl){
		Sample one=samples.get(0);
		
		int z=one.A.length;
		int y=one.A[0].length;
		
		Range r=new Range(1,z,y,1);
		
		Variable Acount =new Variable("ac",true,r);  Acount.setUndef(one.undef);
		Variable A2count=new Variable("a2",true,r); A2count.setUndef(one.undef);
		Variable Ccount =new Variable("cc",true,r);  Ccount.setUndef(one.undef);
		Variable Dcount =new Variable("dc",true,r);  Dcount.setUndef(one.undef);
		Variable Ecount =new Variable("ec",true,r);  Ecount.setUndef(one.undef);
		Variable Tcount =new Variable("tc",true,r);  Tcount.setUndef(one.undef);
		
		float[][][][] Acdata= Acount.getData();  Acount.setComment("D < 0 caused by A < 0 && C > 0");
		float[][][][] A2data=A2count.getData(); A2count.setComment("A < 0 && C < 0");
		float[][][][] Ccdata= Ccount.getData();  Ccount.setComment("D < 0 caused by A > 0 && C < 0");
		float[][][][] Dcdata= Dcount.getData();  Dcount.setComment("D < 0 caused by A > 0 && C > 0");
		float[][][][] Ecdata= Ecount.getData();  Ecount.setComment("D < 0 caused by A < 0 && C < 0");
		float[][][][] Tcdata= Tcount.getData();  Tcount.setComment("D < 0 caused by all");
		
		samples.stream().filter(condSmpl).forEach(s->{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float U=s.undef;
				float a=s.A[k][j];
				float b=s.B[k][j];
				float c=s.C[k][j];
				float d=s.D[k][j];
				
				if(a!=U&&b!=U&&c!=U&&d!=U){
					if(a<0&&c<0) A2data[0][k][j][0]++;
					
					if(d<0){
						Tcdata[0][k][j][0]++;
						
						if(a<0&&c>0) Acdata[0][k][j][0]++;
						if(a>0&&c<0) Ccdata[0][k][j][0]++;
						if(a>0&&c>0) Dcdata[0][k][j][0]++;
						if(a<0&&c<0) Ecdata[0][k][j][0]++;
					}
				}
			}
		});
		
		return new Variable[]{Acount,A2count,Ccount,Dcount,Ecount,Tcount};
	}
	
	static Sample[] mapToSamples(Typhoon ty,DiagnosisFactory df){
		DataDescriptor dd=df.getDataDescriptor();
		
		Range r=new Range("",dd);
		
		int z=dd.getZCount(),y=dd.getYCount();
		
		Variable[] vs=df.getVariables(r,false,"sm","BsinB","CsinB");
		df.setPrinting(false);
		
		Variable A=vs[0];
		Variable B=vs[1];
		Variable C=vs[2];
		Variable D=A.multiply(C).minusEq(B.multiply(B));
		
		Sample[] samples=new Sample[ty.getTCount()];
		
		for(int l=0,L=ty.getTCount();l<L;l++){
			samples[l]=new Sample();
			
			samples[l].undef=A.getUndef();
			
			float[][] Adata=new float[z][y];
			float[][] Bdata=new float[z][y];
			float[][] Cdata=new float[z][y];
			float[][] Ddata=new float[z][y];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				Adata[k][j]=A.getData()[k][j][0][l];
				Bdata[k][j]=B.getData()[k][j][0][l];
				Cdata[k][j]=C.getData()[k][j][0][l];
				Ddata[k][j]=D.getData()[k][j][0][l];
			}
			
			samples[l].r=ty.getRecord(l);
			samples[l].A=Adata;
			samples[l].B=Bdata;
			samples[l].C=Cdata;
			samples[l].D=Ddata;
		}
		
		return samples;
	}
	
	
	static final class Sample{
		//
		float undef=Float.NaN;
		
		Record r=null;
		
		float[][] A=null;
		float[][] B=null;
		float[][] C=null;
		float[][] D=null;
	}
}
