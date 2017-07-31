package pack;
//
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
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
import miniufo.util.GridDataFetcher;


//
public class UnderEstimates{
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
			
			String fname=path+"VortexStat/JMA/"+year+"/"+ty.getID()+"/"+ty.getID()+".ctl";
			
			DataDescriptor dd=DiagnosisFactory.getDataDescriptor(fname);
			
			if(ty.getTime(0)!=dd.getTimes()[0]||ty.getTime(count-1)!=dd.getTimes()[count-1])
			throw new IllegalArgumentException("time ranges are not the same");
			
			samples.addAll(Arrays.asList(mapToSamples(ty,dd)));
		});
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("E:/Data/VortexStat/JMA/Stat/SLPCompare.dat");
		cdws.writeData(mapToVariables(samples)); cdws.closeFile();
	}
	
	static Variable[] mapToVariables(List<Sample> samples){
		Variable[] re=new Variable[4];
		
		Range r=new Range(samples.size(),1,1,1);
		
		re[0]=new Variable("mslp",false,r); re[0].setUndef(-9999); re[0].setComment("mean sea-level pressure (hPa)");
		re[1]=new Variable("prs" ,false,r); re[1].setUndef(-9999); re[1].setComment("best-track sea-level pressure (hPa)");
		re[2]=new Variable("type",false,r); re[2].setUndef(-9999); re[2].setComment("type number (1)");
		re[3]=new Variable("lats",false,r); re[3].setUndef(-9999); re[3].setComment("latitude (degree)");
		
		float[] r0=re[0].getData()[0][0][0];
		float[] r1=re[1].getData()[0][0][0];
		float[] r2=re[2].getData()[0][0][0];
		float[] r3=re[3].getData()[0][0][0];
		
		int ptr=0;
		for(Sample s:samples){
			r0[ptr]=s.mslp;
			r1[ptr]=s.r.getDataValue(3); // prs
			r2[ptr]=s.type.ordinal();
			r3[ptr]=s.r.getLat();
		}
		
		return re;
	}
	
	static Sample[] mapToSamples(Typhoon ty,DataDescriptor dd){
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		Sample[] samples=new Sample[ty.getTCount()];
		
		float[][][] buf=gdf.prepareXYTBuffer("mslp",1);
		
		for(int l=0,L=ty.getTCount();l<L;l++){
			samples[l]=new Sample();
			
			Record r=ty.getRecord(l);
			
			float lon=r.getLon();
			float lat=r.getLat();
			long  tim=r.getTime();
			
			samples[l].r   =r;
			samples[l].mslp=gdf.fetchXYTBuffer(lon,lat,tim,buf);
			samples[l].type=ty.getTypes()[l];
		}
		
		return samples;
	}
	
	
	static final class Sample{
		//
		float mslp=0;
		
		TYPE type=TYPE.OTHERS;
		
		Record r=null;
	}
}
