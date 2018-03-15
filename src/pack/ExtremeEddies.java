//
package pack;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Predicate;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;
import miniufo.mathsphysics.BWLabel;


//
public class ExtremeEddies{
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
		//convertVarToGS("htHFC",15/86400f);
		//convertVarToGS("htVFC",5/86400f  );
		//convertVarToGS("aamHFCR",30/86400f);
		//convertVarToGS("aamVFCR",15/86400f);
		//convertVarToGS("EPy",2e7f);
		//convertVarToGS("EPz",2e6f);
		convertVarToGS("EPDiv",30f);
	}
	
	static void convertVarToGS(String vname,float threshold){
		List<Sample> samples=mapTyphoonsToSamples(all,vname,threshold);
		
		plotSamplesOnGrids(samples,threshold,vname,"XY");
		plotSamplesOnGrids(samples,threshold,vname,"YZ");
		plotSamplesOnStns(samples,threshold,vname);
		
		System.out.println("finished "+vname+" (>"+threshold+")");
	}
	
	static void plotSamplesOnGrids(List<Sample> samples,float threshold,String vname,String plane){
		AtomicInteger counter=new AtomicInteger(0);
		
		StringBuilder sb=new StringBuilder();
		
		Comparator<Sample> compare=(s1,s2)->{return Float.compare(s2.weiCountP.value,s1.weiCountP.value);};
		Consumer<Sample> processSample=null;
		
		if("YZ".equalsIgnoreCase(plane)) processSample=(s)->addYZPlot(s,sb,counter,threshold,vname);
		if("XY".equalsIgnoreCase(plane)) processSample=(s)->addYXPlot(s,sb,counter);
		
		Predicate<Sample> upper=s->s.grdCountP.zPos<500;
		Predicate<Sample> splTD=s->s.type==TYPE.TD;
		Predicate<Sample> splTS=s->s.type==TYPE.TS;
		Predicate<Sample> splTY=s->s.type==TYPE.TY;
		Predicate<Sample> splEC=s->s.type==TYPE.EC;
		
		sb.append("'enable print "+path+"VortexStat/JMA/Stat/"+vname+"/Extreme"+plane+".gmf'\n\n");
		
		counter.set(0);
		samples.stream().filter(upper.and(splTD)).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter(upper.and(splTS)).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter(upper.and(splTY)).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter(upper.and(splEC)).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		sb.append("'disable print'\n");
		sb.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(path+"VortexStat/JMA/Stat/"+vname+"/Extreme"+plane+".gs")){ fw.write(sb.toString());}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void plotSamplesOnStns(List<Sample> samples,float threshold,String vname){
		AtomicInteger counter=new AtomicInteger(0);
		
		StringBuilder sb=new StringBuilder();
		
		Comparator<Sample> compare=(s1,s2)->{return Float.compare(s2.weiCountP.value,s1.weiCountP.value);};
		Consumer<Sample> processSample=(s)->addStnPlot(s,sb,counter,"200");
		
		sb.append("'open "+path+"VortexStat/JMA/etc/grid.ctl'\n");
		sb.append("'enable print "+path+"VortexStat/JMA/Stat/"+vname+"/Stn.gmf'\n\n");
		
		counter.set(0);
		samples.stream().filter((s)->{return s.type==TYPE.TD;}).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter((s)->{return s.type==TYPE.TS;}).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter((s)->{return s.type==TYPE.TY;}).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		counter.set(0);
		samples.stream().filter((s)->{return s.type==TYPE.EC;}).sorted(compare).limit(6).forEach(processSample);
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		sb.append("'disable print'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(path+"VortexStat/JMA/Stat/"+vname+"/stn.gs")){ fw.write(sb.toString());}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void addYZPlot(Sample s,StringBuilder sb,AtomicInteger counter,float threshold,String vname){
		int year=new MDate(s.r.getTime()).getYear();
		
		sb.append("'open "+path+"VortexStat/JMA/"+year+"/"+s.ID+"/"+s.ID+"_model.ctl'\n");
		sb.append("'set t "+s.tstep+"'\n");
		sb.append("'set lev 1000 75'\n");
		sb.append("'set zlog on'\n");
		sb.append("'set csmooth on'\n");
		sb.append("'set grads off'\n");
		sb.append("'setvpage "+StatUtil.getVPageString(2,3,counter.get())+"'\n");
		sb.append("'setlopts 6 0.16 250 100'\n");
		sb.append("'set xaxis 0 1500 250'\n");
		sb.append("'set gxout shade2'\n");
		sb.append("'set cmin "+threshold+"'\n");
		sb.append("'d "+vname+"*86400'\n");
		sb.append("'cbarn'\n");
		sb.append("'set gxout contour'\n");
		sb.append("'d "+vname+"*86400'\n");
		sb.append("'draw title "+s.toString()+" ("+s.grdCountP.value+", "+s.weiCountP.value+")'\n");
		sb.append("'close 1'\n\n");
		
		counter.incrementAndGet();
	}
	
	static void addYXPlot(Sample s,StringBuilder sb,AtomicInteger counter){
		int year=new MDate(s.r.getTime()).getYear();
		
		sb.append("'open "+path+"VortexStat/JMA/"+year+"/"+s.ID+"/"+s.ID+".ctl'\n");
		sb.append("'set t "+s.tstep+"'\n");
		sb.append("'set lev 200'\n");
		sb.append("'set csmooth on'\n");
		sb.append("'set grads off'\n");
		sb.append("'set lon "+(s.r.getXPos()-20)+" "+(s.r.getXPos()+20)+"'\n");
		sb.append("'set lat "+(s.r.getYPos()-18)+" "+(s.r.getYPos()+18)+"'\n");
		sb.append("'setvpage "+StatUtil.getVPageString(2,3,counter.get())+"'\n");
		sb.append("'setlopts 5 0.16 5 5'\n");
		sb.append("'set gxout shade2'\n");
		sb.append("'color 213 231 1 -kind grainbow'\n");
		sb.append("'d t'\n");
		sb.append("'cbarn'\n");
		sb.append("'set arrowhead -0.35'\n");
		sb.append("'set arrscl 0.5 60'\n");
		sb.append("'d skip(u,2);v'\n");
		sb.append("'draw title "+s.toString()+" ("+s.grdCountP.value+", "+s.weiCountP.value+")'\n");
		sb.append("'close 1'\n\n");
		
		counter.incrementAndGet();
	}
	
	static void addStnPlot(Sample s,StringBuilder sb,AtomicInteger counter,String lev){
		int year=new MDate(s.r.getTime()).getYear();
		
		sb.append("'open "+path+"VortexStat/JMA/"+year+"/"+s.ID+"/"+s.ID+"_stn.ctl'\n");
		sb.append("'set csmooth on'\n");
		sb.append("'set grads off'\n");
		sb.append("'set gxout shade2'\n");
		sb.append("'set lon "+(s.r.getXPos()-15)+" "+(s.r.getXPos()+15)+"'\n");
		sb.append("'set lat "+(s.r.getYPos()-15)+" "+(s.r.getYPos()+15)+"'\n");
		sb.append("'define uu=oacres(u,u.2(t="+s.tstep+",lev="+lev+"))'\n");
		sb.append("'define vv=oacres(u,v.2(t="+s.tstep+",lev="+lev+"))'\n");
		sb.append("var='gava'\n");
		sb.append("'setvpage "+StatUtil.getVPageString(2,3,counter.get())+"'\n");
		sb.append("'setlopts 5 0.18 5 5'\n");
		sb.append("'define var=oacres(u,'var'.2(t="+s.tstep+",lev="+lev+"))'\n");
		sb.append("'d var'\n");
		sb.append("'cbarn'\n");
		sb.append("'set arrowhead -0.35'\n");
		sb.append("'set arrscl 0.5 60'\n");
		sb.append("'d skip(uu,5);vv'\n");
		sb.append("'draw title 'var\n");
		sb.append("'close 2'\n\n");
		
		counter.incrementAndGet();
	}
	
	static List<Sample> mapTyphoonsToSamples(List<Typhoon> all,String vname,float threshold){
		int totalLength=all.stream().filter(cond).mapToInt(Typhoon::getTCount).sum();
		
		List<Sample> samples=new ArrayList<>(totalLength);
		
		all.stream().filter(cond).forEach(ty->{
			int year=new MDate(ty.getTime(0)).getYear();
			
			DiagnosisFactory df=DiagnosisFactory.parseFile(path+"VortexStat/JMA/"+year+"/"+ty.getID()+"/"+ty.getID()+"_model.ctl");
			DataDescriptor dd=df.getDataDescriptor();
			df.setPrinting(false);
			
			int t=ty.getTCount(),z=dd.getZCount(),y=dd.getYCount();
			
			Variable v=df.getVariables(new Range("",dd),false,vname)[0];
			
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
				
				double yRatio=Math.toRadians(0.15)*SpatialModel.EARTH_RADIUS/1000.0;
				
				Metric[] tmp=countingOneSample(buf,v.getUndef(),(f)->{return f>threshold;});
				
				s.grdCountP=tmp[0]; tmp[0].yPos*=yRatio; tmp[0].zPos=1000+tmp[0].zPos*dd.getDZDef()[0];
				s.weiCountP=tmp[1]; tmp[1].yPos*=yRatio; tmp[1].zPos=1000+tmp[1].zPos*dd.getDZDef()[0];
				s.weiMeanP =tmp[2]; tmp[2].yPos*=yRatio; tmp[2].zPos=1000+tmp[2].zPos*dd.getDZDef()[0];
				
				tmp=countingOneSample(buf,v.getUndef(),(f)->{return f<-threshold;});
				s.grdCountN=tmp[0]; tmp[0].yPos*=yRatio; tmp[0].zPos=1000+tmp[0].zPos*dd.getDZDef()[0];
				s.weiCountN=tmp[1]; tmp[1].yPos*=yRatio; tmp[1].zPos=1000+tmp[1].zPos*dd.getDZDef()[0];
				s.weiMeanN =tmp[2]; tmp[2].yPos*=yRatio; tmp[2].zPos=1000+tmp[2].zPos*dd.getDZDef()[0];
				
				samples.add(s);
			}
		});
		
		return samples;
	}
	
	static Metric[] countingOneSample(float[][] data,float undef,Predicate<Float> cond){
		int z=data.length,y=data[0].length;
		
		float[][] bw=StatUtil.mapToBWLabelData(data,undef,cond);
		
		BWLabel lb=new BWLabel(bw);
		lb.connComponentlabelling(8);
		
		int label=lb.getLabelsSortedByCount()[0];
		
		if(label==0){
			return new Metric[]{new Metric(0,0,0),new Metric(0,0,0),new Metric(0,0,0)};
			
		}else{
			int count=0;
			double weight=0,c_yPos=0,c_zPos=0,w_yPos=0,w_zPos=0;
			
			float[][] labels=lb.getLabelData();
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) if(labels[k][j]==label){
				count++;   weight+=  data[k][j];
				c_yPos+=j; w_yPos+=j*data[k][j];
				c_zPos+=k; w_zPos+=k*data[k][j];
			}
			
			c_yPos/=count; w_yPos/=weight;
			c_zPos/=count; w_zPos/=weight;
			
			return new Metric[]{
				new Metric((float)count ,(float)c_yPos,(float)c_zPos),
				new Metric((float)weight,(float)w_yPos,(float)w_zPos),
				new Metric((float)(weight/count),(float)w_yPos,(float)w_zPos)
			};
		}
	}
	
	
	static final class Sample{
		//
		int tstep=0;	// start from 1
		
		Metric grdCountP=null;
		Metric weiCountP=null;
		Metric weiMeanP =null;
		
		Metric grdCountN=null;
		Metric weiCountN=null;
		Metric weiMeanN =null;
		
		TYPE type=TYPE.OTHERS;
		
		String name=null;
		String ID=null;
		
		Record r=null;
		
		public String toString(){ return ID+"("+tstep+")";}
	}
	
	static final class Metric{
		//
		float value=0;
		float zPos =0; // hPa
		float yPos =0; // km
		
		//
		public Metric(float value,float yPos,float zPos){
			this.value=value;
			this.yPos =yPos;
			this.zPos =zPos;
		}
	}
}
