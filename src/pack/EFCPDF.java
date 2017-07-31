package pack;
//
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;
import miniufo.mathsphysics.BWLabel;
import pack.ExtremeEddies.Metric;
import pack.ExtremeEddies.Sample;


//
public class EFCPDF{
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
		float threshold=10;
		String vname="aamHFCR";
		
		List<Sample> samples=mapTyphoonsToSamples(all,vname,threshold);
		
		System.out.println(
			"grdCountP   yPos   zPos      weiCountP   yPos   zPos   weiMeanP "+
			"grdCountN   yPos   zPos      weiCountN   yPos   zPos   weiMeanN ID(tstep)   lon    lat"
		);
		
		samples.stream().filter(s->{return s.type==TYPE.EC;}).forEach(s->{
			System.out.println(String.format(
				"%9.1f %6.1f %6.1f %14.7f %6.1f %6.1f %10.6f %9.1f %6.1f %6.1f %14.7f %6.1f %6.1f %10.6f %s %6.1f %6.1f",
				s.grdCountP.value,s.grdCountP.yPos,s.grdCountP.zPos,
				s.weiCountP.value,s.weiCountP.yPos,s.weiCountP.zPos,s.weiMeanP.value,
				s.grdCountN.value,s.grdCountN.yPos,s.grdCountN.zPos,
				s.weiCountN.value,s.weiCountN.yPos,s.weiCountN.zPos,s.weiMeanN.value,
				s,s.r.getLon(),s.r.getLat()
			));
		});
		
		System.out.println("finished "+vname+" ("+threshold+")");
	}
	
	static float[] getCountData(List<Sample> samples,boolean positive){
		Predicate<Sample> cond=s->{return (s.type==TYPE.EC)&&(s.grdCountP.zPos<500);};
		
		double[] tmp=samples.stream().filter(cond).mapToDouble(s->s.grdCountP.value).toArray();
		
		float[] re=new float[tmp.length];
		
		for(int i=0,I=tmp.length;i<I;i++) re[i]=(float)tmp[i];
		
		return re;
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
			
			Variable v=df.getVariables(new Range("",dd),false,vname)[0].multiplyEq(86400);
			
			if(ty.getTime(0)!=dd.getTimes()[0]||ty.getTime(t-1)!=dd.getTimes()[t-1])
			throw new IllegalArgumentException("time ranges are not the same");
			
			for(int l=0;l<t;l++){
				Sample s=new Sample();
				
				float[][] buf=new float[z][y];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=v.getData()[k][j][0][l];
				
				s.name=ty.getName();
				s.type=Typhoon.getType(ty.getRecord(l).getDataValue(4));
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
}
