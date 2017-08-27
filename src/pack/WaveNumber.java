//
package pack;

import java.io.FileWriter;
import java.io.IOException;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.function.Predicate;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.basic.ThermoDynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.application.statisticsModel.FilterMethods;
import miniufo.application.statisticsModel.StatisticsBasicAnalysisMethods;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CsmDataWriteStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


//
public class WaveNumber{
	//
	static final Predicate<Typhoon> cond=ty->{
		int year=new MDate(ty.getTime(0)).getYear();
		//return "0420".equalsIgnoreCase(ty.getID()); // Haima (0420)
		return year==1985&&ty.getName().equalsIgnoreCase("Elena"); // Elena (19855)
		//return (year==2013||year==2014||year==2015)&&Integer.parseInt(ty.getID())==1326;
	};
	
	static final DataSets ds=DataSets.NHC;
	static final String path="F:/Data/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	
	//
	public static void main(String[] args){
		all.stream().filter(cond).forEach(ty->{
			processOneTyphoon(ty);
			generateStationGS(ty);
		});
	}
	
	
	static void processOneTyphoon(Typhoon ty){
		int year=new MDate(ty.getTime(0)).getYear();
		
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			ty.toCSMString(path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/"+ty.getID()+".ctl",144,91,77,0.15f,-12.5f,1000)
		).getDataDescriptor();
		
		// first four are uvwT, and followed by htHFC htVFC aamHFCR aamVFCR
		Variable[] vs=getEddyVariables(ty,csd);
		
		CsmDataWriteStream cdws=new CsmDataWriteStream(path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/WNStation_UVWT.dat");
		cdws.writeData(csd,vs[0],vs[1],vs[2],vs[3]); cdws.closeFile();
		
		List<Variable> ls=new ArrayList<>();
		
		for(Variable v:new Variable[]{vs[4],vs[5],vs[6],vs[7]}){
			Variable[] comps=divideComponents(v,5);
			
			cdws=new CsmDataWriteStream(path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/WNStation_"+v.getName()+".dat");
			cdws.writeData(csd,comps); cdws.closeFile();
			
			Variable[] vars=getVariances(v,comps);
			
			for(Variable vv:vars) ls.add(vv);
		}
		
		CtlDataWriteStream ctws=new CtlDataWriteStream(path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/WNVariance.dat");
		ctws.writeData(csd.getCtlDescriptor(),ArrayUtil.concatAll(Variable.class,
			new Variable[]{vs[0].anomalizeX(),vs[1].anomalizeX(),vs[2].anomalizeX(),vs[3].anomalizeX()},ls.toArray(new Variable[ls.size()])
		)); ctws.closeFile();
	}
	
	static void generateStationGS(Typhoon ty){
		int year=new MDate(ty.getTime(0)).getYear();
		
		StringBuilder sb=new StringBuilder();
		
		sb.append("var='aamHFCR'\n");
		sb.append("'open "+path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/"+ty.getID()+".ctl'\n");
		sb.append("'open "+path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/WNStation_'var'.ctl'\n");
		sb.append("'open "+path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/WNStation_UVWT.ctl'\n");
		sb.append("'enable print "+path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/plot/WNStation.gmf'\n\n");
		
		sb.append("olon="+ArrayUtil.allToString(ty.getXPositions(),"'"," ","'")+"\n");
		sb.append("olat="+ArrayUtil.allToString(ty.getYPositions() ,"'"," ","'")+"\n");
		sb.append("'set lev 200'\n\n");
		
		sb.append("tt=1\n");
		sb.append("while(tt<="+ty.getTCount()+")\n");
		sb.append("lon1=subwrd(olon,tt)-18\n");
		sb.append("lon2=subwrd(olon,tt)+18\n");
		sb.append("lat1=subwrd(olat,tt)-16\n");
		sb.append("lat2=subwrd(olat,tt)+16\n");
		sb.append("'set lon 'lon1' 'lon2\n");
		sb.append("'set lat 'lat1' 'lat2\n");
		sb.append("'set t 'tt\n\n");
		
		sb.append("'setvpage 2 3 2 1'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'set arrowhead -0.35'\n");
		sb.append("'set arrscl 0.5 40'\n");
		sb.append("'d oacres(u,u.3);oacres(u,v.3)'\n");
		sb.append("'drawtime'\n");
		
		sb.append("'setvpage 2 3 2 2'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'1.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'draw title 'var'1'\n");
		
		sb.append("'setvpage 2 3 2 3'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'2.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'draw title 'var'2'\n");
		
		sb.append("'setvpage 2 3 1 1'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'3.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'draw title 'var'3'\n");
		
		sb.append("'setvpage 2 3 1 2'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'4.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'draw title 'var'4'\n");
		
		sb.append("'setvpage 2 3 1 3'\n");
		sb.append("'setlopts 6 0.18 5 5'\n");
		sb.append("'color -180 180 10 -kind grainbow'\n");
		sb.append("'d oacres(u,'var'5.2*86400)'\n");
		sb.append("'cbarnskip 4'\n");
		sb.append("'draw title 'var'5'\n\n");
		
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		
		sb.append("tt=tt+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("'disable print'\n");
		sb.append("'close 3'\n");
		sb.append("'close 2'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(path+"VortexStat/"+ds+"/"+year+"/"+ty.getID()+"/plot/WNStation.gs")){ fw.write(sb.toString());}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static Variable[] getEddyVariables(Typhoon ty,CsmDescriptor csd){
		/***************** building models ******************/
		SphericalSpatialModel   ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		
		/****************** setting range *******************/
		String tstr=new MDate(ty.getTime(0)).toFormattedDate(
			DateTimeFormatter.ofPattern("yyyy.MM.dd.HH").withLocale(Locale.ENGLISH)
		);
		String tend=new MDate(ty.getTime(ty.getTCount()-1)).toFormattedDate(
			DateTimeFormatter.ofPattern("yyyy.MM.dd.HH").withLocale(Locale.ENGLISH)
		);
		
		Range range1=new Range("time("+tstr+","+tend+")",csd);
		
		/************** initializing variables **************/
		Variable u=new Variable("u",false,range1);
		Variable v=new Variable("v",false,range1);
		Variable w=new Variable("w",false,range1);
		Variable T=new Variable("T",false,range1);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(false);
		cdrs.readData(Type.CUBIC_P,u,v,w,T); cdrs.closeFile();
		
		Variable uo=u.copy();
		Variable vo=v.copy();
		Variable wo=w.copy();
		Variable To=T.copy();
		
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		Variable[] vel=ct.reprojectToCylindrical(u,v);
		u=vel[0];
		v=vel[1];
		
		/************** building application ****************/
		EliassenModelInCC emdl=new EliassenModelInCC(csm); //tp.setPrinting(false);
		ThermoDynamicMethodsInCC tdm=new ThermoDynamicMethodsInCC(csm);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		
		emdl.cStormRelativeAziRadVelocity(u,v);
		emdl.cStormRelativeLatLonVelocity(uo,vo);
		T=tdm.cPotentialTemperature(T);
		
		Variable g=dmdl.cAbsoluteAngularMomentum(u);
		
		u.anomalizeX(); u.setName("ua" ); u.setCommentAndUnit("azimuthal wind anomaly (m s^-1)");
		v.anomalizeX(); v.setName("va" ); v.setCommentAndUnit("radial wind anomaly (m s^-1)");
		T.anomalizeX(); T.setName("Ta" ); T.setCommentAndUnit("temperature anomaly (K)");
		g.anomalizeX(); g.setName("aam"); g.setCommentAndUnit("absolute angular momentum (m^2 s^-1)");
		
		/************** eddy calculation ****************/
		Variable[] divT=emdl.cYZDivergence(v.multiply(T),w.multiply(T));
		Variable[] divG=emdl.cYZDivergence(v.multiply(g),w.multiply(g));
		
		divT[0].setName("htHFC"  ); divT[0].setCommentAndUnit("horizontal flux convergence of eddy heat (K s^-1)");
		divT[1].setName("htVFC"  ); divT[1].setCommentAndUnit("vertical   flux convergence of eddy heat (K s^-1)");
		divG[0].setName("aamHFCR"); divG[0].setCommentAndUnit("horizontal flux convergence of eddy AAM  (K s^-1)");
		divG[1].setName("aamVFCR"); divG[1].setCommentAndUnit("vertical   flux convergence of eddy AAM  (K s^-1)");
		
		dmdl.deWeightBSinEq(divG[0].divideEq(EARTH_RADIUS));
		dmdl.deWeightBSinEq(divG[1].divideEq(EARTH_RADIUS));
		
		return new Variable[]{uo,vo,wo,To,divT[0],divT[1],divG[0],divG[1]};
	}
	
	static Variable[] getVariances(Variable v,Variable[] vs){
		Variable[] re=new Variable[vs.length];
		
		re[0]=StatisticsBasicAnalysisMethods.cVarianceAlongX(v); re[0].setName(v.getName()+"0");
		
		for(int i=1,I=vs.length;i<I;i++)
		re[i]=StatisticsBasicAnalysisMethods.cVarianceAlongX(vs[i]);
		
		return re;
	}
	
	static Variable[] divideComponents(Variable v,int no){
		if(no<0) throw new IllegalArgumentException("no should be positive");
		
		Variable[] vs=new Variable[no+1];
		
		vs[0]=v;
		
		for(int i=1;i<=no;i++){
			vs[i]=getComponent(v,i);
			vs[i].setName(v.getName()+i);
		}
		
		return vs;
	}
	
	static Variable getComponent(Variable v,int Ks){
		if(Ks<0) throw new IllegalArgumentException("Ks should be positive");
		
		Variable vClone=v.copy();
		
		FilterMethods.FFTFilter(vClone,Dimension.X,Ks);
		
		return v.copy().minusEq(vClone);
	}
}
