//
package pack;

import java.util.List;
import java.util.function.Predicate;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.application.statisticsModel.FilterMethods;
import miniufo.application.statisticsModel.StatisticsBasicAnalysisMethods;
import miniufo.basic.InterpolationModel.Type;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


//
public class WaveNumberIdealized{
	//
	static final Predicate<Typhoon> cond=ty->{
		//int year=new MDate(ty.getTime(0)).getYear();
		//return (year==2013||year==2014||year==2015)&&Integer.parseInt(ty.getID())==1326;
		return "0420".equalsIgnoreCase(ty.getID());
	};
	
	static final DataSets ds=DataSets.JMA;
	static final String path="F:/Data/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons("d:/Data/Typhoons/"+ds+"/"+ds+".txt","",ds);
	
	
	//
	public static void main(String[] args){
		//generateUniformWind(25f); System.exit(0);
		
		all.stream().filter(cond).forEach(ty->{
			processOneTyphoon(ty);
		});
	}
	
	static void generateUniformWind(float amp){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"VortexStat/JMA/2004/0420/0420.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(1,3);lev(200,200)",dd),true,"u","v");
		
		vs[1].setValue(0);
		
		for(int l=0;l<3;l++){
			float[][] data=vs[0].getData()[l][0];
			
			for(int j=0;j<vs[0].getYCount();j++)
			for(int i=0;i<vs[0].getXCount();i++){
				//double tmp=Math.cos(Math.toRadians(dd.getYDef().getSamples()[j]));
				
				//if(tmp>0.1) data[j][i]=(float)(amp/tmp);
				//else data[j][i]=dd.getUndef(null);
				data[j][i]=amp;
			}
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"VortexStat/JMA/Stat/WN/UniformWind.dat");
		dw.writeData(dd,vs); dw.closeFile();
	}
	
	static void processOneTyphoon(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseFile(path+"VortexStat/JMA/Stat/WN/TC.csm").getDataDescriptor();
		
		/***************** building models ******************/
		SphericalSpatialModel   ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		
		Range range1=new Range("t(1,3)",csd);
		
		/************** initializing variables **************/
		Variable u=new Variable("u",false,range1);
		Variable v=new Variable("v",false,range1);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(false);
		cdrs.readData(Type.CUBIC_P,u,v); cdrs.closeFile();
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		Variable[] vel=ct.reprojectToCylindrical(u,v);
		u=vel[0];
		v=vel[1];
		
		/************** building application ****************/
		EliassenModelInCC emdl=new EliassenModelInCC(csm); //tp.setPrinting(false);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		
		Variable g=dmdl.cAbsoluteAngularMomentum(u);
		
		u.anomalizeX(); u.setName("ua" ); u.setCommentAndUnit("azimuthal wind anomaly (m s^-1)");
		v.anomalizeX(); v.setName("va" ); v.setCommentAndUnit("radial wind anomaly (m s^-1)");
		g.anomalizeX(); g.setName("aam"); g.setCommentAndUnit("absolute angular momentum (m^2 s^-1)");
		
		/************** eddy calculation ****************/
		Variable[] divG=emdl.cYZDivergence(v.multiply(g),u.copy().setValue(0));
		
		divG[0].setName("aamHFCR"); divG[0].setCommentAndUnit("horizontal flux convergence of eddy AAM  (K s^-1)");
		dmdl.deWeightBSinEq(divG[0].divideEq(EARTH_RADIUS));
		
		Variable a1=getComponent(divG[0],1); a1.setName("a1");
		Variable a2=getComponent(divG[0],2); a2.setName("a2");
		Variable a3=getComponent(divG[0],3); a3.setName("a3");
		Variable a4=getComponent(divG[0],4); a4.setName("a4");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"VortexStat/JMA/Stat/WN/WN3DIdealized.dat");
		cdws.writeData(csd.getCtlDescriptor(),u,v,g,divG[0],a1,a2,a3,a4); cdws.closeFile();
		
		Variable GVarAll=StatisticsBasicAnalysisMethods.cVarianceAlongX(divG[0]); GVarAll.setName("Gvar");
		Variable G1=StatisticsBasicAnalysisMethods.cVarianceAlongX(a1); G1.setName("G1");
		Variable G2=StatisticsBasicAnalysisMethods.cVarianceAlongX(a2); G2.setName("G2");
		Variable G3=StatisticsBasicAnalysisMethods.cVarianceAlongX(a3); G3.setName("G3");
		Variable G4=StatisticsBasicAnalysisMethods.cVarianceAlongX(a4); G4.setName("G4");
		
		cdws=new CtlDataWriteStream(path+"VortexStat/JMA/Stat/WN/WNVarIdealized.dat");
		cdws.writeData(csd.getCtlDescriptor(),GVarAll,G1,G2,G3,G4);
		cdws.closeFile();
	}
	
	static Variable getComponent(Variable v,int Ks){
		if(Ks<1) throw new IllegalArgumentException("Ks should be larger positive");
		
		Variable vClone=v.copy();
		
		FilterMethods.FFTFilter(vClone,Dimension.X,Ks);
		
		return v.copy().minusEq(vClone);
	}
}
