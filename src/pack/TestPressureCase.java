//
package pack;

import java.util.concurrent.TimeUnit;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import miniufo.util.TicToc;
import miniufo.diagnosis.Range;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.basic.ThermoDynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;


//
public class TestPressureCase{
	//
	private static final String name="Haima";
	private static final String path="D:/Data/DiagnosisVortex/"+name+"/Interim/Prs/";
	
	private static final Typhoon ty=AccessBestTrack.getTyphoons(
		//"d:/Data/Typhoons/NHC/NHC.txt","name="+name+";time=1Aug1985-30Sep1985;",DataSets.NHC
		"d:/Data/Typhoons/CMA/CMA.txt","name="+name+";time=1sep2004-30Sep2004;",DataSets.CMA
	).get(0);
	
	
	// modelling
	static void modeling(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			ty.toCSMString(path+name+".ctl",144,91,37,0.15f,-25f,1000)
		).getDataDescriptor();
		
		//System.out.println(ty.toCSMString(path+"Haima.ctl",144,91,37,0.15f,-25f,1000));
		
		TicToc.tic("start calculating "+ty.getName()+" ("+ty.getID()+")");
		
		/***************** building models ******************/
		SphericalSpatialModel   ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		
		Range range=new Range("",csd);
		
		
		/************** initializing variables **************/
		Variable u =new Variable("u",false,range);
		Variable v =new Variable("v",false,range);
		Variable w =new Variable("w",false,range);
		Variable T =new Variable("t",false,range);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd);
		cdrs.readData(Type.CUBIC_P,u,v,w,T); cdrs.closeFile();
		
		
		/************** building application ****************/
		EliassenModelInCC  emdl=new EliassenModelInCC(csm); //tp.setPrinting(false);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		ThermoDynamicMethodsInCC tdm=new ThermoDynamicMethodsInCC(csm);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable[] vel=ct.reprojectToCylindrical(u,v);	u=vel[0]; v=vel[1];
		emdl.cStormRelativeAziRadVelocity(u,v);
		
		
		/************** variable calculation ****************/
		Variable g  =dmdl.cAbsoluteAngularMomentum(u);
		Variable th = tdm.cPotentialTemperature(T);
		
		Variable um =  u.anomalizeX();	 um.setName("utm");
		Variable vm =  v.anomalizeX();	 vm.setName("vrm");
		Variable wm =  w.anomalizeX();	 wm.setName("wm" );
		Variable Tm =  T.anomalizeX();	 Tm.setName("Tm" );
		Variable gm =  g.anomalizeX();	 gm.setName("gm" );
		Variable thm= th.anomalizeX();	thm.setName("thm");
		Variable gm2=dmdl.cMeanAbsoluteAngularMomentum(um); gm2.setName("gm2");
		
		
		/************** eddy flux calculation ****************/
		Variable tava=th.multiply(v); tava.setName("tava"); tava.setCommentAndUnit("eddy heat ¦È'v' (K m s^-1)");
		Variable tawa=th.multiply(w); tawa.setName("tawa"); tawa.setCommentAndUnit("eddy heat ¦È'w' (K Pa s^-1)");
		Variable gava= g.multiply(v); gava.setName("gava"); gava.setCommentAndUnit("eddy AAM  g'v' (m^3 s^-2)");
		Variable gawa= g.multiply(w); gawa.setName("gawa"); gawa.setCommentAndUnit("eddy AAM  g'w' (m^2 Pa s^-2)");
		
		Variable tavam=tava.anomalizeX(); tavam.setName("tavam"); tavam.setCommentAndUnit("eddy heat horizontal flux [¦È'v'] (K m s^-1)");
		Variable tawam=tawa.anomalizeX(); tawam.setName("tawam"); tawam.setCommentAndUnit("eddy heat vertical   flux [¦È'w'] (K Pa s^-1)");
		Variable gavam=gava.anomalizeX(); gavam.setName("gavam"); gavam.setCommentAndUnit("eddy AAM  horizontal flux [g'v'] (m^3 s^-2)");
		Variable gawam=gawa.anomalizeX(); gawam.setName("gawam"); gawam.setCommentAndUnit("eddy AAM  vertical   flux [g'w'] (m^2 Pa s^-2)");
		Variable gavamR=emdl.deWeightBSin(gavam).divideEq(SpatialModel.EARTH_RADIUS);
		Variable gawamR=emdl.deWeightBSin(gawam).divideEq(SpatialModel.EARTH_RADIUS);
		gavamR.setName("gavamR"); gavamR.setCommentAndUnit("eddy AAM  horizontal flux divided by r [g'v']/r (m^2 s^-2)");
		gawamR.setName("gawamR"); gawamR.setCommentAndUnit("eddy AAM  vertical   flux divided by r [g'w']/r (m Pa s^-2)");
		
		Variable sfe   =emdl.cEddyInducedStreamfunction(tavam,thm);
		Variable[] vedy=emdl.cEddyInducedVelocity(sfe);
		Variable[] EPVector=emdl.cEPVector(tavam,gavam,gawam,gm,thm,0.8f);
		Variable EPDiv =dmdl.cYZDivergence(EPVector[0],EPVector[1]); EPDiv.setName("EPDiv");
		Variable EPDivR=dmdl.deWeightBSin(EPDiv).divideEq(SpatialModel.EARTH_RADIUS); EPDivR.setName("EPDivR"); EPDivR.setCommentAndUnit("EPDiv divided by R (m s^-2)");
		
		Variable htHFC =emdl.cEddyHeatHFC(tavam);
		Variable htVFC =emdl.cEddyHeatVFC(tawam);
		Variable aamHFC=emdl.cEddyAAMHFC(gavam);
		Variable aamVFC=emdl.cEddyAAMVFC(gawam);
		Variable aamHFCR=dmdl.deWeightBSin(aamHFC).divideEq(SpatialModel.EARTH_RADIUS); aamHFCR.setName("aamHFCR"); aamHFCR.setCommentAndUnit("aamHFC divided by R (m s^-2)");
		Variable aamVFCR=dmdl.deWeightBSin(aamVFC).divideEq(SpatialModel.EARTH_RADIUS); aamVFCR.setName("aamVFCR"); aamVFCR.setCommentAndUnit("aamVFC divided by R (m s^-2)");
		
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+name+"PrsTest.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),
			um,vm,wm,Tm,thm,gm,gm2,tavam,tawam,gavam,gawam,gavamR,gawamR,htHFC,htVFC,aamHFC,aamVFC,aamHFCR,aamVFCR,
			sfe,vedy[0],vedy[1],EPVector[0],EPVector[1],EPVector[2],EPVector[3],EPDiv,EPDivR
		); cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	
	/*** test ***/
	public static void main(String[] args){
		//verticalInterpolation(path+name+"DL.ctl",path+name+".dat",39); System.exit(0);
		
		System.out.println(ty);
		
		modeling(ty);
	}
}
