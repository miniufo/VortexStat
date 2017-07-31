//
package pack;

import java.util.concurrent.TimeUnit;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import miniufo.util.TicToc;
import miniufo.diagnosis.Range;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;


//
public class ExactGWBalance{
	//
	private static final String name="Haima";
	private static final String path="D:/Data/DiagnosisVortex/"+name+"/Interim/Prs/";
	
	private static final Typhoon ty=AccessBestTrack.getTyphoons(
		//"d:/Data/Typhoons/NHC/NHC.txt","name="+name+";time=1Aug1985-30Sep1985;",DataSets.NHC
		"d:/Data/Typhoons/CMA/CMA.txt","name="+name+";time=1sep2004-30Sep2004;",DataSets.CMA
	).get(0);
	
	
	// modelling
	static void modelingBalance(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			ty.toCSMString(path+"Haima.ctl",144,91,37,0.15f,-25f,1000)
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
		Variable T =new Variable("t",false,range);
		Variable h =new Variable("h",false,range);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(false);
		cdrs.readData(Type.CUBIC_P,u,v,T,h); cdrs.closeFile();
		
		
		/************** building application ****************/
		EliassenModelInCC  emdl=new EliassenModelInCC(csm); //tp.setPrinting(false);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable[] vel=ct.reprojectToCylindrical(u,v);	u=vel[0]; v=vel[1];
		emdl.cStormRelativeAziRadVelocity(u,v);
		
		Variable umo=u.copy().anomalizeX(); umo.setName("umo");
		
		
		/************** variable calculation ****************/
		Variable um =  u.anomalizeX();
		Variable vm =  v.anomalizeX();
		Variable Tm =  T.anomalizeX();
		Variable hm =  h.anomalizeX(); hm.multiplyEq(9.8f);
		
		Variable am=dmdl.cMeanAbsoluteAngularMomentum(um);
		
		Variable[] blc1=dmdl.cGradientHydrostaticBalancedFieldsByT(am,Tm,hm);
		Variable[] blc2=dmdl.cGradientHydrostaticBalancedFieldsByH(am,hm);
		
		dmdl.cDeviationFromOuterBoundary(blc1[0]);
		dmdl.cDeviationFromOuterBoundary(blc1[1]);
		dmdl.cDeviationFromOuterBoundary(blc2[0]);
		dmdl.cDeviationFromOuterBoundary(blc2[1]);
		
		Variable gdw=dmdl.cMeanGradientWindByCurvedGWB(hm);
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream("d:/BalancedT.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),umo,um,vm,Tm,hm,am,gdw,blc1[0],blc1[1],blc1[2],blc2[0],blc2[1],blc2[2]);
		cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	
	/*** test ***/
	public static void main(String[] args){
		System.out.println(ty);
		
		modelingBalance(ty);
	}
}
