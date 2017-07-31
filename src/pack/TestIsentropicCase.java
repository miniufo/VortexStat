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
import miniufo.diagnosis.SpatialModel.LevelType;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;


//
public class TestIsentropicCase{
	//
	private static final String name="Haima";
	private static final String path="D:/Data/DiagnosisVortex/"+name+"/Interim/PT/";
	
	private static final Typhoon ty=AccessBestTrack.getTyphoons(
		//"D:/Data/Typhoons/NHC/NHC.txt","name="+name+";time=1Aug1985-30Sep1985;",DataSets.NHC
		"d:/Data/Typhoons/CMA/CMA.txt","name="+name+";time=1sep2004-30Sep2004;",DataSets.CMA
	).get(0);
	
	// modelling
	static void modeling(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			ty.toCSMString(path+name+".ctl",144,91,34,0.15f,5f,265)
		).getDataDescriptor();
		
		//System.out.println(ty.toCSMString(path+"Elena.ctl",144,91,34,0.15f,5f,265));
		
		TicToc.tic("Start calculating "+ty.getName()+" ("+ty.getID()+")");
		
		
		/***************** building models ******************/
		SphericalSpatialModel   ssm=new SphericalSpatialModel(csd.getCtlDescriptor(),LevelType.ISENTROPIC);
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd,LevelType.ISENTROPIC);
		
		Range range=new Range("",csd);
		
		
		/************** initializing variables **************/
		Variable u =new Variable("u"   ,false,range);
		Variable v =new Variable("v"   ,false,range);
		Variable mg=new Variable("mont",false,range);
		Variable pr=new Variable("prs" ,false,range);
		
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(true);
		cdrs.readData(Type.CUBIC_P,u,v,mg,pr); cdrs.closeFile();
		
		
		/************** building application ****************/
		EliassenModelInCC  emdl=new EliassenModelInCC(csm); //tp.setPrinting(false);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable[] vel=ct.reprojectToCylindrical(u,v);	u=vel[0]; v=vel[1];
		emdl.cStormRelativeAziRadVelocity(u,v);
		
		
		/************** variable calculation ****************/
		Variable sig=emdl.cPseudoDensity(pr); Variable sigo=sig.copy();
		Variable sv =sig.multiply(v);
		Variable Mx =emdl.cMontLambda(mg);
		Variable aa =dmdl.cAbsoluteAngularMomentum(u);
		Variable pv =dmdl.cAbsoluteVorticity(u,v).divideEq(sig);
		pv.setName("pv"); pv.setCommentAndUnit("vertical component of potential voticity (m^2 K kg^-1 s^-1)");
		
		
		/************** mean variable calculation ****************/
		Variable  uma= u.copy(); Variable  umm=emdl.weightedAnomalizeX( uma,sig);  uma.setName( "uma");
		Variable  vma= v.copy(); Variable  vmm=emdl.weightedAnomalizeX( vma,sig);  vma.setName( "vma");
		Variable aama=aa.copy(); Variable aamm=emdl.weightedAnomalizeX(aama,sig); aama.setName("aama");
		Variable pvma=pv.copy(); Variable pvmm=emdl.weightedAnomalizeX(pvma,sig); pvma.setName("pvma");
		
		Variable  um=  u.anomalizeX();  um.setName("utm");
		Variable  vm=  v.anomalizeX();  vm.setName("vrm");
		Variable aam= aa.anomalizeX(); aam.setName("aam");
		Variable svm= sv.anomalizeX(); svm.setName("svm");
		Variable pvm= pv.anomalizeX(); pvm.setName("pvm");
		Variable sgm=sig.anomalizeX(); sgm.setName("sgm");
		Variable mgm= mg.anomalizeX(); mgm.setName("mgm");
		Variable prm= pr.anomalizeX(); prm.setName("prm");
		Variable Mxm= Mx.anomalizeX(); Mxm.setName("Mxm");
		
		
		/************** eddy flux calculation ****************/
		Variable PVFlx  =pv.multiply(v).anomalizeX().multiplyEq(sgm);	  PVFlx.setName("pvflx"  );	  PVFlx.setCommentAndUnit("[s][P'v'] (m^2 s^-2)");
		Variable PVmFlx =emdl.cMassWeightedEddyFlux(pvma,vma,sigo);  	 PVmFlx.setName("pvmflx" );  PVmFlx.setCommentAndUnit("[sP*v*]   (m^2 s^-2)");
		Variable sVaMa  =sv.multiply(aa).anomalizeX();					  sVaMa.setName("sVaMa"  );   sVaMa.setCommentAndUnit("[(sv)'M'] (kg m K^-1 s^-2)");
		Variable sVaUa  =sv.multiply(u).anomalizeX();					  sVaUa.setName("sVaUa"  );   sVaUa.setCommentAndUnit("[(sv)'U'] (kg K^-1 s^-2)");
		Variable sVmaMma=emdl.cMassWeightedEddyFlux(vma,aama,sigo);		sVmaMma.setName("sVmaMma"); sVmaMma.setCommentAndUnit("[sv*M*]   (kg m K^-1 s^-2)");
		Variable sVmaUma=emdl.cMassWeightedEddyFlux(vma,uma,sigo);		sVmaUma.setName("sVmaUma"); sVmaUma.setCommentAndUnit("[sv*U*]   (kg K^-1 s^-2)");
		Variable paMxa  =pr.multiply(Mx).anomalizeX();					  paMxa.setName("paMxa"  );	  paMxa.setCommentAndUnit("[p'M¦Ë']   (Pa m^2 s^-2)");
		Variable sWmaMma=new Variable("sWmaMma",sVmaMma);		sWmaMma.setValue(sVmaMma.getUndef());
		Variable sWmaUma=new Variable("sWmaUma",sVmaUma);		sWmaUma.setValue(sVmaUma.getUndef());
		
		Variable[] EPFlxM =emdl.cIsentropicEPVectorByM(sVaMa,paMxa,1f);
		Variable[] EPmFlxM=emdl.cIsentropicFAEPVectorByM(sVmaMma,paMxa,sWmaMma,1f);
		Variable[] EPFlxU =emdl.cIsentropicEPVectorByU(sVaUa,paMxa,1f);
		Variable[] EPmFlxU=emdl.cIsentropicFAEPVectorByU(sVmaUma,paMxa,sWmaUma,1f);
		Variable[] div1=emdl.cYZDivergence( EPFlxU[0], EPFlxU[1]);
		Variable[] div2=emdl.cYZDivergence(EPmFlxU[0],EPmFlxU[1]);
		Variable[] div3=emdl.cYZDivergence( EPFlxM[0], EPFlxM[1]);
		Variable[] div4=emdl.cYZDivergence(EPmFlxM[0],EPmFlxM[1]);
		
		for(int i=0;i<3;i++){
			div1[i].divideEq(sgm); div1[i].setName("EPdivU"+i);		emdl.deWeightBSinEq(div1[i]).divideEq(SpatialModel.EARTH_RADIUS);
			div2[i].divideEq(sgm); div2[i].setName("EPmdivU"+i);	emdl.deWeightBSinEq(div2[i]).divideEq(SpatialModel.EARTH_RADIUS);
			div3[i].divideEq(sgm); div3[i].setName("EPdivM"+i);		emdl.deWeightBSinEq(div3[i]).divideEq(SpatialModel.EARTH_RADIUS);
			div4[i].divideEq(sgm); div4[i].setName("EPmdivM"+i);	emdl.deWeightBSinEq(div4[i]).divideEq(SpatialModel.EARTH_RADIUS);
		}
		
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+name+"IsenTest.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),
			um,vm,aam,pvm,sgm,mgm,prm,Mxm,umm,svm,vmm,aamm,pvmm,PVFlx,PVmFlx,sVaMa,sVaUa,sVmaMma,sVmaUma,paMxa,sWmaMma,sWmaUma,
			div1[0],div1[1],div1[2],div2[0],div2[1],div2[2],div3[0],div3[1],div3[2],div4[0],div4[1],div4[2],
			EPFlxM[0],EPFlxM[1],EPFlxM[2],EPFlxM[3],EPmFlxM[0],EPmFlxM[1],EPmFlxM[2],EPmFlxM[3],EPmFlxM[4],EPmFlxM[5],
			EPFlxU[0],EPFlxU[1],EPFlxU[2],EPFlxU[3],EPmFlxU[0],EPmFlxU[1],EPmFlxU[2],EPmFlxU[3],EPmFlxU[4],EPmFlxU[5]
		);	cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	
	/*** test ***/
	public static void main(String[] args){
		//verticalInterpolation(path+name+"DL.ctl",path+name+".dat",34); System.exit(0);
		
		System.out.println(ty);
		
		modeling(ty);
	}
}
