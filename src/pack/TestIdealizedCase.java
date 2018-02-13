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
import miniufo.application.advanced.EllipticEqSORSolver;
import miniufo.application.advanced.EllipticEqSORSolver.DimCombination;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.basic.ThermoDynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;


//
public class TestIdealizedCase{
	//
	private static final String name="Vicente";
	private static final String path="E:/Data/VortexStat/JMA/2012/1208/";
	
	private static final Typhoon ty=AccessBestTrack.getTyphoons(
		//"d:/Data/Typhoons/NHC/NHC.txt","name="+name+";time=1Aug1985-30Sep1985;",DataSets.NHC
		"d:/Data/Typhoons/JMA/JMA.txt","name="+name+";time=1Jul2012-31Jul2012;",DataSets.JMA
	).get(0);
	
	
	// modelling
	static void modeling(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseFile(path+"1208.csm").getDataDescriptor();
		
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
		EllipticEqSORSolver ees=new EllipticEqSORSolver(csm);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable[] vel=ct.reprojectToCylindrical(u,v);	u=vel[0]; v=vel[1];
		emdl.cStormRelativeAziRadVelocity(u,v);
		
		
		/************** variable calculation ****************/
		Variable g  =dmdl.cAbsoluteAngularMomentum(u);
		Variable th = tdm.cPotentialTemperature(T);
		Variable s  = tdm.cStaticStabilityArgByPT(th);
		
		Variable um =  u.anomalizeX();	 um.setName("utm");
		Variable vm =  v.anomalizeX();	 vm.setName("vrm");
		Variable wm =  w.anomalizeX();	 wm.setName("wm" );
		Variable gm =  g.anomalizeX();	 gm.setName("gm" );
		Variable sm =  s.anomalizeX();	 sm.setName("sm" );
		Variable thm= th.anomalizeX();	thm.setName("thm");
		
		
		/************ calculate stream function *************/
		Variable sf=new Variable("sf",vm); sf.setCommentAndUnit("streamfunction (no)");
		Variable fc=new Variable("fc",vm); fc.setCommentAndUnit("force (no");
		setForcing(fc);
		
		/*** SOR without boundary ***/
		Variable APrime=emdl.cAPrime(sm);  Variable APrimeC=APrime.copy(); APrimeC.setName("Apc");
		Variable BPrime=emdl.cBPrime(thm); Variable BPrimeC=BPrime.copy(); BPrimeC.setName("Bpc");
		Variable CPrime=emdl.cCPrime(gm);  Variable CPrimeC=CPrime.copy(); CPrimeC.setName("Cpc");
		
		ees.setDimCombination(DimCombination.YZ);
		ees.setABC(APrimeC,BPrimeC,CPrimeC);
		ees.setTolerance(1e-10);
		
		ees.solve(sf,fc);
		
		/** calculate radial, vertical velocity components **/
		Variable[] vv=emdl.cVW(sf);
		vv[0].setName("vs");	vv[0].setCommentAndUnit("radial wind (no");
		vv[1].setName("ws");	vv[1].setCommentAndUnit("omega (no");
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+name+"Idealized.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),
			um,vm,wm,thm,gm,sf,fc,APrime,BPrime,CPrime,vv[0],vv[1]
		); cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	static void setForcing(Variable fc){
		int t=fc.getTCount(),z=fc.getZCount(),y=fc.getYCount();
		
		float[][][][] fdata=fc.getData();
		
		double A =-1e-10;
		double j0=35;
		double k0=13;
		double j1=35;
		double k1=65;
		
		if(fc.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				fdata[l][k][j][0] =(float)(A*Math.exp((-Math.pow(j-j0,2)-Math.pow(k-k0,2))/10.0));
				fdata[l][k][j][0]+=(float)(A*Math.exp((-Math.pow(j-j1,2)-Math.pow(k-k1,2))/10.0));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				fdata[k][j][0][l] =(float)(A*Math.exp((-Math.pow(j-j0,2)-Math.pow(k-k0,2))/10.0));
				fdata[k][j][0][l]+=(float)(A*Math.exp((-Math.pow(j-j1,2)-Math.pow(k-k1,2))/10.0));
			}
		}
	}
	
	
	/*** test ***/
	public static void main(String[] args){
		//verticalInterpolation(path+name+"DL.ctl",path+name+".dat",39); System.exit(0);
		
		modeling(ty);
	}
}
