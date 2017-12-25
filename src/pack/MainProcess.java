//
package pack;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataReadStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import miniufo.test.util.OpenGrADS;
import miniufo.util.DataInterpolation;
import miniufo.util.TicToc;
import miniufo.diagnosis.Range;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.advanced.EllipticEqSORSolver;
import miniufo.application.advanced.EllipticEqSORSolver.DimCombination;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.application.basic.ThermoDynamicMethodsInCC;
import miniufo.application.diagnosticModel.EliassenModelInCC;
import miniufo.basic.InterpolationModel.Type;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


//
public final class MainProcess{
	// modelling
	static void modeling(Typhoon ty,String path,boolean deleteOld){
		if(deleteOld)
		try{
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_model.dat"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_model.ctl"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.dat"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.ctl"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.map"));
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			//ty.toCSMString(path+ty.getID()+".ctl",144,91,77,0.15f,-12.5f,1000)
			ty.interpolateToDT(7200).toCSMString(path+ty.getID()+".ctl",180,91,77,0.1f,-12.5f,1000)
		).getDataDescriptor();
		
		TicToc.tic("Start simulating "+String.format("%10s",ty.getName())+" ("+ty.getID()+")");
		
		
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
		Range range2=new Range("time("+tstr+","+tend+");z(1,1)",csd);
		
		
		/************** initializing variables **************/
		Variable u=new Variable("u",false,range1);
		Variable v=new Variable("v",false,range1);
		Variable w=new Variable("w",false,range1);
		Variable z=new Variable("z",false,range1);
		Variable T=new Variable("t",false,range1);
		Variable q=new Variable("q",false,range1);
		
		Variable u10=new Variable("u10" ,false,range2);
		Variable v10=new Variable("v10" ,false,range2);
		Variable sfp=new Variable("psfc",false,range2);
		Variable sfz=new Variable("zsfc",false,range2);
		Variable blh=new Variable("blh" ,false,range2);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(false);
		cdrs.readData(Type.CUBIC_P,q,T,u,v,w,z);
		cdrs.readData(Type.CUBIC_P,u10,v10,sfp,sfz);
		cdrs.readData(Type.LINEAR,blh);	cdrs.closeFile();
		
		
		/************** building application ****************/
		EliassenModelInCC  emdl=new EliassenModelInCC(csm,sfp);
		ThermoDynamicMethodsInCC  tdmdl=new ThermoDynamicMethodsInCC(csm);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		EllipticEqSORSolver ees=new EllipticEqSORSolver(csm); ees.setPrinting(false);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable uo=u.copy(); // not storm-relative azimuthal-radial wind, for tilting term
		Variable vo=v.copy(); // not storm-relative azimuthal-radial wind, for tilting term
		
		Variable[] vel=null;
		vel=ct.reprojectToCylindrical(uo ,vo );	uo =vel[0];	vo =vel[1];
		vel=ct.reprojectToCylindrical(u  ,v  );	u  =vel[0];	v  =vel[1];
		vel=ct.reprojectToCylindrical(u10,v10);	u10=vel[0];	v10=vel[1];
		
		emdl.cStormRelativeAziRadVelocity(u  ,v  );
		emdl.cStormRelativeAziRadVelocity(u10,v10);
		
		
		/************** variable calculation ****************/
		Variable lh    =tdmdl.cLatentHeating(q,u,v,w,T);
		Variable Tv    =tdmdl.cVirtualTemperature(T,q);
		Variable th    =tdmdl.cPotentialTemperature(T);
		Variable Qth   =tdmdl.cDiabaticHeatingRateByPT(th,u,v,w);
		Variable s     =tdmdl.cStaticStabilityArgByPT(th);
		Variable the   =tdmdl.cEquivalentPotentialTemperature(T,q,tdmdl.cLCLTemperature(T,tdmdl.cRelativeHumidity(T,q)));
		Variable g     = dmdl.cAbsoluteAngularMomentum(u);
		Variable[] tau = dmdl.cSurfaceFrictionalStress((mag)->EliassenModelInCC.Cd,u10,v10);
		Variable[] fsfc= dmdl.cSurfaceFriction(tau[0],tau[1],blh);
		Variable fsfcXL= emdl.assignSurfaceToLevels(fsfc[0],sfz.divide(9.8f),blh,z.divide(9.8f));
		Variable fsfcYL= emdl.assignSurfaceToLevels(fsfc[1],sfz.divide(9.8f),blh,z.divide(9.8f));
		
		Variable sm   =   s.anomalizeX();	sm.setName("sm");
		Variable gm   =   g.anomalizeX();	gm.setName("gm");
		Variable qm   =   q.anomalizeX();	qm.setName("qm");
		Variable lhm  =  lh.anomalizeX();	lhm.setName("lhm");
		Variable Qthm = Qth.anomalizeX();	Qthm.setName("Qthm");
		Variable um   =   u.anomalizeX();	um.setName("um");
		Variable vm   =   v.anomalizeX();	vm.setName("vm");
		Variable wm   =   w.anomalizeX();	wm.setName("wm");
		Variable Tm   =   T.anomalizeX();	Tm.setName("Tm");
		Variable Tvm  =  Tv.anomalizeX();	Tvm.setName("Tvm");
		Variable zm   =   z.anomalizeX();	zm.setName("zm");	zm.setCommentAndUnit("geopotential");
		Variable thm  =  th.anomalizeX();	thm.setName("thm");
		Variable them = the.anomalizeX();	them.setName("them");
		Variable gdm  =dmdl.cMeanGradientWindByCurvedGWB(zm); gdm.setName("gdw");
		Variable fsXLm=fsfcXL.anomalizeX();
		Variable fsYLm=fsfcYL.anomalizeX();
		Variable fvhm =emdl.cHorizontalViscousFriction(um);
		Variable fvvm =emdl.cVerticalViscousFriction(um);
		
		Variable Am=emdl.cA(sm );
		Variable Bm=emdl.cB(thm);
		Variable Cm=emdl.cC(gm );
		Variable BsinB=emdl.weightBSin(Bm);
		Variable CsinB=emdl.weightBSin(Cm);
		BsinB.setName("BsinB"); BsinB.setComment("B * sin(beta) (m^2 kg^-1)");
		CsinB.setName("CsinB"); CsinB.setComment("C * sin(beta) (s^-2)");
		
		
		/************** eddy calculation ****************/
		Variable tava=th.multiply(v); tava.setName("tava"); tava.setCommentAndUnit("eddy heat ¦È'v' (K m s^-1)");
		Variable tawa=th.multiply(w); tawa.setName("tawa"); tawa.setCommentAndUnit("eddy heat ¦È'w' (K Pa s^-1)");
		Variable gava= g.multiply(v); gava.setName("gava"); gava.setCommentAndUnit("eddy AAM  g'v' (m^3 s^-2)");
		Variable gawa= g.multiply(w); gawa.setName("gawa"); gawa.setCommentAndUnit("eddy AAM  g'w' (m^2 Pa s^-2)");
		
		/****
		Variable gu=dmdl.cRelativeAngularMomentum(u);
		
		Variable[] divT=emdl.cYZDivergence(v.multiply(th),w.multiply(th));System.gc();
		Variable[] divG=emdl.cYZDivergence(v.multiply(g ),w.multiply(g ));System.gc();
		Variable[] divU=emdl.cYZDivergence(v.multiply(gu),w.multiply(gu));System.gc();
		
		dmdl.deWeightBSinEq(divG[0].divideEq(EARTH_RADIUS)); dmdl.deWeightBSinEq(divU[0].divideEq(EARTH_RADIUS));
		dmdl.deWeightBSinEq(divG[1].divideEq(EARTH_RADIUS)); dmdl.deWeightBSinEq(divU[1].divideEq(EARTH_RADIUS));
		
		divT[0].setName("htHFC"); divG[0].setName("aamHFCR"); divU[0].setName("ramHFCR");
		divT[1].setName("htVFC"); divG[1].setName("aamVFCR"); divU[1].setName("ramVFCR");
		
		divT[0].setComment("eddy heat horizontal flux convergence (K s^-1)");
		divT[1].setComment("eddy heat vertical   flux convergence (K s^-1)");
		divG[0].setComment("eddy AAM  horizontal flux convergence divided by R (m s^-2)");
		divG[1].setComment("eddy AAM  vertical   flux convergence divided by R (m s^-2)");
		divU[0].setComment("eddy RAM  horizontal flux convergence divided by R (m s^-2)");
		divU[1].setComment("eddy RAM  vertical   flux convergence divided by R (m s^-2)");
		
		CsmDataWriteStream cs=new CsmDataWriteStream(path+ty.getID()+"_stn.dat"); cs.setPrinting(false);
		cs.writeData(csd,stnvar[0],stnvar[1],stnvar[2],stnvar[3],divT[0],divT[1],divG[0],divG[1],divU[0],divU[1]);
		cs.closeFile();*/
		
		
		/**************** eddy flux calculation *****************/
		Variable tavam=tava.anomalizeX(); tavam.setName("tavam"); tavam.setCommentAndUnit("eddy heat horizontal flux [¦È'v'] (K m s^-1)");
		Variable tawam=tawa.anomalizeX(); tawam.setName("tawam"); tawam.setCommentAndUnit("eddy heat vertical   flux [¦È'w'] (K Pa s^-1)");
		Variable gavam=gava.anomalizeX(); gavam.setName("gavam"); gavam.setCommentAndUnit("eddy AAM  horizontal flux [g'v'] (m^3 s^-2)");
		Variable gawam=gawa.anomalizeX(); gawam.setName("gawam"); gawam.setCommentAndUnit("eddy AAM  vertical   flux [g'w'] (m^2 Pa s^-2)");
		Variable gavaR=dmdl.deWeightBSin(gavam).divideEq(EARTH_RADIUS); gavaR.setName("gavaR"); gavaR.setCommentAndUnit("eddy AAM horizontal flux divided by R [g'v']/r (m^2 s^-2)");
		Variable gawaR=dmdl.deWeightBSin(gawam).divideEq(EARTH_RADIUS); gawaR.setName("gawaR"); gawaR.setCommentAndUnit("eddy AAM  vertical  flux divided by R [g'w']/r (m Pa s^-2)");
		
		Variable sfe    =emdl.cEddyInducedStreamfunction(tavam,thm);
		Variable sfe2   =emdl.cCoordinateIndependentEddyInducedStreamfunction(tavam,tawam,thm); sfe2.setName("sfe2");
		Variable[] vedy =emdl.cEddyInducedVelocity(sfe);
		Variable[] vedy2=emdl.cEddyInducedVelocity(sfe2); vedy2[0].setName("vedy2"); vedy2[1].setName("wedy2");
		Variable[] EPVector=emdl.cEPVector(tavam,gavam,gawam,gm,thm,0.8f);
		Variable EPDiv  =dmdl.cYZDivergence(EPVector[0],EPVector[1]); EPDiv.setName("EPDiv");
		Variable EPDivR =dmdl.deWeightBSin(EPDiv).divideEq(EARTH_RADIUS); EPDivR.setName("EPDivR"); EPDivR.setCommentAndUnit("EPDiv divided by R (m s^-2)");
		
		
		/**************** force calculation *****************/
		Variable f01=emdl.cDiabaticHeatingForce(Qthm);	f01.setCommentAndUnit("force due to latent heating (m^2 kg^-1 s^-1)");
		Variable f04=emdl.cEddyHeatHFCForce(tavam);
		Variable f05=emdl.cEddyHeatVFCForce(tawam);
		Variable f06=emdl.cEddyAAMHFCForce(gm,gavam);
		Variable f07=emdl.cEddyAAMVFCForce(gm,gawam);
		Variable f08=emdl.cFrictionalTorqueForce(gm,fsXLm); f08.setName("fsfFor"); f08.setCommentAndUnit("force due to surface frictional torque (m^2 kg^-1 s^-1)");
		Variable f09=emdl.cFrictionalTorqueForce(gm,fvhm ); f09.setName("fvhFor"); f09.setCommentAndUnit("force due to horizontal viscous torque (m^2 kg^-1 s^-1)");
		Variable f10=emdl.cFrictionalTorqueForce(gm,fvvm ); f10.setName("fvvFor"); f10.setCommentAndUnit("force due to vertical viscous torque (m^2 kg^-1 s^-1)");
		Variable f11=emdl.cTiltingForce(gm,uo,vo);
		
		Variable dhr   =emdl.cDiabaticHeatingRate(Qthm);
		Variable htHFC =emdl.cEddyHeatHFC(tavam);
		Variable htVFC =emdl.cEddyHeatVFC(tawam);
		Variable aamHFC=emdl.cEddyAAMHFC(gavam);
		Variable aamVFC=emdl.cEddyAAMVFC(gawam);
		Variable aamHFCR=dmdl.deWeightBSin(aamHFC).divideEq(EARTH_RADIUS); aamHFCR.setName("aamHFCR"); aamHFCR.setCommentAndUnit("aamHFC divided by R (m s^-2)");
		Variable aamVFCR=dmdl.deWeightBSin(aamVFC).divideEq(EARTH_RADIUS); aamVFCR.setName("aamVFCR"); aamVFCR.setCommentAndUnit("aamVFC divided by R (m s^-2)");
		Variable frictq=emdl.cFrictionalTorque(fsXLm); frictq.setName("frictq"); frictq.setCommentAndUnit("surface frictional torque (m^2 s^-2)");
		Variable fvhmtq=emdl.cFrictionalTorque(fvhm);  fvhmtq.setName("fvhmtq"); fvhmtq.setCommentAndUnit("horizontal viscous frictional torque (m^2 s^-2)");
		Variable fvvmtq=emdl.cFrictionalTorque(fvvm);  fvvmtq.setName("fvvmtq"); fvvmtq.setCommentAndUnit("vertical viscous frictional torque (m^2 s^-2)");
		Variable tilt  =emdl.cTilting(uo,vo);
		
		Variable f01f=emdl.cHeatFF(dhr);
		Variable f04f=emdl.cHeatFF(htHFC);
		Variable f05f=emdl.cHeatFF(htVFC);
		Variable f06f=emdl.cMomentumFF(gm,aamHFC);
		Variable f07f=emdl.cMomentumFF(gm,aamVFC);
		Variable f08f=emdl.cMomentumFF(gm,frictq); f08f.setName("fsfFF"); f08f.setCommentAndUnit("forcing factor due to surface frictional torque (m s^-3)");
		Variable f09f=emdl.cMomentumFF(gm,fvhmtq); f09f.setName("fhvFF"); f09f.setCommentAndUnit("forcing factor due to horizontal viscous torque (m s^-3)");
		Variable f10f=emdl.cMomentumFF(gm,fvvmtq); f10f.setName("fvvFF"); f10f.setCommentAndUnit("forcing factor due to vertical viscous torque (m s^-3)");
		Variable f11f=emdl.cMomentumFF(gm,tilt);
		
		
		/********** calculate thermodynamic force ***********/
		Variable fat=f01.copy();	fat.setName("fat"); fat.setCommentAndUnit("all thermal dynamic force (m^2 kg^-1 s^-1)");
		fat.plusEq(f04); fat.plusEq(f05);
		
		
		/************* calculate dynamic force **************/
		Variable fad=f06.copy();	fad.setName("fad"); fad.setCommentAndUnit("all dynamic force (m^2 kg^-1 s^-1)");
		fad.plusEq(f07); fad.plusEq(f08); fad.plusEq(f09); fad.plusEq(f10); fad.plusEq(f11);
		
		
		/************** calculate total force ***************/
		Variable faf=f01.copy();	faf.setName("faf"); faf.setCommentAndUnit("total force (m^2 kg^-1 s^-1)");
		faf.plusEq(f04); faf.plusEq(f05);
		faf.plusEq(f06); faf.plusEq(f07); faf.plusEq(f08); faf.plusEq(f09); faf.plusEq(f10); faf.plusEq(f11);
		
		
		/************ calculate stream function *************/
		Variable sffr=new Variable("sffr",vm); sffr.setCommentAndUnit("streamfunction forced by all internal forcings (Pa m s^-1)");
		Variable sf01=new Variable("sf01",vm); sf01.setCommentAndUnit("streamfunction forced by latent heating (Pa m s^-1)");
		Variable sf04=new Variable("sf04",vm); sf04.setCommentAndUnit("streamfunction forced by eddy heat HFC (Pa m s^-1)");
		Variable sf05=new Variable("sf05",vm); sf05.setCommentAndUnit("streamfunction forced by eddy heat VFC (Pa m s^-1)");
		Variable sf06=new Variable("sf06",vm); sf06.setCommentAndUnit("streamfunction forced by eddy AAM  HFC (Pa m s^-1)");
		Variable sf07=new Variable("sf07",vm); sf07.setCommentAndUnit("streamfunction forced by eddy AAM  VFC (Pa m s^-1)");
		Variable sf08=new Variable("sf08",vm); sf08.setCommentAndUnit("streamfunction forced by surface frictional torque (Pa m s^-1)");
		Variable sf09=new Variable("sf09",vm); sf09.setCommentAndUnit("streamfunction forced by horizontal viscous torque (Pa m s^-1)");
		Variable sf10=new Variable("sf10",vm); sf10.setCommentAndUnit("streamfunction forced by vertical viscous torque (Pa m s^-1)");
		Variable sf11=new Variable("sf11",vm); sf11.setCommentAndUnit("streamfunction forced by tilting effect (Pa m s^-1)");
		Variable sftf=new Variable("sftf",vm); sftf.setCommentAndUnit("streamfunction forced by all thermal forces (Pa m s^-1)");
		Variable sfdf=new Variable("sfdf",vm); sfdf.setCommentAndUnit("streamfunction forced by all dynamic forces (Pa m s^-1)");
		Variable sfbd=new Variable("sfbd",vm); sfbd.setCommentAndUnit("streamfunction forced by boundary effect (Pa m s^-1)");
		Variable sfsm=new Variable("sfsm",vm); sfsm.setCommentAndUnit("streamfunction forced by all internal forcings and boundary effect (Pa m s^-1)");
		
		/*** SOR without boundary ***/
		Variable APrime=emdl.cAPrime(sm);  Variable APrimeC=APrime.copy(); APrimeC.setName("Apc");
		Variable BPrime=emdl.cBPrime(thm); Variable BPrimeC=BPrime.copy(); BPrimeC.setName("Bpc");
		Variable CPrime=emdl.cCPrime(gm);  Variable CPrimeC=CPrime.copy(); CPrimeC.setName("Cpc");
		
		ees.setDimCombination(DimCombination.YZ);
		ees.setABC(APrimeC,BPrimeC,CPrimeC);
		ees.setTolerance(1e-10);
		
		ees.solve(sf01,f01,true);
		ees.solve(sf04,f04);
		ees.solve(sf05,f05);
		ees.solve(sf06,f06);
		ees.solve(sf07,f07);
		ees.solve(sf08,f08);
		ees.solve(sf09,f09);
		ees.solve(sf10,f10);
		ees.solve(sf11,f11);
		ees.solve(sftf,fat);
		ees.solve(sfdf,fad);
		ees.solve(sffr,faf);

		/*** SOR with boundary ***/
		emdl.initialSFBoundary(sfbd,vm,wm);
		emdl.initialSFBoundary(sfsm,vm,wm);
		ees.solve(sfbd,null);
		ees.solve(sfsm,faf);
		
		
		/** calculate radial, vertical velocity components **/
		Variable[] vvfr=emdl.cVW(sffr);
		vvfr[0].setName("vsfr");	vvfr[0].setCommentAndUnit("radial wind forced by all internal forcings (m s^-1)");
		vvfr[1].setName("wsfr");	vvfr[1].setCommentAndUnit("omega forced by all internal forcings (Pa s^-1)");
		Variable[] vv01=emdl.cVW(sf01);
		vv01[0].setName("vs01");	vv01[0].setCommentAndUnit("radial wind forced by latent heating (m s^-1)");
		vv01[1].setName("ws01");	vv01[1].setCommentAndUnit("omega forced by latent heating (Pa s^-1)");
		Variable[] vv04=emdl.cVW(sf04);
		vv04[0].setName("vs04");	vv04[0].setCommentAndUnit("radial wind forced by eddy heat HFC (m s^-1)");
		vv04[1].setName("ws04");	vv04[1].setCommentAndUnit("omega forced by eddy heat HFC (Pa s^-1)");
		Variable[] vv05=emdl.cVW(sf05);
		vv05[0].setName("vs05");	vv05[0].setCommentAndUnit("radial wind forced by eddy heat VFC (m s^-1)");
		vv05[1].setName("ws05");	vv05[1].setCommentAndUnit("omega forced by eddy heat VFC (Pa s^-1)");
		Variable[] vv06=emdl.cVW(sf06);
		vv06[0].setName("vs06");	vv06[0].setCommentAndUnit("radial wind forced by eddy AAM HFC (m s^-1)");
		vv06[1].setName("ws06");	vv06[1].setCommentAndUnit("omega forced by eddy AAM HFC (Pa s^-1)");
		Variable[] vv07=emdl.cVW(sf07);
		vv07[0].setName("vs07");	vv07[0].setCommentAndUnit("radial wind forced by eddy AAM VFC (m s^-1)");
		vv07[1].setName("ws07");	vv07[1].setCommentAndUnit("omega forced by eddy AAM VFC (Pa s^-1)");
		Variable[] vv08=emdl.cVW(sf08);
		vv08[0].setName("vs08");	vv08[0].setCommentAndUnit("radial wind forced by surface frictional torque (m s^-1)");
		vv08[1].setName("ws08");	vv08[1].setCommentAndUnit("omega forced by surface frictional torque (Pa s^-1)");
		Variable[] vv09=emdl.cVW(sf09);
		vv09[0].setName("vs09");	vv09[0].setCommentAndUnit("radial wind forced by horizontal viscous torque (m s^-1)");
		vv09[1].setName("ws09");	vv09[1].setCommentAndUnit("omega forced by horizontal viscous torque (Pa s^-1)");
		Variable[] vv10=emdl.cVW(sf10);
		vv10[0].setName("vs10");	vv10[0].setCommentAndUnit("radial wind forced by vertical viscous torque (m s^-1)");
		vv10[1].setName("ws10");	vv10[1].setCommentAndUnit("omega forced by vertical viscous torque (Pa s^-1)");
		Variable[] vv11=emdl.cVW(sf11);
		vv11[0].setName("vs11");	vv11[0].setCommentAndUnit("radial wind forced by tilting effect (m s^-1)");
		vv11[1].setName("ws11");	vv11[1].setCommentAndUnit("omega forced by tilting effect (Pa s^-1)");
		Variable[] vvtf=emdl.cVW(sftf);
		vvtf[0].setName("vstf");	vvtf[0].setCommentAndUnit("radial wind forced by all thermal forces (m s^-1)");
		vvtf[1].setName("wstf");	vvtf[1].setCommentAndUnit("omega forced by all thermal forces (Pa s^-1)");
		Variable[] vvdf=emdl.cVW(sfdf);
		vvdf[0].setName("vsdf");	vvdf[0].setCommentAndUnit("radial wind forced by all dynamic forces (m s^-1)");
		vvdf[1].setName("wsdf");	vvdf[1].setCommentAndUnit("omega forced by all dynamic forces (Pa s^-1)");
		Variable[] vvbd=emdl.cVW(sfbd);
		vvbd[0].setName("vsbd");	vvbd[0].setCommentAndUnit("radial wind forced by boundary effect (m s^-1)");
		vvbd[1].setName("wsbd");	vvbd[1].setCommentAndUnit("omega forced by boundary effect (Pa s^-1)");
		Variable[] vvsm=emdl.cVW(sfsm );
		vvsm[0].setName("vsm" );	vvsm[0].setCommentAndUnit("radial wind forced by all (m s^-1)");
		vvsm[1].setName("wsm" );	vvsm[1].setCommentAndUnit("omega forced by all (Pa s^-1)");
		
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+ty.getID()+"_model.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),
			  sm, gm ,gdm  , qm  , lhm ,Qthm , um , vm , wm , Am , Bm , Cm ,BsinB,CsinB, Tm ,Tvm ,thm ,them, zm ,fsXLm,fsYLm,APrime,BPrime,CPrime,
			fvhm,fvvm,tavam,tawam,gavam,gawam,gavaR,gawaR,dhr ,htHFC,htVFC,aamHFC,aamVFC,aamHFCR,aamVFCR,frictq,fvhmtq,fvvmtq,tilt,APrimeC,BPrimeC,CPrimeC,
			f01f,f04f,f05f,f06f,f07f,f08f,f09f,f10f,f11f,sfe,sfe2,vedy[0],vedy[1],vedy2[0],vedy2[1],EPVector[0],EPVector[1],EPVector[2],EPVector[3],EPDiv,EPDivR,
			 faf, fat, fad, f01, f04, f05, f06, f07, f08, f09, f10, f11,
			sffr,sftf,sfdf,sf01,sf04,sf05,sf06,sf07,sf08,sf09,sf10,sf11,sfbd,sfsm,
			vvfr[0],vv01[0],vv04[0],vv05[0],vv06[0],vv07[0],vv08[0],vv09[0],vv10[0],vv11[0],vvtf[0],vvdf[0],vvbd[0],vvsm[0],
			vvfr[1],vv01[1],vv04[1],vv05[1],vv06[1],vv07[1],vv08[1],vv09[1],vv10[1],vv11[1],vvtf[1],vvdf[1],vvbd[1],vvsm[1]
		);
		cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	static void modelingWithoutBLH(Typhoon ty,String path,boolean deleteOld){
		if(deleteOld)
		try{
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_model.dat"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_model.ctl"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.dat"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.ctl"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+"_stn.map"));
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			//ty.toCSMString(path+ty.getID()+".ctl",144,91,77,0.15f,-12.5f,1000)
			//ty.toCSMString(path+ty.getID()+".ctl",150,161,77,0.05f,-12.5f,1000)
			ty.interpolateToDT(7200).toCSMString(path+ty.getID()+".ctl",180,91,77,0.1f,-12.5f,1000)
		).getDataDescriptor();
		
		TicToc.tic("Start simulating "+String.format("%10s",ty.getName())+" ("+ty.getID()+")");
		
		
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
		Range range2=new Range("time("+tstr+","+tend+");z(1,1)",csd);
		
		
		/************** initializing variables **************/
		Variable u=new Variable("u",false,range1);
		Variable v=new Variable("v",false,range1);
		Variable w=new Variable("w",false,range1);
		Variable z=new Variable("z",false,range1);
		Variable T=new Variable("t",false,range1);
		Variable q=new Variable("q",false,range1);
		
		Variable u10=new Variable("u10" ,false,range2);
		Variable v10=new Variable("v10" ,false,range2);
		
		/****************** reading data ********************/
		CsmDataReadStream  cdrs=new CsmDataReadStream(csd); cdrs.setPrinting(false);
		cdrs.readData(Type.CUBIC_P,q,T,u,v,w,z);
		cdrs.readData(Type.CUBIC_P,u10,v10);
		cdrs.closeFile();
		
		
		/************** building application ****************/
		EliassenModelInCC  emdl=new EliassenModelInCC(csm);
		ThermoDynamicMethodsInCC  tdmdl=new ThermoDynamicMethodsInCC(csm);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		EllipticEqSORSolver ees=new EllipticEqSORSolver(csm); ees.setPrinting(false);
		
		
		/************** vectors reprojection ****************/
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		Variable uo=u.copy(); // not storm-relative azimuthal-radial wind, for tilting term
		Variable vo=v.copy(); // not storm-relative azimuthal-radial wind, for tilting term
		
		Variable[] vel=null;
		vel=ct.reprojectToCylindrical(uo ,vo );	uo =vel[0];	vo =vel[1];
		vel=ct.reprojectToCylindrical(u  ,v  );	u  =vel[0];	v  =vel[1];
		vel=ct.reprojectToCylindrical(u10,v10);	u10=vel[0];	v10=vel[1];
		
		emdl.cStormRelativeAziRadVelocity(u  ,v  );
		emdl.cStormRelativeAziRadVelocity(u10,v10);
		
		
		/************** variable calculation ****************/
		Variable lh    =tdmdl.cLatentHeating(q,u,v,w,T);
		Variable Tv    =tdmdl.cVirtualTemperature(T,q);
		Variable th    =tdmdl.cPotentialTemperature(T);
		Variable Qth   =tdmdl.cDiabaticHeatingRateByPT(th,u,v,w);
		Variable s     =tdmdl.cStaticStabilityArgByPT(th);
		Variable the   =tdmdl.cEquivalentPotentialTemperature(T,q,tdmdl.cLCLTemperature(T,tdmdl.cRelativeHumidity(T,q)));
		Variable g     = dmdl.cAbsoluteAngularMomentum(u);
		
		Variable sm   =   s.anomalizeX();	sm.setName("sm");
		Variable gm   =   g.anomalizeX();	gm.setName("gm");
		Variable qm   =   q.anomalizeX();	qm.setName("qm");
		Variable lhm  =  lh.anomalizeX();	lhm.setName("lhm");
		Variable Qthm = Qth.anomalizeX();	Qthm.setName("Qthm");
		Variable um   =   u.anomalizeX();	um.setName("um");
		Variable vm   =   v.anomalizeX();	vm.setName("vm");
		Variable wm   =   w.anomalizeX();	wm.setName("wm");
		Variable Tm   =   T.anomalizeX();	Tm.setName("Tm");
		Variable Tvm  =  Tv.anomalizeX();	Tvm.setName("Tvm");
		Variable zm   =   z.anomalizeX();	zm.setName("zm");	zm.setCommentAndUnit("geopotential");
		Variable thm  =  th.anomalizeX();	thm.setName("thm");
		Variable them = the.anomalizeX();	them.setName("them");
		Variable gdm  =dmdl.cMeanGradientWindByCurvedGWB(zm); gdm.setName("gdw");
		Variable fvhm =emdl.cHorizontalViscousFriction(um);
		Variable fvvm =emdl.cVerticalViscousFriction(um);
		
		Variable Am=emdl.cA(sm );
		Variable Bm=emdl.cB(thm);
		Variable Cm=emdl.cC(gm );
		Variable BsinB=emdl.weightBSin(Bm);
		Variable CsinB=emdl.weightBSin(Cm);
		BsinB.setName("BsinB"); BsinB.setComment("B * sin(beta) (m^2 kg^-1)");
		CsinB.setName("CsinB"); CsinB.setComment("C * sin(beta) (s^-2)");
		
		
		/************** eddy calculation ****************/
		Variable tava=th.multiply(v); tava.setName("tava"); tava.setCommentAndUnit("eddy heat ¦È'v' (K m s^-1)");
		Variable tawa=th.multiply(w); tawa.setName("tawa"); tawa.setCommentAndUnit("eddy heat ¦È'w' (K Pa s^-1)");
		Variable gava= g.multiply(v); gava.setName("gava"); gava.setCommentAndUnit("eddy AAM  g'v' (m^3 s^-2)");
		Variable gawa= g.multiply(w); gawa.setName("gawa"); gawa.setCommentAndUnit("eddy AAM  g'w' (m^2 Pa s^-2)");
		
		/****
		Variable gu=dmdl.cRelativeAngularMomentum(u);
		
		Variable[] divT=emdl.cYZDivergence(v.multiply(th),w.multiply(th));System.gc();
		Variable[] divG=emdl.cYZDivergence(v.multiply(g ),w.multiply(g ));System.gc();
		Variable[] divU=emdl.cYZDivergence(v.multiply(gu),w.multiply(gu));System.gc();
		
		dmdl.deWeightBSinEq(divG[0].divideEq(EARTH_RADIUS)); dmdl.deWeightBSinEq(divU[0].divideEq(EARTH_RADIUS));
		dmdl.deWeightBSinEq(divG[1].divideEq(EARTH_RADIUS)); dmdl.deWeightBSinEq(divU[1].divideEq(EARTH_RADIUS));
		
		divT[0].setName("htHFC"); divG[0].setName("aamHFCR"); divU[0].setName("ramHFCR");
		divT[1].setName("htVFC"); divG[1].setName("aamVFCR"); divU[1].setName("ramVFCR");
		
		divT[0].setComment("eddy heat horizontal flux convergence (K s^-1)");
		divT[1].setComment("eddy heat vertical   flux convergence (K s^-1)");
		divG[0].setComment("eddy AAM  horizontal flux convergence divided by R (m s^-2)");
		divG[1].setComment("eddy AAM  vertical   flux convergence divided by R (m s^-2)");
		divU[0].setComment("eddy RAM  horizontal flux convergence divided by R (m s^-2)");
		divU[1].setComment("eddy RAM  vertical   flux convergence divided by R (m s^-2)");
		
		CsmDataWriteStream cs=new CsmDataWriteStream(path+ty.getID()+"_stn.dat"); cs.setPrinting(false);
		cs.writeData(csd,stnvar[0],stnvar[1],stnvar[2],stnvar[3],divT[0],divT[1],divG[0],divG[1],divU[0],divU[1]);
		cs.closeFile();*/
		
		
		/**************** eddy flux calculation *****************/
		Variable tavam=tava.anomalizeX(); tavam.setName("tavam"); tavam.setCommentAndUnit("eddy heat horizontal flux [¦È'v'] (K m s^-1)");
		Variable tawam=tawa.anomalizeX(); tawam.setName("tawam"); tawam.setCommentAndUnit("eddy heat vertical   flux [¦È'w'] (K Pa s^-1)");
		Variable gavam=gava.anomalizeX(); gavam.setName("gavam"); gavam.setCommentAndUnit("eddy AAM  horizontal flux [g'v'] (m^3 s^-2)");
		Variable gawam=gawa.anomalizeX(); gawam.setName("gawam"); gawam.setCommentAndUnit("eddy AAM  vertical   flux [g'w'] (m^2 Pa s^-2)");
		Variable gavaR=dmdl.deWeightBSin(gavam).divideEq(EARTH_RADIUS); gavaR.setName("gavaR"); gavaR.setCommentAndUnit("eddy AAM horizontal flux divided by R [g'v']/r (m^2 s^-2)");
		Variable gawaR=dmdl.deWeightBSin(gawam).divideEq(EARTH_RADIUS); gawaR.setName("gawaR"); gawaR.setCommentAndUnit("eddy AAM  vertical  flux divided by R [g'w']/r (m Pa s^-2)");
		
		Variable sfe    =emdl.cEddyInducedStreamfunction(tavam,thm);
		Variable sfe2   =emdl.cCoordinateIndependentEddyInducedStreamfunction(tavam,tawam,thm); sfe2.setName("sfe2");
		Variable[] vedy =emdl.cEddyInducedVelocity(sfe);
		Variable[] vedy2=emdl.cEddyInducedVelocity(sfe2); vedy2[0].setName("vedy2"); vedy2[1].setName("wedy2");
		Variable[] EPVector=emdl.cEPVector(tavam,gavam,gawam,gm,thm,0.8f);
		Variable EPDiv  =dmdl.cYZDivergence(EPVector[0],EPVector[1]); EPDiv.setName("EPDiv");
		Variable EPDivR =dmdl.deWeightBSin(EPDiv).divideEq(EARTH_RADIUS); EPDivR.setName("EPDivR"); EPDivR.setCommentAndUnit("EPDiv divided by R (m s^-2)");
		
		
		/**************** force calculation *****************/
		Variable f01=emdl.cDiabaticHeatingForce(Qthm);	f01.setCommentAndUnit("force due to latent heating (m^2 kg^-1 s^-1)");
		Variable f04=emdl.cEddyHeatHFCForce(tavam);
		Variable f05=emdl.cEddyHeatVFCForce(tawam);
		Variable f06=emdl.cEddyAAMHFCForce(gm,gavam);
		Variable f07=emdl.cEddyAAMVFCForce(gm,gawam);
		Variable f09=emdl.cFrictionalTorqueForce(gm,fvhm ); f09.setName("fvhFor"); f09.setCommentAndUnit("force due to horizontal viscous torque (m^2 kg^-1 s^-1)");
		Variable f10=emdl.cFrictionalTorqueForce(gm,fvvm ); f10.setName("fvvFor"); f10.setCommentAndUnit("force due to vertical viscous torque (m^2 kg^-1 s^-1)");
		Variable f11=emdl.cTiltingForce(gm,uo,vo);
		
		Variable dhr   =emdl.cDiabaticHeatingRate(Qthm);
		Variable htHFC =emdl.cEddyHeatHFC(tavam);
		Variable htVFC =emdl.cEddyHeatVFC(tawam);
		Variable aamHFC=emdl.cEddyAAMHFC(gavam);
		Variable aamVFC=emdl.cEddyAAMVFC(gawam);
		Variable aamHFCR=dmdl.deWeightBSin(aamHFC).divideEq(EARTH_RADIUS); aamHFCR.setName("aamHFCR"); aamHFCR.setCommentAndUnit("aamHFC divided by R (m s^-2)");
		Variable aamVFCR=dmdl.deWeightBSin(aamVFC).divideEq(EARTH_RADIUS); aamVFCR.setName("aamVFCR"); aamVFCR.setCommentAndUnit("aamVFC divided by R (m s^-2)");
		Variable fvhmtq=emdl.cFrictionalTorque(fvhm);  fvhmtq.setName("fvhmtq"); fvhmtq.setCommentAndUnit("horizontal viscous frictional torque (m^2 s^-2)");
		Variable fvvmtq=emdl.cFrictionalTorque(fvvm);  fvvmtq.setName("fvvmtq"); fvvmtq.setCommentAndUnit("vertical viscous frictional torque (m^2 s^-2)");
		Variable tilt  =emdl.cTilting(uo,vo);
		
		Variable f01f=emdl.cHeatFF(dhr);
		Variable f04f=emdl.cHeatFF(htHFC);
		Variable f05f=emdl.cHeatFF(htVFC);
		Variable f06f=emdl.cMomentumFF(gm,aamHFC);
		Variable f07f=emdl.cMomentumFF(gm,aamVFC);
		Variable f09f=emdl.cMomentumFF(gm,fvhmtq); f09f.setName("fhvFF"); f09f.setCommentAndUnit("forcing factor due to horizontal viscous torque (m s^-3)");
		Variable f10f=emdl.cMomentumFF(gm,fvvmtq); f10f.setName("fvvFF"); f10f.setCommentAndUnit("forcing factor due to vertical viscous torque (m s^-3)");
		Variable f11f=emdl.cMomentumFF(gm,tilt);
		
		
		/********** calculate thermodynamic force ***********/
		Variable fat=f01.copy();	fat.setName("fat"); fat.setCommentAndUnit("all thermal dynamic force (m^2 kg^-1 s^-1)");
		fat.plusEq(f04); fat.plusEq(f05);
		
		
		/************* calculate dynamic force **************/
		Variable fad=f06.copy();	fad.setName("fad"); fad.setCommentAndUnit("all dynamic force (m^2 kg^-1 s^-1)");
		fad.plusEq(f07); fad.plusEq(f09); fad.plusEq(f10); fad.plusEq(f11);
		
		
		/************** calculate total force ***************/
		Variable faf=f01.copy();	faf.setName("faf"); faf.setCommentAndUnit("total force (m^2 kg^-1 s^-1)");
		faf.plusEq(f04); faf.plusEq(f05);
		faf.plusEq(f06); faf.plusEq(f07); faf.plusEq(f09); faf.plusEq(f10); faf.plusEq(f11);
		
		
		/************ calculate stream function *************/
		Variable sffr=new Variable("sffr",vm); sffr.setCommentAndUnit("streamfunction forced by all internal forcings (Pa m s^-1)");
		Variable sf01=new Variable("sf01",vm); sf01.setCommentAndUnit("streamfunction forced by latent heating (Pa m s^-1)");
		Variable sf04=new Variable("sf04",vm); sf04.setCommentAndUnit("streamfunction forced by eddy heat HFC (Pa m s^-1)");
		Variable sf05=new Variable("sf05",vm); sf05.setCommentAndUnit("streamfunction forced by eddy heat VFC (Pa m s^-1)");
		Variable sf06=new Variable("sf06",vm); sf06.setCommentAndUnit("streamfunction forced by eddy AAM  HFC (Pa m s^-1)");
		Variable sf07=new Variable("sf07",vm); sf07.setCommentAndUnit("streamfunction forced by eddy AAM  VFC (Pa m s^-1)");
		Variable sf09=new Variable("sf09",vm); sf09.setCommentAndUnit("streamfunction forced by horizontal viscous torque (Pa m s^-1)");
		Variable sf10=new Variable("sf10",vm); sf10.setCommentAndUnit("streamfunction forced by vertical viscous torque (Pa m s^-1)");
		Variable sf11=new Variable("sf11",vm); sf11.setCommentAndUnit("streamfunction forced by tilting effect (Pa m s^-1)");
		Variable sftf=new Variable("sftf",vm); sftf.setCommentAndUnit("streamfunction forced by all thermal forces (Pa m s^-1)");
		Variable sfdf=new Variable("sfdf",vm); sfdf.setCommentAndUnit("streamfunction forced by all dynamic forces (Pa m s^-1)");
		Variable sfbd=new Variable("sfbd",vm); sfbd.setCommentAndUnit("streamfunction forced by boundary effect (Pa m s^-1)");
		Variable sfsm=new Variable("sfsm",vm); sfsm.setCommentAndUnit("streamfunction forced by all internal forcings and boundary effect (Pa m s^-1)");
		
		/*** SOR without boundary ***/
		Variable APrime=emdl.cAPrime(sm);  Variable APrimeC=APrime.copy(); APrimeC.setName("Apc");
		Variable BPrime=emdl.cBPrime(thm); Variable BPrimeC=BPrime.copy(); BPrimeC.setName("Bpc");
		Variable CPrime=emdl.cCPrime(gm);  Variable CPrimeC=CPrime.copy(); CPrimeC.setName("Cpc");
		
		ees.setDimCombination(DimCombination.YZ);
		ees.setABC(APrimeC,BPrimeC,CPrimeC);
		ees.setTolerance(1e-10);
		
		ees.solve(sf01,f01,true);
		ees.solve(sf04,f04);
		ees.solve(sf05,f05);
		ees.solve(sf06,f06);
		ees.solve(sf07,f07);
		ees.solve(sf09,f09);
		ees.solve(sf10,f10);
		ees.solve(sf11,f11);
		ees.solve(sftf,fat);
		ees.solve(sfdf,fad);
		ees.solve(sffr,faf);

		/*** SOR with boundary ***/
		emdl.initialSFBoundary(sfbd,vm,wm);
		emdl.initialSFBoundary(sfsm,vm,wm);
		ees.solve(sfbd,null);
		ees.solve(sfsm,faf);
		
		
		/** calculate radial, vertical velocity components **/
		Variable[] vvfr=emdl.cVW(sffr);
		vvfr[0].setName("vsfr");	vvfr[0].setCommentAndUnit("radial wind forced by all internal forcings (m s^-1)");
		vvfr[1].setName("wsfr");	vvfr[1].setCommentAndUnit("omega forced by all internal forcings (Pa s^-1)");
		Variable[] vv01=emdl.cVW(sf01);
		vv01[0].setName("vs01");	vv01[0].setCommentAndUnit("radial wind forced by latent heating (m s^-1)");
		vv01[1].setName("ws01");	vv01[1].setCommentAndUnit("omega forced by latent heating (Pa s^-1)");
		Variable[] vv04=emdl.cVW(sf04);
		vv04[0].setName("vs04");	vv04[0].setCommentAndUnit("radial wind forced by eddy heat HFC (m s^-1)");
		vv04[1].setName("ws04");	vv04[1].setCommentAndUnit("omega forced by eddy heat HFC (Pa s^-1)");
		Variable[] vv05=emdl.cVW(sf05);
		vv05[0].setName("vs05");	vv05[0].setCommentAndUnit("radial wind forced by eddy heat VFC (m s^-1)");
		vv05[1].setName("ws05");	vv05[1].setCommentAndUnit("omega forced by eddy heat VFC (Pa s^-1)");
		Variable[] vv06=emdl.cVW(sf06);
		vv06[0].setName("vs06");	vv06[0].setCommentAndUnit("radial wind forced by eddy AAM HFC (m s^-1)");
		vv06[1].setName("ws06");	vv06[1].setCommentAndUnit("omega forced by eddy AAM HFC (Pa s^-1)");
		Variable[] vv07=emdl.cVW(sf07);
		vv07[0].setName("vs07");	vv07[0].setCommentAndUnit("radial wind forced by eddy AAM VFC (m s^-1)");
		vv07[1].setName("ws07");	vv07[1].setCommentAndUnit("omega forced by eddy AAM VFC (Pa s^-1)");
		Variable[] vv09=emdl.cVW(sf09);
		vv09[0].setName("vs09");	vv09[0].setCommentAndUnit("radial wind forced by horizontal viscous torque (m s^-1)");
		vv09[1].setName("ws09");	vv09[1].setCommentAndUnit("omega forced by horizontal viscous torque (Pa s^-1)");
		Variable[] vv10=emdl.cVW(sf10);
		vv10[0].setName("vs10");	vv10[0].setCommentAndUnit("radial wind forced by vertical viscous torque (m s^-1)");
		vv10[1].setName("ws10");	vv10[1].setCommentAndUnit("omega forced by vertical viscous torque (Pa s^-1)");
		Variable[] vv11=emdl.cVW(sf11);
		vv11[0].setName("vs11");	vv11[0].setCommentAndUnit("radial wind forced by tilting effect (m s^-1)");
		vv11[1].setName("ws11");	vv11[1].setCommentAndUnit("omega forced by tilting effect (Pa s^-1)");
		Variable[] vvtf=emdl.cVW(sftf);
		vvtf[0].setName("vstf");	vvtf[0].setCommentAndUnit("radial wind forced by all thermal forces (m s^-1)");
		vvtf[1].setName("wstf");	vvtf[1].setCommentAndUnit("omega forced by all thermal forces (Pa s^-1)");
		Variable[] vvdf=emdl.cVW(sfdf);
		vvdf[0].setName("vsdf");	vvdf[0].setCommentAndUnit("radial wind forced by all dynamic forces (m s^-1)");
		vvdf[1].setName("wsdf");	vvdf[1].setCommentAndUnit("omega forced by all dynamic forces (Pa s^-1)");
		Variable[] vvbd=emdl.cVW(sfbd);
		vvbd[0].setName("vsbd");	vvbd[0].setCommentAndUnit("radial wind forced by boundary effect (m s^-1)");
		vvbd[1].setName("wsbd");	vvbd[1].setCommentAndUnit("omega forced by boundary effect (Pa s^-1)");
		Variable[] vvsm=emdl.cVW(sfsm );
		vvsm[0].setName("vsm" );	vvsm[0].setCommentAndUnit("radial wind forced by all (m s^-1)");
		vvsm[1].setName("wsm" );	vvsm[1].setCommentAndUnit("omega forced by all (Pa s^-1)");
		
		
		/******************** data output *******************/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+ty.getID()+"_model.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),
			  sm, gm ,gdm  , qm  , lhm ,Qthm , um , vm , wm , Am , Bm , Cm ,BsinB,CsinB, Tm ,Tvm ,thm ,them, zm ,APrime,BPrime,CPrime,
			fvhm,fvvm,tavam,tawam,gavam,gawam,gavaR,gawaR,dhr ,htHFC,htVFC,aamHFC,aamVFC,aamHFCR,aamVFCR,fvhmtq,fvvmtq,tilt,APrimeC,BPrimeC,CPrimeC,
			f01f,f04f,f05f,f06f,f07f,f09f,f10f,f11f,sfe,sfe2,vedy[0],vedy[1],vedy2[0],vedy2[1],EPVector[0],EPVector[1],EPVector[2],EPVector[3],EPDiv,EPDivR,
			 faf, fat, fad, f01, f04, f05, f06, f07, f09, f10, f11,
			sffr,sftf,sfdf,sf01,sf04,sf05,sf06,sf07,sf09,sf10,sf11,sfbd,sfsm,
			vvfr[0],vv01[0],vv04[0],vv05[0],vv06[0],vv07[0],vv09[0],vv10[0],vv11[0],vvtf[0],vvdf[0],vvbd[0],vvsm[0],
			vvfr[1],vv01[1],vv04[1],vv05[1],vv06[1],vv07[1],vv09[1],vv10[1],vv11[1],vvtf[1],vvdf[1],vvbd[1],vvsm[1]
		);
		cdws.closeFile();
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	/**************************** data pre-process **********************************/
	static void verticalInterpolation(Typhoon ty,String path,boolean deleteOld){
		if(deleteOld)
		try{
			Files.deleteIfExists(Paths.get(path+ty.getID()+".dat"));
			Files.deleteIfExists(Paths.get(path+ty.getID()+".ctl"));
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		// interpolation
		DataInterpolation di=new DataInterpolation(DiagnosisFactory.getDataDescriptor(path+ty.getID()+"dl.ctl"));
		di.verticalInterp(path+ty.getID()+".dat",77,"q","t","u","v","w","z");
		//di.verticalInterp(path+ty.getID()+"_ln.dat",73,"q","t","u","v","w","z");
	}
	
	static void extractGribData(Typhoon ty,String path,int levs){
		MDate tstr=new MDate(ty.getTime(0));
		
		int tt=1;
		
		if(tstr.getHour()==0 ) tt=1;
		else if(tstr.getHour()==6 ) tt=2;
		else if(tstr.getHour()==12) tt=3;
		else if(tstr.getHour()==18) tt=4;
		else throw new IllegalArgumentException("invalid hour for "+tstr);
		
		StringBuilder sb=new StringBuilder();
		
		String id=ty.getID();
		String year="20"+id.substring(0,2);
		
		if(id.length()>4) year=id.substring(0,4);
		
		sb.append("'open "+path+id+"pl.ctl'\n");
		sb.append("'open "+path+id+"sfc.ctl'\n");
		sb.append("'sdfopen "+Paths.get(path).getParent().toString().replaceAll("\\\\","/")+"/blh."+year+".nc'\n\n");
		
		sb.append("'set gxout fwrite'\n");
		sb.append("'set fwrite "+path+id+"dl.dat'\n\n");
		
		sb.append("tt="+tt+"\n");
		sb.append("while(tt<="+(ty.getTCount()+tt-1)+")\n");
		sb.append("'set t 'tt\n");
		sb.append("'d no10Usfc.2'\n");
		sb.append("'d no10Vsfc.2'\n");
		sb.append("'d MSLsfc.2'\n");
		sb.append("'d SPsfc.2'\n");
		sb.append("'d SSTKsfc.2'\n");
		sb.append("'d Zsfc.2'\n");
		sb.append("'d blh.3(z=1)'\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Qprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Tprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Uprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Vprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Wprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("zz=1\n");
		sb.append("while(zz<="+levs+")\n");
		sb.append("'set z 'zz\n");
		sb.append("'d Zprs'\n");
		sb.append("zz=zz+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("tt=tt+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("'disable fwrite'\n");
		sb.append("'close 3'\n");
		sb.append("'close 2'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		sb.append("'quit'\n");
		
		try(FileWriter fw=new FileWriter(new File(path+"extract.gs"))){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		OpenGrADS.runGS(path+"extract.gs");
		
		sb=new StringBuilder();
		
		sb.append("dset ^"+id+"dl.dat\n");
		sb.append("undef -9.99E8\n");
		sb.append("title "+id+" TC\n");
		
		if(id.length()>4){
			sb.append("ydef  97 linear  -9.0 0.75\n");
			sb.append("xdef 169 linear 231.0 0.75\n");
		}else{
			sb.append("ydef 133 linear -15.0 0.75\n");
			sb.append("xdef 189 linear  81.0 0.75\n");
		}
		sb.append("tdef "+ty.getTCount()+" linear "+new MDate(ty.getTime(0)).toGradsDate()+" 6hr\n");
		if(levs==27)
		sb.append("zdef 27 levels 1000 975 950 925 900 875 850 825 800 775 750 700 650 600 550 500 450 400 350 300 250 225 200 175 150 125 100\n");
		else if(levs==29)
		sb.append("zdef 29 levels 1000 975 950 925 900 875 850 825 800 775 750 700 650 600 550 500 450 400 350 300 250 225 200 175 150 125 100 70 50\n");
		else throw new IllegalArgumentException("invalid levs [27 or 29]");
		sb.append("vars 13\n");
		sb.append("u10  0  99  ** 10 metre U wind component [m s**-1]\n");
		sb.append("v10  0  99  ** 10 metre V wind component [m s**-1]\n");
		sb.append("mslp 0  99  ** Mean sea level pressure [Pa]\n");
		sb.append("psfc 0  99  ** Surface pressure [Pa]\n");
		sb.append("sst  0  99  ** Sea surface temperature [K]\n");
		sb.append("zsfc 0  99  ** surface Geopotential [m**2 s**-2]\n");
		sb.append("blh  0  99  ** boundary layer height [m]\n");
		sb.append("q   "+levs+"  99  ** (profile) Specific humidity [kg kg**-1]\n");
		sb.append("T   "+levs+"  99  ** (profile) Temperature [K]\n");
		sb.append("u   "+levs+"  99  ** (profile) U velocity [m s**-1]\n");
		sb.append("v   "+levs+"  99  ** (profile) V velocity [m s**-1]\n");
		sb.append("w   "+levs+"  99  ** (profile) Vertical velocity [Pa s**-1]\n");
		sb.append("z   "+levs+"  99  ** (profile) Geopotential [m**2 s**-2]\n");
		sb.append("endvars\n");
		
		try(FileWriter fw=new FileWriter(new File(path+id+"dl.ctl"))){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void prepareGrib(){
		String path="/lustre/home/qianyk/Data/VortexStat/JMA/";
		
		int[] count=new int[]{23,26,23,21,29,23,23,24,22,22,14,21,20};
		
		StringBuffer sb=new StringBuffer();
		sb.append("#!/bin/bash\n");
		
		for(int yy=2000;yy<2013;yy++){
			sb.append("cd "+yy+"\n");
			
			for(int i=1;i<=count[yy-2000];i++){
				String id=String.format(String.valueOf(yy).substring(2)+"%02d",i);
				String pfile=id+"pl.grb";
				String sfile=id+"sfc.grb";
				String pctl =id+"pl.ctl";
				String sctl =id+"sfc.ctl";
				
				sb.append("mkdir "+id+"\n");
				sb.append("mv "+pfile+" "+id+"\n");
				sb.append("mv "+sfile+" "+id+"\n");
				sb.append("cd "+id+"\n");
				sb.append("grib2ctl.pl "+pfile+" > "+pctl+"\n");
				sb.append("grib2ctl.pl "+sfile+" > "+sctl+"\n");
				sb.append("gribmap -i "+pctl+"\n");
				sb.append("gribmap -i "+sctl+"\n");
				sb.append("cd ../\n\n");
			}
			
			sb.append("cd ../\n\n\n");
		}
		
		try(FileWriter fw=new FileWriter(path+"exGrib.bsh")){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/*** test **
	public static void main(String[] args){
		DataSets ds=DataSets.JMA;
		String yy="2004";
		String id="0420";
		String path="/lustre/home/qianyk/Data/";
		List<Typhoon> all=AccessBestTrack.getTyphoons(path+"/Typhoons/JMA/JMA.txt","",ds);
		
		List<Typhoon> ls=AccessBestTrack.getTyphoons(all,"id="+id);
		if(ls.size()!=1) throw new IllegalArgumentException("find more than one TC");
		Typhoon ty=ls.get(0);
		
		MainProcess.extractGribData(ty,path+"VortexStat/JMA/"+yy+"/"+id+"/",27);
		MainProcess.verticalInterpolation(ty,path+"VortexStat/JMA/"+yy+"/"+id+"/",true);
	}*/
}
