//
package pack;

import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CsmDataWriteStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;


//
public class StormRelative{
	//
	private static final String name="Haima";
	private static final String path="D:/Data/DiagnosisVortex/"+name;
	
	private static final Typhoon ty=AccessBestTrack.getTyphoons(
		//"d:/Data/Typhoons/NHC/NHC.txt","name="+name+";time=1Aug1985-30Sep1985;",DataSets.NHC
		"d:/Data/Typhoons/CMA/CMA.txt","name="+name+";time=1sep2004-30Sep2004;",DataSets.CMA
	).get(0);
	
	
	// modelling
	static void modeling(Typhoon ty){
		CsmDescriptor csd=(CsmDescriptor)DiagnosisFactory.parseContent(
			ty.toCSMString(path+"/Interim/Prs/"+name+".ctl",144,91,37,0.15f,-25f,1000)
		).getDataDescriptor();
		
		SphericalSpatialModel   ssm=new SphericalSpatialModel(csd.getCtlDescriptor());
		CylindricalSpatialModel csm=new CylindricalSpatialModel(csd);
		DynamicMethodsInCC dmdl=new DynamicMethodsInCC(csm);
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		
		for(int l=1;l<ty.getTCount()-1;l++) System.out.println((l+1)+"\t"+csm.getUWhole()[l]+"\t"+csm.getVWhole()[l]+"\t"+
			Math.toRadians(csd.getOLon()[l+1]-csd.getOLon()[l-1])/(2.0*3600*6)*1500000.0*Math.sin(Math.toRadians(csd.getOLat()[l]))
		);
		
		Variable[] Wanaly=dmdl.cInhomoTranslatingVelocityByAnaly();
		Variable[] WDiff =dmdl.cInhomoTranslatingVelocityByDiff();
		
		CsmDataWriteStream cdw=new CsmDataWriteStream(path+"/StormRelative/cylind.dat");
		cdw.writeData(csd,Wanaly[0],Wanaly[1],WDiff[0],WDiff[1]); cdw.closeFile();
		
		Wanaly=ct.reprojectToCylindrical(Wanaly[0],Wanaly[1]);
		WDiff =ct.reprojectToCylindrical( WDiff[0], WDiff[1]);
		
		Variable vtmA=Wanaly[0].anomalizeX(); vtmA.setName("vtmA");
		Variable vrmA=Wanaly[1].anomalizeX(); vrmA.setName("vrmA");
		Variable vtmD= WDiff[0].anomalizeX(); vtmD.setName("vtmD");
		Variable vrmD= WDiff[1].anomalizeX(); vrmD.setName("vrmD");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"/StormRelative/mean.dat"); cdws.setPrinting(false);
		cdws.writeData(csd.getCtlDescriptor(),vtmA,vrmA,vtmD,vrmD); cdws.closeFile();
	}
	
	
	/*** test ***/
	public static void main(String[] args){
		System.out.println(ty);
		
		modeling(ty);
	}
}
