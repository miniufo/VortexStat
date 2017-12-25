package pack;
//
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.stream.Stream;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Typhoon;
import miniufo.util.GridDataFetcher;
import miniufo.util.TCTracker;
import miniufo.util.TicToc;


//
public class TrackAndIntensity{
	//
	static final DataSets ds=DataSets.JMA;
	static final String path="/lustre/home/qianyk/Data/";
	
	static final List<Typhoon> all=AccessBestTrack.getTyphoons(path+"/Typhoons/JMA/JMA.txt","",ds);
	
	//
	public static void main(String[] args){
		//generateBSTFile();
		combineBSTFile();
	}
	
	public static void combineBSTFile(){
		int[] count=new int[]{23,26,23,21,29,23,23,24,22,22,14,21,20};
		
		FileWriter fw=null;
		
		try{ fw=new FileWriter(path+"VortexStat/JMA/All.bst");}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		for(int yy=2000;yy<2013;yy++){
			TicToc.tic("computing tracks for "+yy);
			
			for(int i=1;i<=count[yy-2000];i++){
				String id=String.format(String.valueOf(yy).substring(2)+"%02d",i);
				
				try(BufferedReader br=new BufferedReader(new FileReader(path+"VortexStat/JMA/"+yy+"/"+id+"/"+id+".bst"))){
					Stream<String> lines=br.lines();
					
					final FileWriter tmp=fw;
					
					Consumer<String> write=s->{
						try{ tmp.write(s+"\n");}
						catch(IOException e){ e.printStackTrace(); System.exit(0);}
					};
					
					lines.filter(s->s.indexOf("*")==-1).forEach(write);
					
				}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			}
			
			TicToc.toc(TimeUnit.SECONDS);
		}
		
		try{fw.close();}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void generateBSTFile(){
		int[] count=new int[]{23,26,23,21,29,23,23,24,22,22,14,21,20};
		
		for(int yy=2000;yy<2013;yy++){
			TicToc.tic("computing tracks for "+yy);
			
			for(int i=1;i<=count[yy-2000];i++){
				String id=String.format(String.valueOf(yy).substring(2)+"%02d",i);
				
				List<Typhoon> ls=AccessBestTrack.getTyphoons(all,"id="+id);
				
				if(ls.size()!=1) throw new IllegalArgumentException("find more than one TC");
				
				Typhoon ty=ls.get(0);
				
				DiagnosisFactory df=DiagnosisFactory.parseFile(path+"VortexStat/JMA/"+yy+"/"+id+"/"+id+"_ln.ctl");
				DataDescriptor dd=df.getDataDescriptor();
				
				TCTracker tt=new TCTracker(dd);
				tt.findTrackBySLP("mslp",ty,2);
				
				float[] mslp=getAlongTrackSLP(dd,ty);
				
				toTrackFile(path+"VortexStat/JMA/"+yy+"/"+id+"/"+id+".bst",ty,tt,mslp);
			}
			
			TicToc.toc(TimeUnit.SECONDS);
		}
	}
	
	public static float[] getAlongTrackSLP(DataDescriptor dd,Typhoon ty){
		float[] mslp=new float[ty.getTCount()];
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		for(int l=0,L=ty.getTCount();l<L;l++){
			Variable buf=gdf.prepareXYBuffer("mslp",l+1,1);
			mslp[l]=gdf.fetchXYBuffer(ty.getXPosition(l),ty.getYPosition(l),buf);
		}
		
		gdf.closeFile();
		
		return mslp;
	}
	
	public static void toTrackFile(String fname,Typhoon ty,TCTracker tt,float[] mslp){
		int tcount=ty.getTCount();
		
		if(tcount!=tt.getTCount()) throw new IllegalArgumentException("T-lengths not equals");
		
		StringBuilder sb=new StringBuilder();
		
		float[] lons=tt.getLongitudes();
		float[] lats=tt.getLatitudes();
		
		float[] press=ty.getPressures();
		float[] minP=tt.getMinimumSLPs();
		
		sb.append("*   "+tcount+" "+ty.getName()+" "+ty.getID()+"\n");
		for(int l=0;l<tcount;l++){
			
			String s=String.format("%7.2f  %7.2f  %8.2f  %8.2f  %6.2f  %d    %7.2f  %7.2f  %8.2f    %6.2f  %6.2f  %6.2f\n",
					ty.getXPosition(l),ty.getYPosition(l),press[l],mslp[l]/100f,(mslp[l]/100f-press[l]),ty.getTime(l),
					lons[l],lats[l],minP[l]/100f,
					(lons[l]-ty.getXPosition(l)),(lats[l]-ty.getYPosition(l)),(minP[l]/100f-press[l])
				);
		sb.append(s);
		}
		
		try(FileWriter fw=new FileWriter(fname)){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
}
