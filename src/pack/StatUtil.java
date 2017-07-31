package pack;
//
import java.util.function.Predicate;


public final class StatUtil{
	//
	public static float[][] mapToBWLabelData(float[][] data,float undef,Predicate<Float> condition){
		int z=data.length,y=data[0].length;
		
		float[][] BWlabel=new float[z][y];
		
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++)
		if(condition.test(data[k][j])&&data[k][j]!=undef) BWlabel[k][j]=1;
		
		return BWlabel;
	}
	
	public static String getVPageString(int row,int col,int idx){
		int total=row*col;
		
		if(idx< 0    ) throw new IllegalArgumentException("idx start from 0");
		if(idx>=total) throw new IllegalArgumentException("idx ("+idx+") should be smaller than row*col ("+total+")");
		
		int rowC=idx/col;
		int colC=idx%col+1;
		
		return row+" "+col+" "+(row-rowC)+" "+colC;
	}
}
