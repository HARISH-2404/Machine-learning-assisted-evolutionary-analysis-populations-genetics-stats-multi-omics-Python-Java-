package org.ipgml.transcriptomics;
import java.util.*;

/** Relative isoform usage from transcript abundance matrix. */
public class IsoformUsage {
  public static double[][] proportions(double[][] abundance){
    int n=abundance.length;if(n==0)return new double[0][0];int samples=abundance[0].length;
    double[][] out=new double[n][samples];
    for(int j=0;j<samples;j++){double total=0;for(int i=0;i<n;i++){if(abundance[i].length!=samples||abundance[i][j]<0||!Double.isFinite(abundance[i][j]))throw new IllegalArgumentException("Invalid abundance matrix");total+=abundance[i][j];}
      if(total>0)for(int i=0;i<n;i++)out[i][j]=abundance[i][j]/total;
    }return out;
  }
}
