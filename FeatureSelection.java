package org.ipgml.machinelearning;
import java.util.*;
public final class FeatureSelection {
 public static double[] absoluteCorrelations(double[][] X,double[] y){int n=X.length,p=X[0].length;if(y.length!=n)throw new IllegalArgumentException("Length mismatch");double my=Arrays.stream(y).average().orElseThrow(),vy=0;for(double v:y)vy+=(v-my)*(v-my);double[]s=new double[p];for(int j=0;j<p;j++){double mx=0;for(double[] row:X)mx+=row[j]/n;double vx=0,c=0;for(int i=0;i<n;i++){vx+=(X[i][j]-mx)*(X[i][j]-mx);c+=(X[i][j]-mx)*(y[i]-my);}s[j]=vx==0||vy==0?0:Math.abs(c/Math.sqrt(vx*vy));}return s;}
}
