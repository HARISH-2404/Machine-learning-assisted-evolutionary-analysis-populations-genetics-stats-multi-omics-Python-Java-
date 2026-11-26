package org.ipgml.quantitativegenetics;
import java.util.*;
/** Independent ridge fits for each trait; not a full covariance-aware multivariate model. */
public final class MultiTraitRrblup {
 public static double[][] fit(double[][] X,double[][] Y,double lambda){
  int n=X.length;if(n==0||Y.length!=n)throw new IllegalArgumentException("Row mismatch");int t=Y[0].length;double[][] out=new double[n][t];
  for(int j=0;j<t;j++){double[] y=new double[n];for(int i=0;i<n;i++)y[i]=Y[i][j];double[] pred=(double[])Rrblup.fit(X,y,lambda).get("predictions");for(int i=0;i<n;i++)out[i][j]=pred[i];}return out;
 }
}
