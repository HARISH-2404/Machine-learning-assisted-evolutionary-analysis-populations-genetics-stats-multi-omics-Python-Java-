package org.ipgml.machinelearning;
import java.util.*;
public final class ModelMetrics {
 public static Map<String,Double> regression(double[] y,double[] p){if(y.length!=p.length||y.length==0)throw new IllegalArgumentException("Length mismatch");double mse=0,mae=0,mean=Arrays.stream(y).average().orElseThrow(),sst=0;for(int i=0;i<y.length;i++){double e=y[i]-p[i];mse+=e*e;mae+=Math.abs(e);sst+=(y[i]-mean)*(y[i]-mean);}mse/=y.length;mae/=y.length;Map<String,Double>m=new LinkedHashMap<>();m.put("mse",mse);m.put("rmse",Math.sqrt(mse));m.put("mae",mae);m.put("r_squared",sst==0?Double.NaN:1-mse*y.length/sst);return m;}
 public static Map<String,Double> classification(int[] y,int[] p){if(y.length!=p.length||y.length==0)throw new IllegalArgumentException("Length mismatch");int tp=0,tn=0,fp=0,fn=0;for(int i=0;i<y.length;i++){if(y[i]==1&&p[i]==1)tp++;else if(y[i]==0&&p[i]==0)tn++;else if(y[i]==0)fp++;else fn++;}double prec=tp+fp==0?0:tp/(double)(tp+fp),rec=tp+fn==0?0:tp/(double)(tp+fn);Map<String,Double>m=new LinkedHashMap<>();m.put("accuracy",(tp+tn)/(double)y.length);m.put("precision",prec);m.put("recall",rec);m.put("f1",prec+rec==0?0:2*prec*rec/(prec+rec));return m;}
}
