package org.ipgml.genomics;
/** VanRaden-style additive genomic relationship matrix for diploid 0/1/2 dosages. */
public class GenomicRelationship {
  public static double[][] vanRaden(double[][] g){
    int n=g.length,m=g[0].length;double[] p=new double[m];boolean[] keep=new boolean[m];double denom=0;
    for(int j=0;j<m;j++){double sum=0;int c=0;for(int i=0;i<n;i++)if(Double.isFinite(g[i][j])){sum+=g[i][j];c++;}
      p[j]=c==0?Double.NaN:sum/(2*c);keep[j]=Double.isFinite(p[j])&&p[j]>0&&p[j]<1;
      if(keep[j])denom+=2*p[j]*(1-p[j]);
    }
    if(denom<=0)throw new IllegalArgumentException("No polymorphic markers");
    double[][] k=new double[n][n];
    for(int a=0;a<n;a++)for(int b=a;b<n;b++){double s=0;
      for(int j=0;j<m;j++)if(keep[j]){
        double ga=Double.isFinite(g[a][j])?g[a][j]:2*p[j];
        double gb=Double.isFinite(g[b][j])?g[b][j]:2*p[j];
        s+=(ga-2*p[j])*(gb-2*p[j]);
      }k[a][b]=k[b][a]=s/denom;
    }return k;
  }
}
