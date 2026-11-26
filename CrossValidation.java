package org.ipgml.machinelearning;
import java.util.*;
public final class CrossValidation {
 public static List<int[][]> kfold(int n,int k,long seed){if(k<2||k>n)throw new IllegalArgumentException("Invalid folds");List<Integer> ids=new ArrayList<>();for(int i=0;i<n;i++)ids.add(i);Collections.shuffle(ids,new Random(seed));List<int[][]> out=new ArrayList<>();for(int f=0;f<k;f++){List<Integer>tr=new ArrayList<>(),te=new ArrayList<>();for(int i=0;i<n;i++)(i%k==f?te:tr).add(ids.get(i));out.add(new int[][]{tr.stream().mapToInt(Integer::intValue).toArray(),te.stream().mapToInt(Integer::intValue).toArray()});}return out;}
}
