/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import java.util.HashSet;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to compute various forms of location/distance based features
 */
public class LocationSimilarity {
    HashMap<String,Pair<ArrayList<Double>,ArrayList<Double>>> locationN=new HashMap<>();
    HashMap<String,Integer> COGIndex=new HashMap<>();
    HashMap<Integer,String> IndexCOG=new HashMap<>();
    
    void initializeSimilarity(COGGOMap cgmap){
        Iterator<String> iterator=cgmap.CogGOmap.keySet().iterator();
        int ind=0;
        
        while(iterator.hasNext()){
            String cog=iterator.next();
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
            ArrayList<Double> tmp1=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
            Pair<ArrayList<Double>,ArrayList<Double>> pr=new Pair(tmp,tmp1);
            locationN.put(cog, pr);
            COGIndex.put(cog, ind);
            IndexCOG.put(ind, cog);
            ind++;
        }
    }
    
    void computeSimilarityB(COG cog, COGGOMap cgmap, CircularDistance distance){  
        
        HashMap<Integer,Pair<HashSet<Integer>,ArrayList<Double>>> tmpDist=new HashMap<>();
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;    
       
        if(!tmpDist.containsKey(COGIndex.get(ck.getValue2()))){
            HashSet<Integer> tmp=new HashSet<>();
            ArrayList<Double> tA=new ArrayList<>(Collections.nCopies(COGIndex.keySet().size(), Double.POSITIVE_INFINITY));
            tmpDist.put(COGIndex.get(ck.getValue2()),new Pair<HashSet<Integer>,ArrayList<Double>>(tmp,tA));
        }
       
       for(int j=i+1;j<cog.cogs.size()-1;j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
            double dist=tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().get(COGIndex.get(ckN.getValue2()));
            double tmpDist1=distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate);
            if(tmpDist1<dist){
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().set(COGIndex.get(ckN.getValue2()),tmpDist1);
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue0().add(COGIndex.get(ckN.getValue2()));
                
                if(!tmpDist.containsKey(COGIndex.get(ckN.getValue2()))){
                    HashSet<Integer> tmp=new HashSet<>();
                     ArrayList<Double> tA=new ArrayList<>(Collections.nCopies(COGIndex.keySet().size(), Double.POSITIVE_INFINITY));
                     tmpDist.put(COGIndex.get(ckN.getValue2()),new Pair<HashSet<Integer>,ArrayList<Double>>(tmp,tA));
                  }
                
                tmpDist.get(COGIndex.get(ckN.getValue2())).getValue1().set(COGIndex.get(ck.getValue2()),tmpDist1);
                tmpDist.get(COGIndex.get(ckN.getValue2())).getValue0().add(COGIndex.get(ck.getValue2()));
                
            }
       }
      }
      
      for(int k:tmpDist.keySet()){
          HashSet<Integer> tmp=tmpDist.get(k).getValue0();
          for(int i:tmp){
              double dist=tmpDist.get(k).getValue1().get(i);
              double dist1=locationN.get(IndexCOG.get(k)).getValue0().get(i);
              double oldNum=locationN.get(IndexCOG.get(k)).getValue1().get(i);
              oldNum++;
              dist1+=dist;
              Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(IndexCOG.get(k));
       element.getValue1().set(i, oldNum);//.setAt1(number);
       element.getValue0().set(i, dist1);
       locationN.put(IndexCOG.get(k), element);
          }
      }
      
      System.out.println("Number of COGs in LS: "+locationN.keySet().size());
      
    }
    
    
     void computeLogSimilarityB(COG cog, COGGOMap cgmap, CircularDistance distance){  
        HashMap<Integer,Pair<HashSet<Integer>,ArrayList<Double>>> tmpDist=new HashMap<>();
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;
       
        if(!tmpDist.containsKey(COGIndex.get(ck.getValue2()))){
            HashSet<Integer> tmp=new HashSet<>();
            ArrayList<Double> tA=new ArrayList<>(Collections.nCopies(COGIndex.keySet().size(), Double.POSITIVE_INFINITY));
            tmpDist.put(COGIndex.get(ck.getValue2()),new Pair<HashSet<Integer>,ArrayList<Double>>(tmp,tA));
        }
       
       for(int j=0;j<cog.cogs.size();j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
          continue;
     
        
         if(tmpDist.containsKey(COGIndex.get(ck.getValue2()))){
            HashSet<Integer> tS=tmpDist.get(COGIndex.get(ck.getValue2())).getValue0();
            if(tS.contains(COGIndex.get(ckN.getValue2()))){
            double dist=tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().get(COGIndex.get(ckN.getValue2()));
            double tmpDist1=Math.log(distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate));
            if(tmpDist1<dist)
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().set(COGIndex.get(ckN.getValue2()),tmpDist1);
            }
            else{
             tS.add(COGIndex.get(ckN.getValue2()));
             double tmpDist1=Math.log(distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate));
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().set(COGIndex.get(ckN.getValue2()),tmpDist1);
            }   
        }
       }
      }
      
      
      for(int k:tmpDist.keySet()){
          HashSet<Integer> tmp=tmpDist.get(k).getValue0();
          for(int i:tmp){
              double dist=tmpDist.get(k).getValue1().get(i);
              double dist1=locationN.get(IndexCOG.get(k)).getValue0().get(i);
              double oldNum=locationN.get(IndexCOG.get(k)).getValue1().get(i);
              oldNum++;
              dist1+=dist;
              Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(IndexCOG.get(k));
       element.getValue1().set(i, oldNum);//.setAt1(number);
       element.getValue0().set(i, dist1);
       locationN.put(IndexCOG.get(k), element);
          }
      }  
      
    }
    
    void computeLogSimilarityLocB(COG cog, COGGOMap cgmap, CircularDistance distance){  
        
        HashMap<Integer,Pair<HashSet<Integer>,ArrayList<Double>>> tmpDist=new HashMap<>();
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;
       
        if(!tmpDist.containsKey(COGIndex.get(ck.getValue2()))){
            HashSet<Integer> tmp=new HashSet<>();
            ArrayList<Double> tA=new ArrayList<>(Collections.nCopies(COGIndex.keySet().size(), Double.POSITIVE_INFINITY));
            tmpDist.put(COGIndex.get(ck.getValue2()),new Pair<HashSet<Integer>,ArrayList<Double>>(tmp,tA));
        }
       
       for(int j=0;j<cog.cogs.size();j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
        
         if(tmpDist.containsKey(COGIndex.get(ck.getValue2()))){
            HashSet<Integer> tS=tmpDist.get(COGIndex.get(ck.getValue2())).getValue0();
            if(tS.contains(COGIndex.get(ckN.getValue2()))){
            double dist=tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().get(COGIndex.get(ckN.getValue2()));
            double tmpDist1=distance.computeDistance(Math.log(ck.getValue0()), Math.log(ck.getValue1()), Math.log(ckN.getValue0()), Math.log(ckN.getValue1()),Math.log(cog.maxCoordinate));
            if(tmpDist1<dist)
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().set(COGIndex.get(ckN.getValue2()),tmpDist1);
            }
            else{
             tS.add(COGIndex.get(ckN.getValue2()));
             double tmpDist1=distance.computeDistance(Math.log(ck.getValue0()), Math.log(ck.getValue1()), Math.log(ckN.getValue0()), Math.log(ckN.getValue1()),Math.log(cog.maxCoordinate));
                tmpDist.get(COGIndex.get(ck.getValue2())).getValue1().set(COGIndex.get(ckN.getValue2()),tmpDist1);
            }   
        }
       }
      }
      
      
      for(int k:tmpDist.keySet()){
          HashSet<Integer> tmp=tmpDist.get(k).getValue0();
          for(int i:tmp){
              double dist=tmpDist.get(k).getValue1().get(i);
              double dist1=locationN.get(IndexCOG.get(k)).getValue0().get(i);
              double oldNum=locationN.get(IndexCOG.get(k)).getValue1().get(i);
              oldNum++;
              dist1+=dist;
              Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(IndexCOG.get(k));
       element.getValue1().set(i, oldNum);//.setAt1(number);
       element.getValue0().set(i, dist1);
       locationN.put(IndexCOG.get(k), element);
          }
      }     
    }
    
    void computeSimilarity(COG cog, COGGOMap cgmap,GeneOGMapping geneOGMap ,CircularDistance distance){  
        
        HashSet<String> contCogs=new HashSet<>();
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.geneOGsMap.containsKey(ck.getValue2())))
           continue;
       
       HashSet<String> ogs=geneOGMap.geneOGsMap.get(ck.getValue2());
       
       for(int j=0;j<cog.cogs.size();j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.geneOGsMap.containsKey(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
        HashSet<String> ogs1=geneOGMap.geneOGsMap.get(ckN.getValue2());
     
      for(String ogS:ogs){ 
          if(!COGIndex.containsKey(ogS))
              continue;
          for(String ogS1:ogs1){
             if(!COGIndex.containsKey(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       dist+=distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate);
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
         }
        }
       }
      }
    }
    
    void computeLogSimilarity(COG cog, COGGOMap cgmap, GeneOGMapping geneOGMap, CircularDistance distance){  
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.geneOGsMap.containsKey(ck.getValue2())))
           continue;
       
       HashSet<String> ogs=geneOGMap.geneOGsMap.get(ck.getValue2());
       
       for(int j=0;j<cog.cogs.size();j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.geneOGsMap.containsKey(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
        HashSet<String> ogs1=geneOGMap.geneOGsMap.get(ckN.getValue2());

         for(String ogS:ogs){ 
          if(!COGIndex.containsKey(ogS))
              continue;
          for(String ogS1:ogs1){
             if(!COGIndex.containsKey(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       dist+=Math.log10(distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate))/Math.log10(2);//ln vs log2
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
         }
       }
       }
      }
    }
    
    void computeLogSimilarityLoc(COG cog, COGGOMap cgmap, CircularDistance distance){  
        
      for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;
       for(int j=0;j<cog.cogs.size();j++){
          Triplet<Integer,Integer,String> ckN=cog.cogs.get(j);
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
       double dist=locationN.get(ck.getValue2()).getValue0().get(COGIndex.get(ckN.getValue2()));
       double number=locationN.get(ck.getValue2()).getValue1().get(COGIndex.get(ckN.getValue2()));
       number=number+1.0;
       dist+=distance.computeDistance(Math.log(ck.getValue0()), Math.log(ck.getValue1()), Math.log(ckN.getValue0()), Math.log(ckN.getValue1()),Math.log(cog.maxCoordinate));
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ck.getValue2());
       element.getValue1().set(COGIndex.get(ckN.getValue2()), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ckN.getValue2()), dist);
       locationN.put(ck.getValue2(), element);
       }
      }
    }
    
    void normalize(){
        Iterator<String> cogIt=locationN.keySet().iterator();
        
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(element.getValue1().get(i)!=0){
                 dist/=element.getValue1().get(i);
                 element.getValue0().set(i, dist);
                 }
                 else if(i!=COGIndex.get(cog))
                   element.getValue0().set(i,Double.POSITIVE_INFINITY);  
                 else if(i==COGIndex.get(cog))
                     element.getValue0().set(i, 0.0);
            }
             locationN.put(cog, element);
        }
    }
    
    void normalize1(){
        Iterator<String> cogIt=locationN.keySet().iterator();
        
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(element.getValue1().get(i)!=0){
                 dist/=element.getValue1().get(i);
                 element.getValue0().set(i, dist);
                 } 
                 else 
                     element.getValue0().set(i, 0.0);
            }
             locationN.put(cog, element);
        }
    }
    
    void transformToWeights(){
        Iterator<String> cogIt=locationN.keySet().iterator();
        double maxDist=Double.NEGATIVE_INFINITY;
        
         while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
               if(dist>maxDist)
                   maxDist=dist;
            }
        }
        
         cogIt=locationN.keySet().iterator();
         
         while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(dist!=0.0){
                     dist/=maxDist;
                     dist=1.0-dist+Math.pow(10, -10);
                     element.getValue0().set(i, dist);
                 }
                if(i==COGIndex.get(cog))
                     element.getValue0().set(i, 1.0);
            }
        }   
    }
    
     void normalize(int type){
        Iterator<String> cogIt=locationN.keySet().iterator();
        
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(element.getValue1().get(i)!=0){
                 dist/=element.getValue1().get(i);
                 element.getValue0().set(i, dist);
                 }
                 else if(i!=COGIndex.get(cog))
                   element.getValue0().set(i,Double.POSITIVE_INFINITY);  
                 else if(i==COGIndex.get(cog) && (type==0 || type>1))
                     element.getValue0().set(i, 0.0);
                 else if(i==COGIndex.get(cog) && type==1)
                     element.getValue0().set(i, 0.0/*Math.log(Math.pow(10.0, -10.0))*/);                     
            }
             locationN.put(cog, element);
        }
    }
    
}
