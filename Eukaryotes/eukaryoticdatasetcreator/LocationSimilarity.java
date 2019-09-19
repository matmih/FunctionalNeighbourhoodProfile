/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import java.util.HashSet;


/**
 *
 * @author matej
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to containing information and functions to compute distance - based OG similarity
 */
public class LocationSimilarity {
    HashMap<String,Pair<ArrayList<Double>,ArrayList<Double>>> locationN=new HashMap<>();
    HashMap<String,Integer> COGIndex=new HashMap<>();
    HashMap<Integer,String> IndexCOG=new HashMap<>();
    int fullSize=-1;
    
    void initializeSimilarity(OGGOMapping cgmap){
        Iterator<String> iterator=cgmap.CogGOmap.keySet().iterator();
        int ind=0;
        
        while(iterator.hasNext()){
            String cog=iterator.next();
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
            ArrayList<Double> tmp1=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
            Pair<ArrayList<Double>,ArrayList<Double>> pr=new Pair(tmp,tmp1);
            locationN.put(new String(cog), pr);
            COGIndex.put(new String(cog), ind);
            IndexCOG.put(ind, new String(cog));
            ind++;
        }
        fullSize=ind;
    }
    
    void computeSimilarity(Gene gene, OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation, int taxID){  
        
        HashSet<String> contCogs=new HashSet<>();
        
      for(int i=0;i<gene.genes.size();i++){
       Triplet<Integer,Integer,String> ck=gene.genes.get(i);
       //System.out.println("Ime COG-a: "+ck.getValue2());
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ck.getValue2())))
           continue;
       
       HashSet<Pair<String,String>> OGsT=geneOGMap.get(ck.getValue2());
       
        HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
       
     //  contCogs.add(ck.getValue2());
       for(int j=i+1;j<gene.genes.size();j++){
          Triplet<Integer,Integer,String> ckN=gene.genes.get(j);
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
        HashSet<Pair<String,String>> OGsT1=geneOGMap.get(ckN.getValue2());
     
        
         HashSet<String> OGs1 = new HashSet<>();
            
            for(Pair<String,String> p:OGsT1){
                if(OGs.contains(p.getValue0()))
                    continue;
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs1.add(p.getValue0());
            }
       
        
      for(String ogS:OGs){ 
          if(!COGIndex.containsKey(ogS))
              continue;
          for(String ogS1:OGs1){
             if(!COGIndex.containsKey(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       double t = Math.abs(ckN.getValue1()-ck.getValue1());
       if(t!=Double.POSITIVE_INFINITY && t!=Double.NEGATIVE_INFINITY)
            dist+=Math.abs(ckN.getValue1()-ck.getValue1());
       //dist+=distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate);
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
        dist=locationN.get(ogS1).getValue0().get(COGIndex.get(ogS));
       number=locationN.get(ogS1).getValue1().get(COGIndex.get(ogS));
       number=number+1.0;
       t=Math.abs(ckN.getValue1()-ck.getValue1());;
       if(t!=Double.POSITIVE_INFINITY && t!=Double.NEGATIVE_INFINITY)
             dist+=Math.abs(ckN.getValue1()-ck.getValue1());
       element=locationN.get(ogS1);
       element.getValue1().set(COGIndex.get(ogS), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS), dist);
       locationN.put(ogS1, element);
         }
        }
       }
      }
    }
    
    
     
    void computeSimilaritynonCH(Gene gene, OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation, int taxID){  
        
        HashSet<String> contCogs=new HashSet<>();
      
     for(int g=0;g<gene.genesContig.size();g++){
         gene.genes=gene.genesContig.get(g);
      for(int i=0;i<gene.genes.size();i++){
       Triplet<Integer,Integer,String> ck=gene.genes.get(i);
       //System.out.println("Ime COG-a: "+ck.getValue2());
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ck.getValue2())))
           continue;
       
       HashSet<Pair<String,String>> OGsT=geneOGMap.get(ck.getValue2());
       
        HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
       
     //  contCogs.add(ck.getValue2());
       for(int j=i+1;j<gene.genes.size();j++){
          Triplet<Integer,Integer,String> ckN=gene.genes.get(j);
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
           continue;
     
        HashSet<Pair<String,String>> OGsT1=geneOGMap.get(ckN.getValue2());
     
        
         HashSet<String> OGs1 = new HashSet<>();
            
            for(Pair<String,String> p:OGsT1){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs1.add(p.getValue0());
            }
       
        
      for(String ogS:OGs){ 
          if(!COGIndex.containsKey(ogS))
              continue;
          for(String ogS1:OGs1){
             if(!COGIndex.containsKey(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       double t=Math.abs(ckN.getValue1()-ck.getValue1());
       if(t!=Double.POSITIVE_INFINITY && t!=Double.NEGATIVE_INFINITY)
       dist+=Math.abs(ckN.getValue1()-ck.getValue1());
       //dist+=distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate);
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
       
         dist=locationN.get(ogS1).getValue0().get(COGIndex.get(ogS));
       number=locationN.get(ogS1).getValue1().get(COGIndex.get(ogS));
       number=number+1.0;
        dist+=Math.abs(ckN.getValue1()-ck.getValue1());
       element=locationN.get(ogS1);
       element.getValue1().set(COGIndex.get(ogS), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS), dist);
       locationN.put(ogS1, element);
       
         }
        }
       }
      }
     }
    }
    
    void computeLogSimilarity(Gene gene, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation, int taxID){  
        
      for(int i=0;i<gene.genes.size();i++){
       Triplet<Integer,Integer,String> ck=gene.genes.get(i);
       //System.out.println("Ime COG-a: "+ck.getValue2());
       
      /* if(ck.getValue2().equals("MGG_03650")){
           System.out.println("gene found: "+ck.getValue2()+" taxId: "+taxID);
       }*/
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ck.getValue2())))
           continue;
       
      // HashSet<Pair<String,String>> ogs=geneOGMap.get(ck.getValue2());
       
       HashSet<Pair<String,String>> OGsT=geneOGMap.get(ck.getValue2());
       
        HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
            
            
           /*  if(ck.getValue2().equals("MGG_03650")){
           System.out.println("OGs Size:"+OGs.size()+" taxId: "+taxID);
       }*/
                
       /*if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;*/
       for(int j=i+1;j<gene.genes.size();j++){
          Triplet<Integer,Integer,String> ckN=gene.genes.get(j);
          
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ckN.getValue2())))
           continue;
       /* if(ck.getValue2().equals(ckN.getValue2()))
           continue;*/
     
       /* if(ck.getValue2().equals("MGG_03650")){
           System.out.println("N gene found: "+ckN.getValue2()+" taxId: "+taxID);
           
       }*/
        
       
       
        //HashSet<Pair<String,String>> ogs1=geneOGMap.get(ckN.getValue2());
          
         HashSet<Pair<String,String>> OGsT1=geneOGMap.get(ckN.getValue2());
     
        
         HashSet<String> OGs1 = new HashSet<>();
            
            for(Pair<String,String> p:OGsT1){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs1.add(p.getValue0());
            }
            
             /* if(ck.getValue2().equals("MGG_03650")){
           System.out.println("OGsN Size:"+OGs1.size()+" taxId: "+taxID);
       }*/
            
            /* if(ck.getValue2().equals("FOXG_05859") || (ck.getValue2().equals("FVEG_03736T0")) || (ck.getValue2().equals("FVEG_03736")) || (ck.getValue2().equals("FOXG_05859P0")) || (ck.getValue2().equals("NechaG86793")) || (ck.getValue2().equals("NechaG82508")) || (ck.getValue2().equals("NechaP86793")) || (ck.getValue2().equals("NechaP82508"))){
                 System.out.println(ck.getValue2()+" "+OGs.size()+" "+OGs1.size());
            }*/
          
          /*if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
          continue;*/
             
            /* if(OGs.contains("ENOG410PSIY")){
                 System.out.println("taxId: "+taxID);
                 System.out.println("gene: "+ck.getValue2());
                 System.out.println("neighbour gene: "+ckN.getValue2());
                 System.out.println("NOgs: ");
                 for(String og:OGs1)
                     System.out.print(og+" ");
                 System.out.println();
             }*/
             
             /*if((OGs.size()!=0 && OGs1.size()==0) || OGs.size()==0 && OGs1.size()!=0){
                 System.out.println("Empty nn: ");
                 for(String ss:OGs)
                     System.out.print(ss+" ");
                 System.out.println();
                 for(String ss:OGs1)
                     System.out.print(ss+" ");
                 System.out.println();
             }*/
     
         for(String ogS:OGs){ 
          //if(!COGIndex.containsKey(ogS))
             if(!cgmap.CogGOmap.keySet().contains(ogS))
              continue;
          for(String ogS1:OGs1){
           //  if(!COGIndex.containsKey(ogS1))
              if(!cgmap.CogGOmap.keySet().contains(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       if(dist==Double.POSITIVE_INFINITY || dist == Double.NEGATIVE_INFINITY){
           System.out.println("NFdist: "+dist);
           System.out.println("OG: "+ogS);
           System.out.println("OG1: "+ogS1);
       }
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       if(Math.abs(ckN.getValue1()-ck.getValue1())>=2)
        dist+=Math.log10(Math.abs(ckN.getValue1()-ck.getValue1()))/Math.log10(2);
       else dist+=2;
       //dist+=Math.log10(distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate))/Math.log10(2);//ln vs log2
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
       
        dist=locationN.get(ogS1).getValue0().get(COGIndex.get(ogS));
       number=locationN.get(ogS1).getValue1().get(COGIndex.get(ogS));
       number=number+1.0;
       if(Math.abs(ckN.getValue1()-ck.getValue1())>=2)
        dist+=Math.log10(Math.abs(ckN.getValue1()-ck.getValue1()))/Math.log10(2);
       else dist+=2;
       element=locationN.get(ogS1);
       element.getValue1().set(COGIndex.get(ogS), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS), dist);
       locationN.put(ogS1, element);
       
         }
       }

       }
      }
    }
    
    
     void computeLogSimilaritynonCH(Gene gene, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation, int taxID){  
     
    for(int z=0;z<gene.genesContig.size();z++){
        gene.genes=gene.genesContig.get(z);
    
      for(int i=0;i<gene.genes.size();i++){
       Triplet<Integer,Integer,String> ck=gene.genes.get(i);
       //System.out.println("Ime COG-a: "+ck.getValue2());
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ck.getValue2())))
           continue;
       
      // HashSet<Pair<String,String>> ogs=geneOGMap.get(ck.getValue2());
       
       HashSet<Pair<String,String>> OGsT=geneOGMap.get(ck.getValue2());
       
        HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
       
       /*if(!(cgmap.CogGOmap.keySet().contains(ck.getValue2())))
           continue;*/
       for(int j=i+1;j<gene.genes.size();j++){
          Triplet<Integer,Integer,String> ckN=gene.genes.get(j);
          
          
          if(ckN.getValue2().equals("-"))
           continue;
       if(!(geneOGMap.containsKey(ckN.getValue2())))
           continue;
       /* if(ck.getValue2().equals(ckN.getValue2()))
           continue;*/
     
        //HashSet<Pair<String,String>> ogs1=geneOGMap.get(ckN.getValue2());
          
         HashSet<Pair<String,String>> OGsT1=geneOGMap.get(ckN.getValue2());
     
        
         HashSet<String> OGs1 = new HashSet<>();
            
            for(Pair<String,String> p:OGsT1){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs1.add(p.getValue0());
            }
        
             /*if(ck.getValue2().equals("FOXG_05859") || (ck.getValue2().equals("FVEG_03736T0")) || (ck.getValue2().equals("FVEG_03736")) || (ck.getValue2().equals("FOXG_05859P0")) || (ck.getValue2().equals("NechaG86793")) || (ck.getValue2().equals("NechaG82508")) || (ck.getValue2().equals("NechaP86793")) || (ck.getValue2().equals("NechaP82508"))){
                 System.out.println(ck.getValue2()+" "+OGs.size()+" "+OGs1.size());
            }*/
             
             
             /*if(OGs.contains("ENOG410PSIY")){
                 System.out.println("taxId: "+taxID);
                 System.out.println("gene: "+ck.getValue2());
                 System.out.println("neighbour gene: "+ckN.getValue2());
                 System.out.println("NOgs: ");
                 for(String og:OGs1)
                     System.out.print(og+" ");
                 System.out.println();
             }*/
            
          /*if(ckN.getValue2().equals("-"))
           continue;
       if(!(cgmap.CogGOmap.keySet().contains(ckN.getValue2())))
           continue;
        if(ck.getValue2().equals(ckN.getValue2()))
          continue;*/
     
         for(String ogS:OGs){ 
          if(!COGIndex.containsKey(ogS))
              continue;
          for(String ogS1:OGs1){
             if(!COGIndex.containsKey(ogS1))
                 continue;
       double dist=locationN.get(ogS).getValue0().get(COGIndex.get(ogS1));
       double number=locationN.get(ogS).getValue1().get(COGIndex.get(ogS1));
       number=number+1.0;
       /*if(Math.abs(ckN.getValue1()-ck.getValue1())==0){
           System.out.println("Something wrong, dist 0");
           System.out.println(ckN.getValue1()+" | "+ck.getValue1());
           System.out.println(ckN.getValue2()+" | "+ck.getValue2());
       }*/
       if(Math.abs(ckN.getValue1()-ck.getValue1())>=2)
        dist+=Math.log10(Math.abs(ckN.getValue1()-ck.getValue1()))/Math.log10(2);
       else dist+=2;
       //dist+=Math.log10(distance.computeDistance(ck.getValue0(), ck.getValue1(), ckN.getValue0(), ckN.getValue1(),cog.maxCoordinate))/Math.log10(2);//ln vs log2
       Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(ogS);
       element.getValue1().set(COGIndex.get(ogS1), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS1), dist);
       locationN.put(ogS, element);
       
         dist=locationN.get(ogS1).getValue0().get(COGIndex.get(ogS));
       number=locationN.get(ogS1).getValue1().get(COGIndex.get(ogS));
       number=number+1.0;
       if(Math.abs(ckN.getValue1()-ck.getValue1())>=2)
        dist+=Math.log10(Math.abs(ckN.getValue1()-ck.getValue1()))/Math.log10(2);
       else dist+=2;
       element=locationN.get(ogS1);
       element.getValue1().set(COGIndex.get(ogS), number);//.setAt1(number);
       element.getValue0().set(COGIndex.get(ogS), dist);
       locationN.put(ogS1, element);
       
         }
       }

       }
      }
    }
  }
    
  /*  void computeLogSimilarityLoc(Gene gene, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, int taxId){  
        
      for(int i=0;i<gene.genes.size();i++){
       Triplet<Integer,Integer,String> ck=gene.genes.get(i);
       //System.out.println("Ime COG-a: "+ck.getValue2());
       
       if(ck.getValue2().equals("-"))
           continue;
      
        if(!(geneOGMap.containsKey(ck.getValue2())))
           continue;
       
       HashSet<Pair<String,String>> ogs=geneOGMap.get(ck.getValue2());
       
       
       
       for(int j=0;j<gene.genes.size();j++){
          Triplet<Integer,Integer,String> ckN=gene.genes.get(j);
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
    }*/
    
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
    
    
    void normalize1(int type){
        Iterator<String> cogIt=locationN.keySet().iterator();
        System.out.println("locationN size: "+locationN.keySet().size());
        int count=0;
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            if(count==0){
            System.out.println("element size: "+element.getValue0().size());
            count=1;
            }
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(element.getValue1().get(i)!=0 && i!=COGIndex.get(cog)){
                 dist/=element.getValue1().get(i);
                // System.out.println("dist: "+dist);
                 element.getValue0().set(i, dist);
                 }
                 else if(i!=COGIndex.get(cog)){
                     if(dist!=0)
                         System.out.println("dist: "+dist);
                     element.getValue0().set(i, 0.0);
                 //  element.getValue0().set(i,Double.POSITIVE_INFINITY);  
                 }
                 else if(i==COGIndex.get(cog) && (type==0 || type>1))
                     element.getValue0().set(i, 0.0);
                 else if(i==COGIndex.get(cog) && type==1)
                     element.getValue0().set(i, 0.0/*Math.log(Math.pow(10.0, -10.0))*/);                     
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
        System.out.println("locationN size: "+locationN.keySet().size());
        int count=0;
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            if(count==0){
            System.out.println("element size: "+element.getValue0().size());
            count=1;
            }
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 if(element.getValue1().get(i)!=0 && i!=COGIndex.get(cog)){
                 dist/=element.getValue1().get(i);
                // System.out.println("dist: "+dist);
                 element.getValue0().set(i, dist);
                 }
                 else if(i!=COGIndex.get(cog)){
                     if(dist!=0)
                         System.out.println("dist: "+dist);
                   element.getValue0().set(i,Double.POSITIVE_INFINITY);  
                 }
                 else if(i==COGIndex.get(cog) && (type==0 || type>1))
                     element.getValue0().set(i, 0.0);
                 else if(i==COGIndex.get(cog) && type==1)
                     element.getValue0().set(i, 0.0/*Math.log(Math.pow(10.0, -10.0))*/);                     
            }
             locationN.put(cog, element);
        }
    }
     
     void removeEmpty(){
         
          Iterator<String> cogIt=locationN.keySet().iterator();
          ArrayList<Integer> removed = new ArrayList<>();
          HashSet<String> removedOG = new HashSet<>();
          
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            int numOK=0;
            
            for(int i=0;i<element.getValue0().size();i++){
                 double dist=element.getValue0().get(i);
                 
                 if(Math.abs(dist)<=Double.MAX_VALUE)
                     numOK++;
            }
            if(numOK<=1){
                removedOG.add(cog);
              //  System.out.println("removed OG-ind: "+cog+" "+COGIndex.get(cog));
             removed.add(COGIndex.get(cog));
             cogIt.remove();
            }
        }
       
       Collections.sort(removed, Collections.reverseOrder());
        
        cogIt=locationN.keySet().iterator();
        
        while(cogIt.hasNext()){
            String cog=cogIt.next();
            Pair<ArrayList<Double>,ArrayList<Double>> element=locationN.get(cog);
            int numOK=0;
            
            for(int i=0;i<removed.size();i++){
                ArrayList<Double> tmp=element.getValue0();//.remove(removed.get(i).intValue());
                tmp.remove(removed.get(i).intValue());
                element=element.setAt0(tmp);
                tmp=element.getValue1();//.remove(removed.get(i).intValue());
                tmp.remove(removed.get(i).intValue());
                element=element.setAt1(tmp);
                //System.out.println("removed index: "+removed.get(i));
            }   
            locationN.put(cog, element);
        }
        
    //HashMap<String,Integer> COGIndex=new HashMap<>();
    //HashMap<Integer,String> IndexCOG=new HashMap<>();
      
        HashMap<String,Integer> COGIndexT=new HashMap<>();
        HashMap<Integer,String> IndexCOGT=new HashMap<>();
        
        for(int i=0;i<removed.size();i++)
            IndexCOG.remove(removed.get(i));
        
        for(String s:removedOG)
            COGIndex.remove(s);
        
        int tc=0, maxInd=-1;   
        
        for(Integer i:IndexCOG.keySet()){
            String cog = IndexCOG.get(i);
            if(i>maxInd)
                maxInd=i;
        }
        
         for(int i=0;i<maxInd+1;i++){
             if(!IndexCOG.containsKey(i))
                 continue;
            String cog = IndexCOG.get(i);
            IndexCOGT.put(tc, cog);
            COGIndexT.put(cog,tc++);
         }
        
        IndexCOG.clear();
        COGIndex.clear();
        
        IndexCOG=IndexCOGT;
        COGIndex=COGIndexT;
        
     }
}

