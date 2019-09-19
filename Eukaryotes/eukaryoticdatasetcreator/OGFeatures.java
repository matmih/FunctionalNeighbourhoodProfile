/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store OG features
 */
public class OGFeatures {
    /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

    HashMap<String,Pair<ArrayList<Double>,Double>> COGbaselinefeaturesmap=new HashMap<String,Pair<ArrayList<Double>,Double>>();

     void createFeatures(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<String>>  geneOGMap, Mappings map){//add mappings

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            HashSet<String> OGs=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2()))
                 OGs=geneOGMap.get(key.getValue2());
            
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                }
            }
         }     

        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            HashSet<String> neighOGs=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2()))
                 neighOGs.addAll(geneOGMap.get(value.get(i).getValue2()));
           
      for(String og:neighOGs){     
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
     
     //guava
      void createFeaturesNew(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, int taxid){//add mappings

          System.out.println("Computing features!");
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        System.out.println("Num genes in neigh: "+ng.neighbours.keySet().size());
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
                    
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                }
            }
         }  
         
         System.out.println("OGs added!");

        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            
            HashSet<String> neighOGs=new HashSet<>();
            HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
            
            if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid))){
                      if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
             
      for(String og:neighOGs){     
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
        System.out.println("Neighbourhood computed!");
     }
          System.out.println("Organism completed!");
    }
     
     //TODO
     //need taxID ->COG mapping
     void createFeaturesOrganism(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<String>>  geneOGMap, Mappings map){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<String> OGs=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2()))
                 OGs=geneOGMap.get(key.getValue2());
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
            }
         }

        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<String> neighOGs=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2()))
                 neighOGs.addAll(geneOGMap.get(value.get(i).getValue2()));
           
      for(String og:neighOGs){  
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
     
     
      void createFeaturesOrganismNew(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, int taxID){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();//not exactly OK, should be divided  by contigs
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
            }
         }
          
        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
             
             HashSet<String> neighOGs = new HashSet<>();
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
            
      for(String og:neighOGs){  
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
      
       void createFeaturesNewCount(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, HashMap<Integer,HashSet<String>> ogOcc ,int taxid){//add mappings

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
                       
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
                }
            }
         }  
         
        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            //add map also
            
            HashSet<String> neighOGs=new HashSet<>();
            HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
            
            if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid))){
                      if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
             
      for(String og:neighOGs){     
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
          System.out.println("Organism completed!");
    }
       
       
        void createFeaturesNewCountLeaveOut(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, HashMap<Integer,HashSet<String>> ogOcc, HashSet<String> leaveoutOGs, int taxid){//add mappings

          System.out.println("Computing features!");
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        System.out.println("Num genes in neigh: "+ng.neighbours.keySet().size());
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
                }
            }
         }  
         
         System.out.println("OGs added!");

        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            
            HashSet<String> neighOGs=new HashSet<>();
            HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
            
            if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid))){
                      if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }

      for(String og:neighOGs){
          if(leaveoutOGs.contains(og))
              continue;
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
        System.out.println("Neighbourhood computed!");
     }
          System.out.println("Organism completed!");
    }
     
      
      void createFeaturesOrganismNewCount(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation,HashMap<Integer, HashSet<String>> ogOcc , int taxID){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();//not exactly OK, should be divided  by contigs
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                    ogOcc.get(taxID).add(cog);                         
                }
            }
            else if(!ogOcc.get(taxID).contains(og)){
                 String cog=new String(og);
                 Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(cog);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(cog,pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                }
                ogOcc.get(taxID).add(og);
            }
         }
           
        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
             
             HashSet<String> neighOGs = new HashSet<>();
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
             
      for(String og:neighOGs){  
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
      
      
       void createFeaturesOrganismNewCountLeaveOut(GeneNeighbours ng, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation,HashMap<Integer, HashSet<String>> ogOcc , HashSet<String> leaveoutOGs ,int taxID){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();//not exactly OK, should be divided  by contigs
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                    ogOcc.get(taxID).add(cog);                         
                }
            }
            else if(!ogOcc.get(taxID).contains(og)){
                 String cog=new String(og);
                 Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(cog);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(cog,pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                }
                ogOcc.get(taxID).add(og);
            }
         }

        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
             
             HashSet<String> neighOGs = new HashSet<>();
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }

      for(String og:neighOGs){  
          if(leaveoutOGs.contains(og))
              continue;
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
     
         void createFeaturesNewCountNonCH(ArrayList<GeneNeighbours> ngList, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, HashMap<Integer,HashSet<String>> ogOcc ,int taxid){//add mappings

        for(int ngInd =0; ngInd<ngList.size();ngInd++){
            GeneNeighbours ng = ngList.get(ngInd);
        
          System.out.println("Computing features!");
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        System.out.println("Num genes in neigh: "+ng.neighbours.keySet().size());
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
                }
            }
         }  
         
         System.out.println("OGs added!");

        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            //add map also
            
            HashSet<String> neighOGs=new HashSet<>();
            HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
            
            if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid))){
                      if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
             
      for(String og:neighOGs){     
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
        System.out.println("Neighbourhood computed!");
     }
          System.out.println("Organism completed!");
        }
    }
   
         void createFeaturesNewCountNonCHLeaveOut(ArrayList<GeneNeighbours> ngList, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation, HashMap<Integer,HashSet<String>> ogOcc, HashSet<String> leaveoutOGs, int taxid){//add mappings

        for(int ngInd =0; ngInd<ngList.size();ngInd++){
            GeneNeighbours ng = ngList.get(ngInd);
        
          System.out.println("Computing features!");
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        System.out.println("Num genes in neigh: "+ng.neighbours.keySet().size());
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
         for(String og:OGs){
              if(!cgmap.CogGOmap.keySet().contains(og)){
                continue;
            }
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
            }
            else{
                if(!usedOGs.contains(og)){
                   
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(og);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(og,pr);
                usedOGs.add(og);
                if(!ogOcc.containsKey(taxid)){
                    ogOcc.put(taxid, new HashSet<String>());
                }
                ogOcc.get(taxid).add(og); 
                }
            }
         }  
         
         System.out.println("OGs added!");

        ArrayList<Triplet<Integer,Integer,String>> value=ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            //add map also
            
            HashSet<String> neighOGs=new HashSet<>();
            HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
            
            if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid))){
                      if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
             
      for(String og:neighOGs){  
          if(leaveoutOGs.contains(og))
              continue;
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
     
         for(String ogR:OGs){
             if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
        for(int j=0;j<value1.size();j++){
            if(funcs.contains(value1.get(j)))
                continue;
            double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
     }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
        System.out.println("Neighbourhood computed!");
     }
          System.out.println("Organism completed!");
        }
    }
      
      void createFeaturesOrganismNewCountnonCH(ArrayList<GeneNeighbours> ngList, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation,HashMap<Integer, HashSet<String>> ogOcc , int taxID){

          for(int ccInd = 0; ccInd<ngList.size();ccInd++){
              
          GeneNeighbours ng = ngList.get(ccInd);
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();//not exactly OK, should be divided  by contigs
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                    ogOcc.get(taxID).add(cog);                         
                }
            }
            else if(!ogOcc.get(taxID).contains(og)){
                 String cog=new String(og);
                 Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(cog);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(cog,pr);
                ogOcc.get(taxID).add(og);
            }
         }

        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
             
             HashSet<String> neighOGs = new HashSet<>();
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
            
      for(String og:neighOGs){  
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
          }
    }
      
       void createFeaturesOrganismNewCountnonCHLeaveOut(ArrayList<GeneNeighbours> ngList, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>>  geneOGMap, Mappings map, HashMap<Integer,Integer> taxTranslation,HashMap<Integer, HashSet<String>> ogOcc , HashSet<String> leaveoutOGs ,int taxID){

          for(int ccInd = 0; ccInd<ngList.size();ccInd++){
              
          GeneNeighbours ng = ngList.get(ccInd);
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = ng.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

             HashSet<Pair<String,String>> OGsT=new HashSet<>();
            if(geneOGMap.containsKey(key.getValue2())){
                 OGsT=geneOGMap.get(key.getValue2());
            }
            
            HashSet<String> OGs = new HashSet<>();//not exactly OK, should be divided  by contigs
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID)))
                    OGs.add(p.getValue0());
            }
           
         for(String og:OGs){ 
            if(!cgmap.CogGOmap.keySet().contains(og))
                continue;
            if(!COGbaselinefeaturesmap.containsKey(og)){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(og);
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
                if(!ogOcc.containsKey(taxID)){
                    ogOcc.put(taxID, new HashSet<String>());
                    ogOcc.get(taxID).add(cog);                         
                }
            }
            else if(!ogOcc.get(taxID).contains(og)){
                 String cog=new String(og);
                 Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(cog);
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(cog,pr);
                ogOcc.get(taxID).add(og);
            }
         }

        ArrayList<Triplet<Integer,Integer,String>> value= ng.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

             HashSet<Pair<String,String>> neighOGsT=new HashSet<>();
           
             if(geneOGMap.containsKey(value.get(i).getValue2())){
                 neighOGsT.addAll(geneOGMap.get(value.get(i).getValue2()));
             }
             
             HashSet<String> neighOGs = new HashSet<>();
            
            for(Pair<String,String> p:neighOGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    if(!OGs.contains(p.getValue0()))
                          neighOGs.add(p.getValue0());
                }
            }
            
      for(String og:neighOGs){  
          if(leaveoutOGs.contains(og))
              continue;
        ArrayList<String> value1=cgmap.CogGOmap.get(og);
        if(!cgmap.CogGOmap.keySet().contains(og)){
            continue;
        }
        for(String ogR:OGs){
            if(!cgmap.CogGOmap.keySet().contains(ogR)){
                continue;
             }
             for(int j=0;j<value1.size();j++){
                 if(funcs.contains(value1.get(j)))
                        continue;
                double elem=COGbaselinefeaturesmap.get(ogR).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
                COGbaselinefeaturesmap.get(ogR).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
            }
        }
        funcs.addAll(value1);
       }
      funcs.clear();
      }
     }
    }
  }
      
      
     void normalize(){
         
          Iterator<String> keySetIterator = COGbaselinefeaturesmap.keySet().iterator();
         
          while(keySetIterator.hasNext()){
              String COG=keySetIterator.next();
              Pair<ArrayList<Double>, Double> pr=COGbaselinefeaturesmap.get(COG);
              ArrayList<Double> tmp=pr.getValue0();
              
              for(int i=0;i<tmp.size();i++){
                  tmp.set(i, (tmp.get(i)/pr.getValue1()));
              }
          }
     }
}