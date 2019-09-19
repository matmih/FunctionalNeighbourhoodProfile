/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;

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
import static recursivefilesearch.COG.ENCODING;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store features for the FNP approach
 */
public class COGbaselinefeatures {
    HashMap<String,Pair<ArrayList<Double>,Double>> COGbaselinefeaturesmap=new HashMap<String,Pair<ArrayList<Double>,Double>>();

     void createFeatures(COGSimilarity sim, COGGOMap cgmap){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            if(!COGbaselinefeaturesmap.containsKey(key.getValue2())){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
                String cog=new String(key.getValue2());
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
            }
            else{
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(key.getValue2());
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(key.getValue2(),pr);
            }
                

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
      
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
            continue;
        }
        for(int j=0;j<value1.size();j++){
            double elem=COGbaselinefeaturesmap.get(key.getValue2()).getValue0().get(cgmap.GOtoIndex.get(value1.get(j)));
            COGbaselinefeaturesmap.get(key.getValue2()).getValue0().set(cgmap.GOtoIndex.get(value1.get(j)),elem+1);
        }
       }
     }

    }
     
     
     void createFeatures(COGSimilarity1 sim, COGGOMap cgmap, GeneOGMapping geneOGMap){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        HashSet<String> usedOGs=new HashSet<>();
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            HashSet<String> OGs=geneOGMap.geneOGsMap.get(key.getValue2());
           
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

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){
            
           HashSet<String> neighOGs = geneOGMap.geneOGsMap.get(value.get(i).getValue2());

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
     
     void loadFeatures(String input, COGGOMap cgmap){
         BufferedReader reader;
         try {
      Path path =Paths.get(input);
       reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          String[] feat=line.split(" ");
          if(!COGbaselinefeaturesmap.containsKey(feat[0])){
                ArrayList<Double> tmp=new ArrayList<Double>(Collections.nCopies(feat.length-1, 0.0));
                String cog=new String(feat[0]);
                for(int i=0;i<tmp.size();i++)
                    tmp.set(i, Double.parseDouble(feat[i+1]));
                Pair<ArrayList<Double>,Double> pr=new Pair(tmp,1.0);
                COGbaselinefeaturesmap.put(cog, pr);
            }
      }
      }
         catch(Exception e){
             e.printStackTrace();
         }
     }

     void createFeaturesOrganism(COGSimilarity1 sim, COGGOMap cgmap, GeneOGMapping geneOGMap){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

            System.out.println("Obradujem kromosom!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            HashSet<String> OGs=geneOGMap.geneOGsMap.get(key.getValue2());
           
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

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        HashSet<String> funcs=new HashSet<>();
      
        for(int i=0;i<value.size();i++){

            HashSet<String> neighOGs= geneOGMap.geneOGsMap.get(value.get(i).getValue2());

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
