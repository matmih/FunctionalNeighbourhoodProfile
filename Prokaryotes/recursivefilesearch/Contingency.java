/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.FileWriter;
import java.io.IOException;
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
  *@description class to store contingency tables
 */
public class Contingency {
   HashMap<Integer,HashSet<Integer>> GOtoCogs=new HashMap<>();
   HashMap<Integer,HashMap<Integer,ArrayList<Integer>>> contTable=new HashMap<>(); 
   HashMap<String, Integer> CogToIndex=new HashMap<>();
   
   void createGOCog(COGGOMap map){
       Iterator<String> iter=map.CogGOmap.keySet().iterator();
       
       int count=0;
       while(iter.hasNext()){
           String cog=iter.next();
           CogToIndex.put(cog, count++);
       }
       iter=map.CogGOmap.keySet().iterator();
       
       while(iter.hasNext()){
           String cog=iter.next();
           ArrayList<String> gos=map.CogGOmap.get(cog);
           
           for(int i=0;i<gos.size();i++){
               if(!GOtoCogs.containsKey(map.GOtoIndex.get(gos.get(i)))){
                   HashSet<Integer> tmp=new HashSet<>();
                   tmp.add(CogToIndex.get(cog));
                   GOtoCogs.put(map.GOtoIndex.get(gos.get(i)), tmp);
               }
               else{
                   HashSet<Integer> tmp=GOtoCogs.get(map.GOtoIndex.get(gos.get(i)));
                   tmp.add(CogToIndex.get(cog));
                   GOtoCogs.put(map.GOtoIndex.get(gos.get(i)), tmp);
               }
           }
               
       }
   }
   
   void initializeContingency(COGGOMap map){

           for(int i:map.IndexToGO.keySet()){
               for(int j:map.IndexToGO.keySet())
                   if(!contTable.containsKey(i)){
                       HashMap<Integer,ArrayList<Integer>> tHash=new HashMap<>();
                       ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                       tHash.put(j, tA);
                       contTable.put(i, tHash);
                   }
                   else{
                      HashMap<Integer,ArrayList<Integer>> tHash=contTable.get(i);
                      ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                      tHash.put(j, tA);
                      contTable.put(i, tHash);
                   }
       }
       
   }
   
    void initializeContingencyPairs(COGGOMap map, ArrayList<Pair<String,String>> GOPairs){

        for(int k=0;k<GOPairs.size();k++){
           //for(int i:map.IndexToGO.keySet()){
            String GO=GOPairs.get(k).getValue0();
            System.out.println("GO1: "+GO);
            int i=map.GOtoIndex.get(GO);
            GO=GOPairs.get(k).getValue1();
            System.out.println("GO2: "+GO);
            int j=map.GOtoIndex.get(GO);
               //for(int j:map.IndexToGO.keySet())
                   if(!contTable.containsKey(i)){
                       HashMap<Integer,ArrayList<Integer>> tHash=new HashMap<>();
                       ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                       tHash.put(j, tA);
                       contTable.put(i, tHash);
                   }
                   else{
                      HashMap<Integer,ArrayList<Integer>> tHash=contTable.get(i);
                      ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                      tHash.put(j, tA);
                      contTable.put(i, tHash);
                   }
       }     
   }
   
     void initializeContingency(COGGOMap map, String COG){

           //for(int i:map.IndexToGO.keySet()){
         int i=map.GOtoIndex.get(COG);
                 {
               for(int j:map.IndexToGO.keySet())
                   if(!contTable.containsKey(i)){
                       HashMap<Integer,ArrayList<Integer>> tHash=new HashMap<>();
                       ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                       tHash.put(j, tA);
                       contTable.put(i, tHash);
                   }
                   else{
                      HashMap<Integer,ArrayList<Integer>> tHash=contTable.get(i);
                      ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                      tHash.put(j, tA);
                      contTable.put(i, tHash);
                   }
       }   
   }
     
     void initializeContingency(COGGOMap map, ArrayList<String> COGs, ArrayList<HashSet<String>> categories){

           for(int i1=0;i1<COGs.size();i1++){
         int i=map.GOtoIndex.get(COGs.get(i1));
               for(int j:map.IndexToGO.keySet()){
                   if(!categories.get(0).contains(map.IndexToGO.get(j)))
                       continue;
                   if(!contTable.containsKey(i)){
                       HashMap<Integer,ArrayList<Integer>> tHash=new HashMap<>();
                       ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                       tHash.put(j, tA);
                       contTable.put(i, tHash);
                   }
                   else{
                      HashMap<Integer,ArrayList<Integer>> tHash=contTable.get(i);
                      ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                      tHash.put(j, tA);
                      contTable.put(i, tHash);
                   }
               }
       }   
   }
     
     void initializeContingency(COGGOMap map, String COG, ArrayList<HashSet<String>> categories){

           //for(int i:map.IndexToGO.keySet()){
         int i=map.GOtoIndex.get(COG);
                 {
               for(int j:map.IndexToGO.keySet()){
                   if(!categories.get(0).contains(map.IndexToGO.get(j)))
                       continue;
                   if(!contTable.containsKey(i)){
                       HashMap<Integer,ArrayList<Integer>> tHash=new HashMap<>();
                       ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                       tHash.put(j, tA);
                       contTable.put(i, tHash);
                   }
                   else{
                      HashMap<Integer,ArrayList<Integer>> tHash=contTable.get(i);
                      ArrayList<Integer> tA=new ArrayList<>(Collections.nCopies(4, 0));
                      tHash.put(j, tA);
                      contTable.put(i, tHash);
                   }
               }          
       }  
   }
   
   void addContingency(COGSimilarity1 sim, COGGOMap cgmap, GeneOGMapping geneOGMap){
       
       Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        
       int totalSize=sim.neighbours.size();
       int complited=0;
       double percent=0.0; 
         long startTime = System.currentTimeMillis();      
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
            
            HashSet<String> OGs=geneOGMap.geneOGsMap.get(key.getValue2());
            
            /*if(!cgmap.CogGOmap.containsKey(key.getValue2()))
                continue;*/
            /*else{
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(key.getValue2());
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
               // System.out.println("tmp: "+tmp);
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(key.getValue2(),pr);
                //COGbaselinefeaturesmap.get(key.getValue2()).setAt1((COGbaselinefeaturesmap.get(key.getValue2()).getValue1()+1.0));
            }*/

    
        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        HashSet<Integer> cogsGOsSet=new HashSet<>();
        
        for(String og:OGs){
                if(!cgmap.CogGOmap.keySet().contains(og))
                        continue;
    
        ArrayList<String> cogsGOs=cgmap.CogGOmap.get(og/*key.getValue2()*/);
        //HashSet<String> cogsGOsSet=new HashSet<>();
        
        for(String s:cogsGOs)
            cogsGOsSet.add(cgmap.GOtoIndex.get(s));
        }
        //HashSet<String> neighbGOs=new HashSet<>();
        HashSet<Integer> neighbGOs=new HashSet<>();
        //for(String cg:cgmap.GOtoIndex.keySet()){
            
        for(int i=0;i<value.size();i++){
            
            HashSet<String> OGsN=geneOGMap.geneOGsMap.get(value.get(i).getValue2());
            
            for(String og:OGsN){
                    if(!cgmap.CogGOmap.keySet().contains(og))
                        continue;
                    
        ArrayList<String> value1=cgmap.CogGOmap.get(og/*value.get(i).getValue2()*/);
        /*if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
          //  System.out.println("Nemamo GO funkcije za COG: "+value.get(i).getValue2());
            continue;
        }*/
        for(int j=0;j<value1.size();j++){
             //if(cogsGOs.contains(value1.get(j)))
              neighbGOs.add(cgmap.GOtoIndex.get(value1.get(j)));
        }
       }
      // }
      }
        //have all GOs in the neighbourhood
        
        for(Integer cg:cgmap.IndexToGO.keySet()){
            for(Integer cg1:cgmap.IndexToGO.keySet())
            if(cogsGOsSet.contains(cg)){
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(/*cgmap.GOtoIndex.get(cg)*/ cg).get(cg1/*cgmap.GOtoIndex.get(cg1)*/).get(0);
                   val=val+1;
                   contTable.get(cg/*cgmap.GOtoIndex.get(cg)*/).get(cg1/*cgmap.GOtoIndex.get(cg1)*/).set(0,val);
                }
                else{
                   int val=contTable.get(cg/*cgmap.GOtoIndex.get(cg)*/).get(cg1/*cgmap.GOtoIndex.get(cg1)*/).get(1);
                   val=val+1;
                   contTable.get(cg/*cgmap.GOtoIndex.get(cg)*/).get(cg1/*cgmap.GOtoIndex.get(cg1)*/).set(1,val); 
                }
            }
            else{
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(/*cgmap.GOtoIndex.get(cg)*/cg).get(/*cgmap.GOtoIndex.get(cg1)*/cg1).get(2);
                   val=val+1;
                   contTable.get(/*cgmap.GOtoIndex.get(cg)*/cg).get(/*cgmap.GOtoIndex.get(cg1)*/cg1).set(2,val);
                }
                else{
                   int val=contTable.get(/*cgmap.GOtoIndex.get(cg)*/cg).get(/*cgmap.GOtoIndex.get(cg1)*/cg1).get(3);
                   val=val+1;
                   contTable.get(/*cgmap.GOtoIndex.get(cg)*/cg).get(/*cgmap.GOtoIndex.get(cg1)*/cg1).set(3,val); 
                }
            }
        }
        
        complited++;
        percent=(((double)complited)/totalSize)*100;
        if(complited%100==0 || percent>99){
        System.out.println("complited: "+complited);
        System.out.println("percent: "+percent+"%");
        }
    }
      long estimatedTime = System.currentTimeMillis() - startTime;
      System.out.println("Time: "+(estimatedTime/1000)+" s");
       
   }
   
    void addContingency(COGSimilarity1 sim, COGGOMap cgmap, String GOIn){
       
       Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        
       int totalSize=sim.neighbours.size();
       int complited=0;
       double percent=0.0; 
               
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            /*if(!cgmap.CogGOmap.containsKey(key.getValue2()))
                continue;*/
            /*else{
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(key.getValue2());
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
               // System.out.println("tmp: "+tmp);
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(key.getValue2(),pr);
                //COGbaselinefeaturesmap.get(key.getValue2()).setAt1((COGbaselinefeaturesmap.get(key.getValue2()).getValue1()+1.0));
            }*/

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        ArrayList<String> cogsGOs=cgmap.CogGOmap.get(key.getValue2());
        HashSet<String> neighbGOs=new HashSet<>();
        //for(String cg:cgmap.GOtoIndex.keySet()){
        
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
          //  System.out.println("Nemamo GO funkcije za COG: "+value.get(i).getValue2());
            continue;
        }
        for(int j=0;j<value1.size();j++){
             //if(cogsGOs.contains(value1.get(j)))
              neighbGOs.add(value1.get(j));
        }
      // }
      }
        //have all GOs in the neighbourhood
        
        //for(String cg:cgmap.GOtoIndex.keySet()){
        String cg=GOIn;
              {
            for(String cg1:cgmap.GOtoIndex.keySet())
            if(cogsGOs.contains(cg)){
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0);
                   val=val+1;
                   //System.out.println("value pos pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(0,val);
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1);
                   val=val+1;
                  // System.out.println("value pos neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(1,val); 
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1));
                }
            }
            else{
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2);
                   val=val+1;
                  // System.out.println("value neg pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(2,val);
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3);
                   val=val+1;
                 //  System.out.println("value neg neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(3,val); 
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3));
                }
            }
        }
        
      /*  complited++;
        percent=(((double)complited)/totalSize)*100;
        System.out.println("complited: "+complited);
        System.out.println("percent: "+percent+"%");*/
    }            
   }
            
     void addContingency(COGSimilarity1 sim, COGGOMap cgmap, String GOIn, ArrayList<HashSet<String>> categories){
       
       Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        
       int totalSize=sim.neighbours.size();
       int complited=0;
       double percent=0.0; 
               
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            /*if(!cgmap.CogGOmap.containsKey(key.getValue2()))
                continue;*/
            /*else{
                Pair<ArrayList<Double>,Double> pr=COGbaselinefeaturesmap.get(key.getValue2());
                double tmp=pr.getValue1();
                tmp=tmp+1.0;
               // System.out.println("tmp: "+tmp);
                pr=pr.setAt1(tmp);
                COGbaselinefeaturesmap.put(key.getValue2(),pr);
                //COGbaselinefeaturesmap.get(key.getValue2()).setAt1((COGbaselinefeaturesmap.get(key.getValue2()).getValue1()+1.0));
            }*/

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        ArrayList<String> cogsGOs=cgmap.CogGOmap.get(key.getValue2());
        HashSet<String> neighbGOs=new HashSet<>();
        //for(String cg:cgmap.GOtoIndex.keySet()){
        
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
          //  System.out.println("Nemamo GO funkcije za COG: "+value.get(i).getValue2());
            continue;
        }
        for(int j=0;j<value1.size();j++){
             //if(cogsGOs.contains(value1.get(j)))
              neighbGOs.add(value1.get(j));
        }
      // }
      }
        //have all GOs in the neighbourhood
        
        //for(String cg:cgmap.GOtoIndex.keySet()){
        String cg=GOIn;
              {
            for(String cg1:cgmap.GOtoIndex.keySet()){
                if(!categories.get(0).contains(cg1))
                    continue;
            if(cogsGOs.contains(cg)){
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0);
                   val=val+1;
                   //System.out.println("value pos pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(0,val);
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1);
                   val=val+1;
                  // System.out.println("value pos neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(1,val); 
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1));
                }
            }
            else{
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2);
                   val=val+1;
                  // System.out.println("value neg pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(2,val);
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3);
                   val=val+1;
                 //  System.out.println("value neg neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(3,val); 
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3));
                }
            }
        }
       }
        
      /*  complited++;
        percent=(((double)complited)/totalSize)*100;
        System.out.println("complited: "+complited);
        System.out.println("percent: "+percent+"%");*/
    }            
   }  
     
     
     void addContingency(COGSimilarity1 sim, COGGOMap cgmap, ArrayList<String> GOs, ArrayList<HashSet<String>> categories){
       
       Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        
       int totalSize=sim.neighbours.size();
       int complited=0;
       double percent=0.0; 
               
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        ArrayList<String> cogsGOs=cgmap.CogGOmap.get(key.getValue2());
        HashSet<String> neighbGOs=new HashSet<>();
        //for(String cg:cgmap.GOtoIndex.keySet()){
        
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
          //  System.out.println("Nemamo GO funkcije za COG: "+value.get(i).getValue2());
            continue;
        }
        for(int j=0;j<value1.size();j++){
             //if(cogsGOs.contains(value1.get(j)))
              neighbGOs.add(value1.get(j));
        }
      // }
      }

       for(int cnum=0;cnum<GOs.size();cnum++){ 
        String cg=GOs.get(cnum);
              {
            for(String cg1:cgmap.GOtoIndex.keySet()){
                if(!categories.get(0).contains(cg1))
                    continue;
            if(cogsGOs.contains(cg)){
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0);
                   val=val+1;
                   //System.out.println("value pos pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(0,val);
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1);
                   val=val+1;
                  // System.out.println("value pos neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(1,val); 
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1));
                }
            }
            else{
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2);
                   val=val+1;
                  // System.out.println("value neg pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(2,val);
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3);
                   val=val+1;
                 //  System.out.println("value neg neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(3,val); 
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3));
                }
            }
        }
       }
      }
    }            
   }        
            
    void addContingencyRandom(COGSimilarity1 sim, COGGOMap cgmap, ArrayList<Pair<String,String>> GOPairs){
       
       Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
        
       int totalSize=sim.neighbours.size();
       int complited=0;
       double percent=0.0; 
               
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        ArrayList<String> cogsGOs=cgmap.CogGOmap.get(key.getValue2());
        HashSet<String> neighbGOs=new HashSet<>();
        
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
            continue;
        }
        for(int j=0;j<value1.size();j++){
              neighbGOs.add(value1.get(j));
        }
      // }
      }
        //have all GOs in the neighbourhood
     for(int k1=0;k1<GOPairs.size();k1++){
        String cg=GOPairs.get(k1).getValue0();//COG;
              {
           // for(String cg1:cgmap.GOtoIndex.keySet())
         String cg1=GOPairs.get(k1).getValue1();
            if(cogsGOs.contains(cg)){
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0);
                   val=val+1;
                   //System.out.println("value pos pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(0,val);
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(0));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1);
                   val=val+1;
                  // System.out.println("value pos neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(1,val); 
                  // System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(1));
                }
            }
            else{
                if(neighbGOs.contains(cg1)){
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2);
                   val=val+1;
                  // System.out.println("value neg pos: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(2,val);
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(2));
                }
                else{
                   int val=contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3);
                   val=val+1;
                 //  System.out.println("value neg neg: "+val);
                   contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).set(3,val); 
                   //System.out.println("Value in the table: "+contTable.get(cgmap.GOtoIndex.get(cg)).get(cgmap.GOtoIndex.get(cg1)).get(3));
                }
            }
        }
    }    
   }
 }
   
   void writeContingency(COGGOMap map, String outPath){
       
       FileWriter fw;
        try
             {
                 //the true will append the new data
                for(int i:contTable.keySet()){
                    
                    String tmpOut=outPath+map.IndexToGO.get(i)+".txt";
                    fw = new FileWriter(tmpOut,true);
                    System.out.println("key in output table: "+map.IndexToGO.get(i));
                        HashMap<Integer,ArrayList<Integer>> tmp=contTable.get(i);
   
                        for(int j:tmp.keySet()){
                            ArrayList<Integer> ct=tmp.get(j);
                             //System.out.println("key2 in output table: "+map.IndexToGO.get(j));
                             fw.write("Table: "+map.IndexToGO.get(i)+"-"+map.IndexToGO.get(j)+"\n\n");
                             for(int k=0;k<ct.size();k++){
                                  //System.out.println("value in table: "+ct.get(k));
                                 fw.write(ct.get(k)+" ");
                                 if(k==1 || k==3)
                                     fw.write("\n");
                             }
                             fw.write("\n");
                        }
                         fw.close();
                  }   
             }
        catch(IOException e){
            e.printStackTrace();
        }
   }
   
   void writeContingency(COGGOMap map, ArrayList<String> GOs ,String outPath){
       
       FileWriter fw;
       
      for(int fint=0;fint<GOs.size();fint++){ 
        try
             {
                 String outPTmp=outPath.substring(0, outPath.length()-4);
                 outPTmp+=GOs.get(fint)+"base.txt";
                fw = new FileWriter(outPTmp,true); //the true will append the new data
       
                //for(int i:contTable.keySet()){
                int i=map.GOtoIndex.get(GOs.get(fint));
                {
                    System.out.println("key in output table: "+map.IndexToGO.get(i));
                        HashMap<Integer,ArrayList<Integer>> tmp=contTable.get(i);
   
                        for(int j:tmp.keySet()){
                            ArrayList<Integer> ct=tmp.get(j);
                             System.out.println("key2 in output table: "+map.IndexToGO.get(j));
                             fw.write("Table: "+map.IndexToGO.get(i)+"-"+map.IndexToGO.get(j)+"\n\n");
                             for(int k=0;k<ct.size();k++){
                                  System.out.println("value in table: "+ct.get(k));
                                 fw.write(ct.get(k)+" ");
                                 if(k==1 || k==3)
                                     fw.write("\n");
                             }
                             fw.write("\n");
                        }
                  }
                fw.close();
             }
        catch(IOException e){
            e.printStackTrace();
        }
     }
  } 
}
