/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

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
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class storing FNP phenotype feature information.
 */
public class FunctionalFeatures {
    HashMap<String,Pair<ArrayList<ArrayList<Double>>,Double>> COGfeatureMap=new HashMap<>();
    
    public FunctionalFeatures(Mappings map){
        Iterator<String> ogIt = map.ogFunctionMapping.keySet().iterator();
        
        while(ogIt.hasNext()){
            ArrayList<ArrayList<Double>> t = new ArrayList<>();
            double c = 0.0;
            Pair<ArrayList<ArrayList<Double>>,Double> p = new Pair<>(t,c);
            COGfeatureMap.put(ogIt.next(),p);
        }
    }
    
    
    void createFeatures(COGSimilarity1 sim, Mappings map){
        
        
        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

        Iterator<String> ogIt = map.ogFunctionMapping.keySet().iterator();
       
        while(ogIt.hasNext()){
            String og = ogIt.next();
            Pair<ArrayList<ArrayList<Double>>,Double> p = COGfeatureMap.get(og);
            ArrayList<ArrayList<Double>> t = p.getValue0();
            t.add(new ArrayList<Double>(Collections.nCopies(map.GOtoIndex.keySet().size(), 0.0)));
            p = p.setAt0(t);
            COGfeatureMap.put(og, p);
        }
        
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
        for(int i=0;i<value.size();i++){

          HashSet<String> ogs = map.geneOgMapping.get(value.get(i).getValue2());
          HashSet<String> functions = new HashSet<>();
             
          for(String ogS: ogs){
          
        HashSet<String> value1=map.ogFunctionMapping.get(ogS);
        
        if(!map.ogFunctionMapping.keySet().contains(ogS)){
            continue;
        }
        for(String s:value1){
            if(!functions.contains(s)){
            double elem=COGfeatureMap.get(ogS).getValue0().get(0).get(map.GOtoIndex.get(s));
            COGfeatureMap.get(ogS).getValue0().get(0).set(map.GOtoIndex.get(s),elem+1);
            }
        }
        functions.addAll(value1);
       }
        }
     }
        
    }
    
}
