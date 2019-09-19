/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import org.javatuples.Triplet;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store features for the FNP approach
 */
public class COGfeatures {
    HashMap<String,ArrayList<Integer>> COGfeaturesmap=new HashMap<String,ArrayList<Integer>>();
    int NumFunctionalCOGs=0;
    
    void createFeatures(COGSimilarity sim, COGGOMap cgmap){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

        while(keySetIterator.hasNext()){
        Triplet<Integer,Integer,String> key = keySetIterator.next();
         if(cgmap.CogGOmap.keySet().contains(key))
            NumFunctionalCOGs++; 
        }
        
        keySetIterator = sim.neighbours.keySet().iterator();
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            if(!COGfeaturesmap.containsKey(key.getValue2())){
                ArrayList<Integer> tmp=new ArrayList<Integer>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0));
                String cog=new String(key.getValue2());
                COGfeaturesmap.put(cog, tmp);
            }

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
      
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
            continue;
        }
        for(int j=0;j<value1.size();j++){
            int elem=COGfeaturesmap.get(key.getValue2()).get(cgmap.GOtoIndex.get(value1.get(j)));
            COGfeaturesmap.get(key.getValue2()).set(cgmap.GOtoIndex.get(value1.get(j)),++elem);
        }
       }
     }

    }
    
    void createFeatures(COGSimilarity1 sim, COGGOMap cgmap){

        Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

        while(keySetIterator.hasNext()){
        Triplet<Integer,Integer,String> key = keySetIterator.next();
         if(cgmap.CogGOmap.keySet().contains(key))
            NumFunctionalCOGs++; 
        }
        
        keySetIterator = sim.neighbours.keySet().iterator();
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            if(!COGfeaturesmap.containsKey(key.getValue2())){
                ArrayList<Integer> tmp=new ArrayList<Integer>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0));
                String cog=new String(key.getValue2());
                COGfeaturesmap.put(cog, tmp);
            }

        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
      
        for(int i=0;i<value.size();i++){

        ArrayList<String> value1=cgmap.CogGOmap.get(value.get(i).getValue2());
        if(!cgmap.CogGOmap.keySet().contains(value.get(i).getValue2())){
            continue;
        }
        for(int j=0;j<value1.size();j++){
            int elem=COGfeaturesmap.get(key.getValue2()).get(cgmap.GOtoIndex.get(value1.get(j)));
            COGfeaturesmap.get(key.getValue2()).set(cgmap.GOtoIndex.get(value1.get(j)),++elem);
        }
       }
     }

    }

}
