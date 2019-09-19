 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

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
  *@description class to acount gene arrangements
 */
 
public class CountGeneArrangements {
    HashMap<Integer,ArrayList<String>> genesCount=new HashMap<>();
    
    public CountGeneArrangements(){
        genesCount=new HashMap<>();
    }
    
    public void computeArrangements(int taxId, COGSimilarity1 sim, COGGOMap cgmap, GeneOGMapping geneOGMap, String OGTarget, String GOx, String GOy, COGOccurence cogOc){
        
        
        if(!genesCount.containsKey(taxId)){
            ArrayList<String> orgAr = new ArrayList<>();
           genesCount.put(taxId,orgAr);
        }
        
            Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();

            System.out.println("Obradujem dokument!");
        
          while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();

            HashSet<String> OGs=geneOGMap.geneOGsMap.get(key.getValue2());
           
         if(!OGs.contains(OGTarget))
             continue;
         
        ArrayList<Triplet<Integer,Integer,String>> value=sim.neighbours.get(key);
      
        String order = OGTarget+"_C:GOx"; //"C"+key.getValue2()+"C";
        
        if(cgmap.CogGOmap.get(OGTarget).contains(GOy))
            order+=":GOy";
        
        for(int i=0;i<value.size();i++){
            
            HashSet<String> ogsN = geneOGMap.geneOGsMap.get(value.get(i).getValue2());
            
            if(ogsN == null){
                if(i%2==0 && i==0)
                order+="|"+"NA";//value.get(i).getValue2();
            else if((i%2)==0)
             order+="|"+"NA";//value.get(i).getValue2();
             else if((i%2)==1)
                 order = /*value.get(i).getValue2()*/"NA"+"|"+order;
                continue;
            }  
            else if(ogsN.isEmpty()){
                if(i%2==0 && i==0)
                order+="|"+"NA";//value.get(i).getValue2();
            else if((i%2)==0)
             order+="|"+"NA";//value.get(i).getValue2();
             else if((i%2)==1)
                 order = /*value.get(i).getValue2()*/"NA"+"|"+order;
                continue;
            }
            else{
            String ogCN = cogOc.findMinimumOccuring(ogsN);
            String ocgTemp = ogCN;
            if(cgmap.CogGOmap.get(ocgTemp).contains(GOx))
                ogCN+=":GOx";
            System.out.println("ocgTemp: "+ocgTemp);
            if(cgmap.CogGOmap.get(ocgTemp).contains(GOy))
                ogCN+=":GOy";
          
            if(i%2==0 && i==0)
                order+="|"+ogCN;//value.get(i).getValue2();
            else if((i%2)==0)
             order+="|"+ogCN;//value.get(i).getValue2();
             else if((i%2)==1)
                 order = /*value.get(i).getValue2()*/ogCN+"|"+order;
            }
        }
        
        System.out.println("Order: "+order);
        genesCount.get(taxId).add(order);
       }
    }
    
    public ArrayList<ArrayList<String>> getTopKArrangements(int k){
        ArrayList<ArrayList<String>> ar = new ArrayList<>();
        Iterator<Integer> it= genesCount.keySet().iterator();
        
        while(it.hasNext()){
            int tax = it.next();

                HashSet<String> checked = new HashSet<>();            
           ar.add(genesCount.get(tax));
                
        }
        return ar;
    }
    
}
