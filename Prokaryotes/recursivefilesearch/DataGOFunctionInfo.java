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
import java.util.Random;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute similarity information about pairs of functions (Jaccard, F-score)
 */
public class DataGOFunctionInfo {   
         public HashMap<Integer,ArrayList<Integer>> datasetOccurence=new HashMap<>();

         void computeOccurence(COGGOMap cgmap){
            Iterator<String> it=cgmap.CogGOmap.keySet().iterator();
            
            HashMap<String,Integer> cogMapping=new HashMap<>();
           int count=0;
            while(it.hasNext()){
                    cogMapping.put(it.next(), count++);
           }
            it=cgmap.CogGOmap.keySet().iterator();
            
            while(it.hasNext()){
                String cog=it.next();
                
                ArrayList<String> gos=cgmap.CogGOmap.get(cog);
                
                for(int i=0;i<gos.size();i++){
                    if(gos.get(i).equals("GO0051641") || gos.get(i).equals("GO0006928"))
                        System.out.println(cog+" "+gos.get(i));
                    int go=cgmap.GOtoIndex.get(gos.get(i));
                    if(!datasetOccurence.containsKey(go)){
                        datasetOccurence.put(go, new ArrayList<Integer>(Collections.nCopies(cgmap.CogGOmap.keySet().size(),0)));
                        datasetOccurence.get(go).set(cogMapping.get(cog), 1);
                    }
                    else{
                        datasetOccurence.get(go).set(cogMapping.get(cog), 1);
                    }   
                }
            }
         }
         
         void randomizedOccurence(COGGOMap cgmap,DataGOFunctionInfo datInfo){
             
             Iterator<Integer> itmap = datInfo.datasetOccurence.keySet().iterator();
            
            Random rand = new Random();
            
             while(itmap.hasNext()){
                int go=itmap.next();
                
                ArrayList<Integer> cogs=datInfo.datasetOccurence.get(go);
                datasetOccurence.put(go, new ArrayList<Integer>(Collections.nCopies(cgmap.CogGOmap.keySet().size(),0)));
                
                int count=0;
                
                for(int j=0;j<cogs.size();j++){
                    if(cogs.get(j)==1)
                        count++;
                }
                
                  int randomNum; 
                  
                  HashSet<Integer> locs = new HashSet();
                  ArrayList<Integer> distLoc = new ArrayList<>();
                  
                  for(int j=0;j<count;j++){
                        randomNum= rand.nextInt((cgmap.CogGOmap.keySet().size()-1 - 0) + 1) + 0;
                        while(locs.contains(randomNum)){
                            randomNum= rand.nextInt((cgmap.CogGOmap.keySet().size()-1 - 0) + 1) + 0;
                        }
                        locs.add(randomNum);
                        distLoc.add(randomNum);
                  }
                  
                  for(int j=0;j<distLoc.size();j++){
                      datasetOccurence.get(go).set(distLoc.get(j), 1);
                  }
            }
             
         }
         
         double computeJaccard(int go1, int go2){
             
             ArrayList<Integer> vecgo1=datasetOccurence.get(go1);
             ArrayList<Integer> vecgo2=datasetOccurence.get(go2);
             
             double intersection=0.0;
             double union=0.0;
             
             for(int i=0;i<vecgo1.size();i++){
                 if(vecgo1.get(i)==1 && vecgo2.get(i)==1)
                     intersection=intersection+1.0;
                 if(vecgo1.get(i)==1 || vecgo2.get(i)==1)
                     union=union+1.0;
             }
                          
             return intersection/union;
         }
         
         double computeF1(int go1, int go2){
             double tp=0.0, fn=0.0, fp=0.0,tn=0.0;
             
             ArrayList<Integer> vecgo1=datasetOccurence.get(go1);
             ArrayList<Integer> vecgo2=datasetOccurence.get(go2);
             
             for(int i=0;i<vecgo1.size();i++){
                 if(vecgo1.get(i)==1 && vecgo2.get(i)==1)
                     tp=tp+1.0;
                 else if(vecgo1.get(i)==0 && vecgo2.get(i)==1)
                     fp=fp+1.0;
                 else if(vecgo1.get(i)==1 && vecgo2.get(i)==0)
                     fn=fn+1.0;
                 else if(vecgo1.get(i)==0 && vecgo2.get(i)==0)
                     tn=tn+1.0;
             }
             
             double prec=tp/(tp+fp);
             double rec=tp/(tp+fn);
             if(tp==0)
                 rec=1;
             
             return 2*prec*rec/(prec+rec);   
         }      
}
