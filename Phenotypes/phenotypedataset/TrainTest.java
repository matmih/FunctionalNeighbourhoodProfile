/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create dataset for phenotype prediction
 */
public class TrainTest {
    
    FileWriter fw;
    
  public  TrainTest(String path){
      try{
        fw = new FileWriter(path);
      }
      catch(Exception e){
          e.printStackTrace();
      }
    }
    
   public void createTrainHeader(PhenotypesAndStrains selection,Mappings map){
         try{
 
            int numIter = 1;

            
            fw.write("@RELATION   phenDataFN\n\n");            
            fw.write("@ATTRIBUTE   ID"+"    "+"string"+"\n");
                     
            int numBact = selection.strains.size();
            
           for(int i=0;i<map.cogIndex.keySet().size();i++) 
            for(int j=0;j<map.indexToGO.size();j++){
              
                fw.write("@ATTRIBUTE   "+map.indexToGO.get(j)+"cog"+(map.cogIndex.get(i))+"    "+"numeric"+"\n");
            }
                        
             for(int i=0;i<selection.phenotypes.size();i++){
                 if(i+1<selection.phenotypes.size())
                      fw.write("@ATTRIBUTE   pheno"+(i+1)+"   {1,0}\n");
                 else fw.write("@ATTRIBUTE   pheno"+(i+1)+"   {1,0}\n\n");
             }
             
            fw.write("@DATA\n");
    
   }catch (Exception e) {
            e.printStackTrace();
        }
    }
   
   public void writeBact(String bact){
       try{
             fw.write(bact+",");
       }
       catch(IOException e){
           e.printStackTrace();
       }
   }
    
    public void createTrainBody(FunctionalFeatures cbf, HashMap<String, ArrayList<Integer>> bacterialTargets, Mappings map, String organism){
        
        try{
            Iterator<String> bactIt = cbf.COGfeatureMap.keySet().iterator();
                        
            System.out.println("Num cogs mapping: "+map.cogIndex.size());
             System.out.println("Num cogs keyset: "+map.cogIndex.keySet().size());
            
            for(int iC=0;iC<map.cogIndex.size();iC++){
                
                String og = map.cogIndex.get(iC);
                                
                ArrayList<ArrayList<Double>> f = cbf.COGfeatureMap.get(og).getValue0();//write targets
               System.out.println("f size: "+f.size());
          
                for(int i=0;i<f.size();i++){
                     System.out.print("f"+i+" size: "+f.get(i).size()+" ");
                    for(int j=0;j<f.get(i).size();j++){
                       
                                fw.write(f.get(i).get(j)+",");
                    }
                    System.out.println();
                }
            }
          
            ArrayList<Integer> targets = bacterialTargets.get(organism);
            
            for(int i=0;i<targets.size();i++){
                if((i+1)<targets.size())
                     fw.write(targets.get(i)+",");
                else fw.write(targets.get(i)+"\n");
            }
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    
    public void Close(){
        try{
            fw.close();
        
        }
         catch(IOException e){
            e.printStackTrace();
        }
    }    
}
