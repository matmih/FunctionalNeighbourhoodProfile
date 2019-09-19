/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to create labeled GFP classifier dataset.
 */
 
public class AddLabelsToGFPData {
    
   static public void main(String [] args){
          int OrganismGroup = 0;//0 fungy, 1 metazoa
         
         if(args.length>0)
            OrganismGroup=Integer.parseInt(args[0]);
          
         System.out.println("num args: "+args.length);
         
       // File ogFuncFile=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_1_30F.txt");
        File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        if(OrganismGroup == 1)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        //File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
        OGGOMapping oggomap=new OGGOMapping();
        
        oggomap.createCOGGOMapping(ogFuncFile);
        System.out.println("Num functions: "+oggomap.GOtoIndex.keySet().size());
        System.out.println("Num ogs: "+oggomap.CogGOmap.keySet().size());
         
         File input=null;
         
         if(OrganismGroup==0)
            input=new File("GFPResultsFungi.txt");
         else if(OrganismGroup==1)
              input=new File("GFPResultsMetazoa.txt");
         
         ArrayList<String> cogsOrdered=new ArrayList<>();
         HashMap<String,Pair<Double,Double>> scores=new HashMap<>();
                 
         try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(input)))
        {
            String line;
            int count=0;
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                String tmp[]=line.split(",");
                String goName=oggomap.IndexToGO.get(count);
                Pair p;
                
                    String tmp1[]=tmp[0].split(":");
                    Double score=Double.parseDouble(tmp1[1]);
                    tmp1=tmp[1].split(":");
                    Double score1=Double.parseDouble(tmp1[1]);
                    p=new Pair(score,score1);
                    scores.put(goName, p);
                    
                count++;
            }
            bufRdr1.close();
        }
       catch(Exception e){
           e.printStackTrace();
       }
         
         File output=null;
         
         if(OrganismGroup==0)
            output=new File("GFPResultsLabFungi.txt");
         else if(OrganismGroup == 1)
            output=new File("GFPResultsLabMetazoa.txt"); 
         
         try
            {
                FileWriter fw = new FileWriter(output);
                
                for(String go:scores.keySet()){
                    Pair p=scores.get(go);
                    
                    fw.write(go+"\t"+p.getValue0()+"\t"+p.getValue1()+"\n");
                }
                   
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            }  
   } 
}
