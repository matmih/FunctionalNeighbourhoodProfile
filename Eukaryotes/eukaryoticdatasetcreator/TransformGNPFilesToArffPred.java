/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to transform GFP files to Arff prediction files for eukaryotic organisms
 */
public class TransformGNPFilesToArffPred {
    
    public static void main(String [] args){
         
        int OrganismGroup = 0;//0 fungy, 1 metazoa
         
         if(args.length>0)
            OrganismGroup=Integer.parseInt(args[0]);
          
         System.out.println("num args: "+args.length);
         
        File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        if(OrganismGroup == 1)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        OGGOMapping oggomap=new OGGOMapping();
        
        oggomap.createCOGGOMapping(ogFuncFile);
        System.out.println("Num functions: "+oggomap.GOtoIndex.keySet().size());
        System.out.println("Num ogs: "+oggomap.CogGOmap.keySet().size());
      
        HashMap<String,ArrayList<Double>> realLabels=new HashMap<>();
        HashMap<String,ArrayList<Double>> predictedScores=new HashMap<>();
        ArrayList<String> cogOrder=new ArrayList<>();
        
        FileReader r = null;
        
       try{ 
        if(OrganismGroup==0)
            r=new FileReader("AssociationTargetsFungi.txt");
        else if(OrganismGroup==1)
            r = new FileReader("AssociationTargetsMetazoa.txt");
       }
       catch(IOException e){
           e.printStackTrace();
       }
        
         try (BufferedReader bufRdr = new BufferedReader(r))
        {
            String line;

            while ((line = bufRdr.readLine()) != null)
            {
                String tmp[]=line.split("\t");
                String cog=tmp[0].trim();
                cogOrder.add(cog);
                
                realLabels.put(cog, new ArrayList<Double>());
                
                for(int i=1;i<tmp.length;i++){
                    double lab=Double.parseDouble(tmp[i]);
                    realLabels.get(cog).add(lab);
                }
                
            }
            bufRdr.close();
        }
        catch(Exception e){e.printStackTrace();}
        
         FileReader  rb = null;
         
         try{
             if(OrganismGroup == 0)
                 rb = new FileReader("GFPScoreMatrixFungi.txt");
             else if(OrganismGroup == 1)
                 rb = new FileReader("GFPScoreMatrixMetazoa.txt");
         }catch(IOException e){          
         }
         
        try (BufferedReader bufRdr = new BufferedReader(rb))
        {
            String line;
            int count=0;
            while ((line = bufRdr.readLine()) != null)
            {
                String tmp[]=line.split("\t");
                String go=oggomap.IndexToGO.get(count);
                
                for(int i=0;i<tmp.length;i++){
                    if(!predictedScores.containsKey(cogOrder.get(i)))
                             predictedScores.put(cogOrder.get(i), new ArrayList<Double>(Collections.nCopies(oggomap.GOtoIndex.keySet().size(), -1.0)));
                    double lab=Double.parseDouble(tmp[i]);
                    predictedScores.get(cogOrder.get(i)).set(count,lab);
                }
              count++;  
            }
            bufRdr.close();
        }
        catch(Exception e){e.printStackTrace();} 
         
         
        File output=null;
        
        if(OrganismGroup==0)
            output=new File("GFPFungi.train.pred.arff");
        else if(OrganismGroup==1)
             output=new File("GFPMetazoa.train.pred.arff");
        
         try
            {
                FileWriter fw = new FileWriter(output);
                
                if(OrganismGroup == 0)
                  fw.write("@RELATION 'FungiGenomes-predictions'\n\n");
                else if(OrganismGroup == 1)
                   fw.write("@RELATION 'MetazoaGenomes-predictions'\n\n"); 
                fw.write("@ATTRIBUTE ID                                                               key\n");
                fw.write("@ATTRIBUTE class-a                                                          string\n");
                
                for(int i=0;i<oggomap.IndexToGO.keySet().size();i++){
                  fw.write("@ATTRIBUTE class-a-"+oggomap.IndexToGO.get(i)+"                                                          {1,0}\n");
                  
                }
                
                for(int i=0;i<oggomap.IndexToGO.keySet().size();i++){
                  fw.write("@ATTRIBUTE Original-p-"+oggomap.IndexToGO.get(i)+"                                                          numeric\n");
                  
                }
                
                fw.write("@ATTRIBUTE Original-models                                                          string\n");
                
                fw.write("\n");
                fw.write("@DATA\n");
                
                for(int i=0;i<cogOrder.size();i++){
                    fw.write("\""+cogOrder.get(i)+"\",\"\",");
                    ArrayList<Double> rl=realLabels.get(cogOrder.get(i));
                    for(int j=0;j<rl.size();j++)
                        fw.write(rl.get(j).intValue()+",");
                    
                    ArrayList<Double> pc=predictedScores.get(cogOrder.get(i));
                    
                    for(int j=0;j<pc.size();j++)
                        fw.write(pc.get(j)+",");
                    
                    fw.write("\"\"\n");
                            
                }
                
                fw.close();
            }
           catch(Exception e){
                e.printStackTrace();
            }  
    }
    
}
