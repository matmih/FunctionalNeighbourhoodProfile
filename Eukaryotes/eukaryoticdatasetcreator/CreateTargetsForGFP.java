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
import java.util.HashSet;
import java.util.Iterator;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create targets for GENMANIA software (GFP algorithm)
 */
 
public class CreateTargetsForGFP {
    static public void main(String args[]){
          int OrganismGroup = 0;//0 fungy, 1 metazoa
          
       // File ogFuncFile=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_1_30F.txt");
        File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        if(OrganismGroup == 1)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        //File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
        OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);
        System.out.println("Num functions: "+cgmap.GOtoIndex.keySet().size());
        System.out.println("Num ogs: "+cgmap.CogGOmap.keySet().size());
 
         File input=new File("FungyGFPBB.txt");
         ArrayList<String> cogsOrdered=new ArrayList<>();
         try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(input)))
        {
            String line;
            int count=0;
            while ((line = bufRdr1.readLine()) != null)
            {
                if(count==0){
                    count++;
                    continue;
                }
                line = line.trim();
                String tmp[]=line.split("\t");
                String tmp1[]=tmp[0].split(",");
                tmp1[0]=tmp1[0].replaceAll("\"", "");
                tmp1[0]=tmp1[0].trim();
                cogsOrdered.add(tmp1[0]);
            }
            bufRdr1.close();
        }
       catch(Exception e){
           e.printStackTrace();
       }
         
         
         File targets=new File("AssociationTargetsFungy.txt");
         
         try
            {
                FileWriter fw = new FileWriter(targets);
                
                for(int j=0;j<cogsOrdered.size();j++){
                    String COG=cogsOrdered.get(j);
                    System.out.println("COG: "+COG);
                    //System.out.println("Ex: "+COG);
                  //  System.out.println("Real: "+cgmap.CogGOmap.keySet().toArray()[0]);
                    ArrayList<String> GOs=cgmap.CogGOmap.get(COG);
                    fw.write(COG+"\t");
                    for(int i=0;i<cgmap.IndexToGO.keySet().size();i++){
                        if((i+1)<cgmap.IndexToGO.keySet().size()){
                            if(GOs.contains(cgmap.IndexToGO.get(i)))
                                fw.write(1+"\t");
                            else fw.write(0+"\t");
                        }
                        else{
                            if(GOs.contains(cgmap.IndexToGO.get(i)))
                                fw.write(1+"\n");
                            else fw.write(0+"\n");
                        }
                    }                       
                }
                   
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            }
    }         
}
