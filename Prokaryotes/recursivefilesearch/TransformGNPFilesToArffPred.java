/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

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

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create comparative model performance file (contains comparative AUC/AUPRC scores)
 */
public class TransformGNPFilesToArffPred {
    
    public static void main(String [] args){
         
           //ClusPredictionInfo predInfo1=new ClusPredictionInfo();
        File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
        File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        
       // File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
       // File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        String redCog="COGRed150.txt"; String redGO="GORed5.txt";


        GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(OgGOs);
       
        ReducedGOTranslator rgt=new ReducedGOTranslator();
        rgt.ReadAndTranslate(new File(redGO));

    OntologyTools.GeneOntology myGO=null;

    try{
    myGO= new GeneOntology("go_201401-termdb.obo-xml");
    }
    catch(IOException e){
        e.printStackTrace();
    }

    HashSet<Integer> extendedGOs=new HashSet<>();

    for(int go:rgt.goI){
        Collection<GOTerm> parents = myGO.get(go).getAllParents();
    for ( GOTerm curPar : parents ) {
         if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
           extendedGOs.add(curPar.getId());
    }

    extendedGOs.add(go);

    }

    rgt.translate(extendedGOs);

        ReducedOgsGOs red=new ReducedOgsGOs();
        red.LoadReducedOGGos(new File(redCog), new File(redGO));
        red.redGOs=rgt.goS;
        System.out.println("Loading complete");
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);
        
        File cogPredScores=new File("GFPScoreMatrix.txt");
        File cogRealPreds=new File("AssociationTargetsFullAnot.txt");
        
        HashMap<String,ArrayList<Double>> realLabels=new HashMap<>();
        HashMap<String,ArrayList<Double>> predictedScores=new HashMap<>();
        ArrayList<String> cogOrder=new ArrayList<>();
        
         try (BufferedReader bufRdr = new BufferedReader(new FileReader("AssociationTargets.txt")))
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
        
        try (BufferedReader bufRdr = new BufferedReader(new FileReader("GFPScoreMatrix.txt")))
        {
            String line;
            int count=0;
            while ((line = bufRdr.readLine()) != null)
            {
                String tmp[]=line.split("\t");
                //System.out.println("tmp size: "+tmp.length);
                String go=cgmap.IndexToGO.get(count);
                
                for(int i=0;i<tmp.length;i++){
                    if(!predictedScores.containsKey(cogOrder.get(i)))
                             predictedScores.put(cogOrder.get(i), new ArrayList<Double>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), -1.0)));
                    double lab=Double.parseDouble(tmp[i]);
                    predictedScores.get(cogOrder.get(i)).set(count,lab);
                }
              count++;  
            }
            bufRdr.close();
        }
        catch(Exception e){e.printStackTrace();} 
         
         
        File output=new File("GFP.train.pred.arff");
         
         try
            {
                FileWriter fw = new FileWriter(output);
                
                fw.write("@RELATION 'BacterialGenomes-predictions'\n\n");
                fw.write("@ATTRIBUTE ID                                                               key\n");
                fw.write("@ATTRIBUTE class-a                                                          string\n");
                
                for(int i=0;i<cgmap.IndexToGO.keySet().size();i++){
                  fw.write("@ATTRIBUTE class-a-"+cgmap.IndexToGO.get(i)+"                                                          {1,0}\n");
                  
                }
                
                for(int i=0;i<cgmap.IndexToGO.keySet().size();i++){
                  fw.write("@ATTRIBUTE Original-p-"+cgmap.IndexToGO.get(i)+"                                                          numeric\n");
                  
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
                
                  //  fw.write(outputS);
                
                   
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            }  
    }
    
}
