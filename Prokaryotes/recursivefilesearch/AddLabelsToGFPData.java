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
import java.util.HashMap;
import java.util.HashSet;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to add GO name to AUC-AUPRC scores computed by Gaussian Field Propagation
 */
 
public class AddLabelsToGFPData {
    
   static public void main(String [] args){
       
       File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
       String redCog="COGRed150.txt"; String redGO="GORed5.txt";
         COGGOMap cgmap=new COGGOMap();
          cgmap.createCOGGOMapping(OgGOs);
          
           OntologyTools.GeneOntology myGO=null;

    try{
    myGO= new GeneOntology("go_201401-termdb.obo-xml");
    }
    catch(IOException e){
        e.printStackTrace();
    }
          
          HashSet<Integer> extendedGOs=new HashSet<>();

            ReducedGOTranslator rgt=new ReducedGOTranslator();
        rgt.ReadAndTranslate(new File(redGO));
        
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
          
         
         File input=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\DifferentKResults\\GFPResultsk=5.txt");
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
                String goName=cgmap.IndexToGO.get(count);
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
         
         File output=new File("C:\\Users\\matej\\Downloads\\GeneMANIACode_v2\\GFPResultsLab.txt");
         
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
