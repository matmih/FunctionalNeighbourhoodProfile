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
         File prokaryotInput=new File("prokaryotic.txt");
        HashSet<String> pkGOs=new HashSet<>();
        COGGOMap cgmap=new COGGOMap();
           
        
         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         String redGO="GORed5.txt"; String redCog="COGRed150.txt";//CogReduced50,GOReducedDataset5Occ
         
          GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

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
               
         File input=new File("AssociationLogDFullAnot.txt");
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
             
         File targets=new File("AssociationTargets.txt");
         
         try
            {
                FileWriter fw = new FileWriter(targets);
                
                for(int j=0;j<cogsOrdered.size();j++){
                    String COG=cogsOrdered.get(j);
                    System.out.println("COG: "+COG);
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
