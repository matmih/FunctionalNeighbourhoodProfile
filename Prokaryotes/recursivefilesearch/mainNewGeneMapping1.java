/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to create reduced mappings 
 */
public class mainNewGeneMapping1 {
    
    static public void main(String args[]){
         
        File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
        File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        String redCog="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\COGRed150.txt"; String redGO="GOAll.txt";
        
        GeneOGMapping geneOGMap=new GeneOGMapping();
        
        geneOGMap.loadGOGMapping(geneOgs);
        
        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(OgGOs);
        
        ReducedOgsGOs red=new ReducedOgsGOs();
        red.LoadReducedOGGos(new File(redCog), new File(redGO));
        
        HashMap<String,Integer> GOCount=new HashMap<>();
        
        for(String s:red.redCogs){
            System.out.println(s+" \n");
            if(!cgmap.CogGOmap.containsKey(s))
                continue;
            for(String go:cgmap.CogGOmap.get(s))
                if(!GOCount.containsKey(go))
                    GOCount.put(go, 1);
                else{
                    int count=GOCount.get(go);
                    count=count+1;
                    GOCount.put(go, count);
                }
             
       }
      
        FileWriter fw;
        String outputCounts="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\GOReducedDatasetTmp.txt";
        try{
        fw = new FileWriter(outputCounts); 
        
            
            Iterator<String> it1=GOCount.keySet().iterator();
            
            while(it1.hasNext()){
                String c=it1.next();
                int count=GOCount.get(c);
                fw.write(c+"\t"+count+"\n");
            }
            fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    
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

        red.LoadReducedOGGos(new File(redCog), new File(redGO));
        red.redGOs=rgt.goS;
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);   
        
             outputCounts="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\COGReducedDatasetNew.txt";
        try{
        fw = new FileWriter(outputCounts); 
        
            
            Iterator<String> it1=cgmap.CogGOmap.keySet().iterator();
            
            while(it1.hasNext()){
                String c=it1.next();
                fw.write(c+"\n");
            }
            fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
           outputCounts="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\GOReducedDatasetNew.txt";
        try{
        fw = new FileWriter(outputCounts); 
        
            
            Iterator<String> it1=cgmap.CogGOmap.keySet().iterator();
               HashMap<String,Integer> goCount=new HashMap<>();
            while(it1.hasNext()){
                String c=it1.next();
                
                ArrayList<String> gos=cgmap.CogGOmap.get(c);
                HashSet<String> gs=new HashSet(gos);
                
                Iterator<String> sIt=gs.iterator();
                
                while(sIt.hasNext()){
                    String gt=sIt.next();
                    if(!goCount.containsKey(gt)){
                        goCount.put(gt, 1);
                    }
                    else{
                        int i=goCount.get(gt);
                        ++i;
                        goCount.put(gt,i);
                    }                 
                }               
            }
            
            it1=goCount.keySet().iterator();
            while(it1.hasNext()){
                String c=it1.next();
                
                 fw.write(c+"\t"+goCount.get(c)+"\n");
                
            }
            fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    
}
