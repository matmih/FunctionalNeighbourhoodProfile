/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package semanticsimilarity;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import recursivefilesearch.COGGOMap;
import recursivefilesearch.GeneOGMapping;
import recursivefilesearch.ReducedGOTranslator;
import recursivefilesearch.ReducedOgsGOs;
/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to compute different measures of semantic similarity between all pairs of GO terms in prokaryotes
 */
public class Main {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
     
        ArrayList<ArrayList<Integer>> paths=new ArrayList<>();
         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");//pid2OG mapping
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");//OGGOmapping
         String redGO="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\GORed5.txt"; String redCog="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\COGRed150.txt";//reductions used in our work
         
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
              
        System.out.println("Number of functions in cgmap: "+cgmap.GOtoIndex.keySet().size());
        
        GOMap map=new GOMap();
        map.CreateGOMap(cgmap);
        map.printGOMapNumeric();
        File output=new File("C:\\Users\\matej\\Desktop\\Gosim.txt");
        map.saveGO(output);
        File fileWithUniprotFrequencyCounts=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");//Uniprot-freqs-from idmapping-2014-06-11
        map.loadFrequencies(fileWithUniprotFrequencyCounts);
        map.printFrequencies();
        File input1=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PP.arff");//File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\PP.arff");
        GOMap mapFull=new GOMap();
        mapFull.CreateGOMap(input1);
        mapFull.printGOMap();
        mapFull.loadFrequencies(fileWithUniprotFrequencyCounts);
        mapFull.printFrequencies();
        
        Graph g=new Graph(mapFull.GOmap.keySet().size());
        g.createAdjacency(mapFull);
        File outputM=new File("C:\\Users\\matej\\Desktop\\AdjMat.txt");
        g.saveAdjacency(outputM, mapFull.GOmap.keySet().size());
        
        mapFull.printGOMap();
        ArrayList<Path> Paths=new ArrayList<Path>();
        Path visited=new Path();
          
          Paths.clear();
         visited.path.clear();
      
        Wang w=new Wang();
     
        GOSimilarity sim=new GOSimilarity((map.GOmap.keySet().size()));//-1
       
        Object GO[]=map.GOmap.keySet().toArray();
                
        int missing=0;
        for(int i=0;i<GO.length;i++){
            System.out.println(GO[i]+" "+map.GOmap.get(GO[i]));
            if(mapFull.GOmap.containsKey(GO[i]) && !GO[i].equals("root")){
           sim.addSimilarity(((int)map.GOmap.get(GO[i]))+1, ((int)map.GOmap.get(GO[i]))+1, g.computeSimilarity(mapFull.GOmap.get(GO[i]), mapFull.GOmap.get(GO[i]), mapFull));
            }
            else{
                missing++;
              }
        }
               
        for(int i=0;i<GO.length-1;i++){
            for(int j=i+1;j<GO.length;j++){
                if(mapFull.GOmap.containsKey(GO[i]) && mapFull.GOmap.containsKey(GO[j]) && !GO[i].equals("root") && !GO[j].equals("root")){
        sim.addSimilarity(((int)map.GOmap.get(GO[i])+1), ((int)map.GOmap.get(GO[j])+1), g.computeSimilarity(mapFull.GOmap.get(GO[i]), mapFull.GOmap.get(GO[j]), mapFull));   
                }
            }
            if(i%100==0)
                 System.out.println("index: "+i);
        }
           
        for(int i=0;i<GO.length;i++)
            if(mapFull.GOmap.containsKey(GO[i]) && !GO[i].equals("root")){
                sim.addWangSimilarity(((int)map.GOmap.get(GO[i]))+1, ((int)map.GOmap.get(GO[i]))+1, w.sGO(mapFull.GOmap.get(GO[i]), mapFull.GOmap.get(GO[i]), 1.0, g));
            }
      
        for(int i=0;i<GO.length;i++){
            for(int j=i+1;j<GO.length;j++){
        if(mapFull.GOmap.containsKey(GO[i]) && mapFull.GOmap.containsKey(GO[j]) && !GO[i].equals("root") && !GO[j].equals("root"))
          sim.addWangSimilarity(((int)map.GOmap.get(GO[i]))+1, ((int)map.GOmap.get(GO[j]))+1, w.sGO(mapFull.GOmap.get(GO[i]), mapFull.GOmap.get(GO[j]), 1.0, g));  
            }
             if(i%100==0)
                 System.out.println("index: "+i);
        }
        
        File []outSim=new File[4];
         outSim[0]=new File("C:\\Users\\matej\\Desktop\\PairwiseGOSim.txt");
         outSim[1]=new File("C:\\Users\\matej\\Desktop\\PairwiseGOSimLin.txt");
         outSim[2]=new File("C:\\Users\\matej\\Desktop\\PairwiseGOSimRel.txt");
         outSim[3]=new File("C:\\Users\\matej\\Desktop\\PairwiseGOSimWang.txt");
        sim.writeSimilarity(outSim);
        
    }    

}
