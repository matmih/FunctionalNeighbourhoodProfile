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
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import org.apache.commons.io.FileUtils;
import static recursivefilesearch.COG.ENCODING;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create a file containing gene arrengement information
 */
 
public class CountGeneArrengementsCOG {
    
   public static void main(String [] args){
       
       String taxIDFilePath  = "/home/mmihelcic/GOdatasetGenerator/GOanotacije/TaxIDs.txt";
       String StrainsPath = "AllE.ColiStrains.txt";
       File root=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/tmpFiles/proba");
         File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
       File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        String redCog="COGRed150.txt"; String redGO="GORed5.txt";
        String[] extensions = {"ptt"};
        Boolean recursive = true;
        int k=Integer.parseInt(args[4].trim()), taxIdT=Integer.parseInt(args[0].trim());
        
        
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
        
        String COGTn=args[1].trim(), GOx=args[2].trim(), GOy=args[3].trim();
        
        
        COGOccurence cogOc = new COGOccurence();
        
        HashSet<String> cogsOfInterest = new HashSet<>();
        
        Iterator<String> itcogs = cgmap.CogGOmap.keySet().iterator();
        
        while(itcogs.hasNext()){
           String ct = itcogs.next();
           if(cgmap.CogGOmap.get(ct).contains(GOx)){
               cogsOfInterest.add(ct);
               System.out.println("Interesting: "+ct);
           }
        }
       
       for(String COGT:cogsOfInterest){ 
        CountGeneArrangements cga = new CountGeneArrangements();
       
        BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
         
         File StrainsFile = new File(StrainsPath);
         HashSet<Integer> strainTax = new HashSet<>();
         
         try {
             Path path =Paths.get(StrainsFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          strainTax.add(Integer.parseInt(line.trim()));
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
         
         
    for(int TTax:strainTax){
        if(!taxIDMap.containsValue(TTax))
            continue;
        cogOc.initialize(cgmap);
        if(!cgmap.CogGOmap.containsKey(COGT))
            System.out.println("Without function, choose another!");
        System.out.println("map size: "+cogOc.occurence.keySet().size());
        taxIdT = TTax;
         cogOc.computeCoocurence(root, taxIdT, extensions, recursive, taxIDFilePath , geneOGMap);
         System.out.println("COG occ: "+cogOc.occurence.get(COGT));
         if(cogOc.occurence.get(COGT) == 0)
             continue;
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGbaselinefeatures cbf=new COGbaselinefeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         
         int numGenomes=0, numFiles=0;
        
         //traverse main genome and plasmids for a given taxid, count COG frequency in genes of that organism
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                File input = (File) iterator.next();

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                if(taxid!=taxIdT)
                    continue;
                 
                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap); //change to 0 to get COGs from file
                  
                  
                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboorsVisual(cg,cgmap,geneOGMap,k);
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                
                cga.computeArrangements(taxid, sim, cgmap, geneOGMap, COGT, GOx, GOy,cogOc);
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }   
    }
    
    if(cga.genesCount.keySet().size()>0){
       try{ 
          FileWriter fw=new FileWriter("Arrangements"+COGT+".txt");
          Iterator<Integer> taxIt = cga.genesCount.keySet().iterator();
          while(taxIt.hasNext()){
             int tax = taxIt.next();
             fw.write("Organism taxId: "+tax+"\n");
             ArrayList<String> atemp = cga.genesCount.get(tax);
             
             for(int i1=0;i1<atemp.size();i1++){
                 fw.write(atemp.get(i1)+"\n");
             }
          }
          fw.close();
         }
         catch(Exception e){
             e.printStackTrace();
         }
       }
     }   
   }
}
