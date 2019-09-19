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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute information about difference of LOR between original and randomized dataset
 */
public class FindLOGORDiff {
    static public void main(String arg[]){

         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File OgGOs30=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_30_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File semFile=new File("PairwiseGOSim.txt");
         
         String redCog="COGRed150.txt"; String redGO="GORed5.txt";

        GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(OgGOs);

        
        COGGOMap cgmap30=new COGGOMap();
        cgmap30.createCOGGOMapping(OgGOs30);
        
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
        cgmap30.reduceMap(red);
        geneOGMap.reduceMap(red);

        
        HashMap<String,String> goDescription=new HashMap<>();
       Iterator<String> itDesc=cgmap.GOtoIndex.keySet().iterator();
       
       while(itDesc.hasNext()){
           String go=itDesc.next();
           int goNum=Integer.parseInt(go.replaceAll("GO", ""));
           String desc=myGO.get(goNum).getName();
           
           goDescription.put(go, desc);
       }
       
         semanticsimilarity.GOMap mapFull=new semanticsimilarity.GOMap();
         mapFull.CreateGOMap(cgmap);
         System.out.println("Num GOs: "+cgmap.GOtoIndex.keySet().size());
         System.out.println("GO on index 0 "+cgmap.IndexToGO.get(0));
         GOSemanticSimilarity goSem=new GOSemanticSimilarity(mapFull.GOmap.keySet().size());
         goSem.readSimilarity(semFile);
        
          ArrayList<HashSet<String>> categories=new ArrayList<>();
        GODevider gd=new GODevider();
        try{
            categories=gd.ExtractCategories(new File("go_201401-termdb.obo-xml.gz"),cgmap);
        }
        catch(Exception e){
            e.printStackTrace();
        }  

        semanticsimilarity.GOMap mapMod=new semanticsimilarity.GOMap();
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
     
          String path="";

        FisherExactTest myFisher = new FisherExactTest();   
            
        ArrayList<String> functions=new ArrayList<>();      
     
           
        HashMap<String,Pair<Double,Double>> average=new HashMap<>();
        HashMap<String,Pair<Double,Double>> baseline=new HashMap<>();
        HashMap<String,Pair<Double,Double>> locNeigh=new HashMap<>();
        File functionFile= new File("BaselineOGs_600_.oob");
        File locNeighFile= new File("DistanceLocD_600_.oob");
        File baselineFile = new File("genome5k=5logDist.out");
        HashMap<String,Double> realFreq=new HashMap<>();
        
        try (BufferedReader bufRdr1 = new BufferedReader(new FileReader("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt")))
        {
            String line;
            
            while ((line = bufRdr1.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");

                int go = Integer.parseInt(parts[0]);
                String gString=go+"";
                while(gString.length()<7){
                    gString="0"+gString;
                }
                gString="GO"+gString;
                double generality = Double.parseDouble(parts[1]);
                realFreq.put(gString, generality);
                System.out.println("GO: "+gString);
            }
        }
        catch(Exception e){}

        BufferedReader reader;
         try {
             
      Path path1 =Paths.get(functionFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              average.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
         
      try {
             
      Path path1 =Paths.get(locNeighFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              locNeigh.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }     

         try {
      Path path1 =Paths.get(baselineFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int summary=0, dataSection=0, test=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Training error")){
              summary=1;
              continue;
          }
          
           if(line.contains("P100R") && summary==1){
               System.out.println("Data section reached baseline");
                  dataSection=1;
                 summary=0;
                  continue;
              }
          
           if(line.contains("Testing error") || line.equals(""))
              dataSection=0;

          if(dataSection==1){
              System.out.println("Entered data section");
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first baseline: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              baseline.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)));
          }
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            } 
        
         System.out.println("Average: "+average.keySet().size());
         System.out.println("DistLoc: "+locNeigh.keySet().size());
         System.out.println("Baseline: "+baseline.keySet().size());
         DecimalFormat df = new DecimalFormat("#.###");
         HashMap<String,Integer> GOCount=new HashMap<>();
        
         HashSet<String> functionsCond = new HashSet<>();
         
          
         Iterator<String> itFunc=average.keySet().iterator();
         
         while(itFunc.hasNext()){
             String go=itFunc.next();

                     functions.add(go);
             
         }
         
         System.out.println("Broj funkcija: "+functions.size());
         
        DataGOFunctionInfo datInfo=new DataGOFunctionInfo();
        datInfo.computeOccurence(cgmap);
        DataGOFunctionInfo datInfo1=new DataGOFunctionInfo();
        datInfo1.computeOccurence(cgmap30);
        
             double maxJS=0.0; String maxGO1="",maxGO2="";
        HashMap<String,ArrayList<String>> sigPairs=new HashMap<>();
          
        int criteria=1;//0-strict (logOR(GOx,GOy)>logOR(GOx,GOx), sigg.), 1-medium logOR(GOx,GOy)>0 sigg.
        int countPairs=0;
           
        try{
        
            FileWriter fw = new FileWriter("LORDifferenceProkaryots.txt");
            
            HashSet<String> bpOntol = new HashSet<>();
            
        for(int iter=0;iter<functions.size();iter++){ 
            ArrayList<String> tmp=new ArrayList<>();
            sigPairs.put(functions.get(iter), tmp);
            
            if(categories.get(0).contains(functions.get(iter)))
                bpOntol.add(functions.get(iter));
            
        int aInd=mapFull.GOmap.get(functions.get(iter));
        int aInd1=cgmap.GOtoIndex.get(functions.get(iter));
        int aInd2=cgmap30.GOtoIndex.get(functions.get(iter));
      
         LoadOR lor=new LoadOR();
         LoadOR lor1=new LoadOR();
         
     File ORInp=null;
     File ORInp1=null;
     
        ORInp = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
        ORInp1 = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesRandomNew\\"+"RandomContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
     lor.loadContValue(ORInp);
     lor1.loadContValue(ORInp1);
    
    for(int goindex=0;goindex<mapFull.GOmap.keySet().size();goindex++){
        String GOname=mapFull.GOmapNumeric.get(goindex);
        int iGO=mapFull.GOmap.get(GOname);
        int iGO1=cgmap.GOtoIndex.get(GOname);
        int iGO2=cgmap30.GOtoIndex.get(GOname);    
        
        ArrayList<Double> ct=lor.cont.get(GOname);
         ArrayList<Double> ct1=lor1.cont.get(GOname);

        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
        Double logOr=Math.log10(or)/Math.log10(2);
        
        double orR=((ct1.get(0)+0.5)/(ct1.get(1)+0.5))/((ct1.get(2)+0.5)/(ct1.get(3)+0.5));
        double logOrR=Math.log10(orR)/Math.log10(2);
       
        double diff;
        
      if(logOr>0){
         diff = logOr-logOrR;
         fw.write(diff+"\n");
      }
       
    }

     }

        fw.close();
        }
        catch(IOException e){
            e.printStackTrace();         
        }
    }
}
