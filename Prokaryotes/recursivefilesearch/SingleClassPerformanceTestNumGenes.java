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
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import semanticsimilarity.*;
/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create original features (BP ontology) and the corresponding dataset containing randomized features
 */
public class SingleClassPerformanceTestNumGenes {
     
    public static void main(String[] args) {
        
        File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
        File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        File semFile=new File("PairwiseGOSim.txt");
         
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
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);
        
        semanticsimilarity.GOMap mapFull=new semanticsimilarity.GOMap();
         mapFull.CreateGOMap(cgmap);
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
        
        ArrayList<String> functions=new ArrayList<>();
        functions.add("GO0045184");
        functions.add("GO0006418");
        functions.add("GO0006400");
        functions.add("GO0005975");
        functions.add("GO0008360");
        functions.add("GO0006281");
        functions.add("GO0006865");
        functions.add("GO0006260");
        functions.add("GO0006457");
        functions.add("GO0000271");
        functions.add("GO0009306");
        functions.add("GO0000902");
        functions.add("GO0006310");
        functions.add("GO0006508");
        functions.add("GO0019439");
        functions.add("GO0016051");
        functions.add("GO0006974");
        functions.add("GO0046365");
     
        HashMap<String,HashSet<String>> GOParents=new HashMap<>();

        for(int i=0;i<functions.size();i++){

        HashSet<String> gos=new HashSet<>();
        gos.add(functions.get(i));
        
        int go1=rgt.translateToIntSingle(functions.get(i));
        HashSet<Integer> pars=new HashSet<>();
        Collection<GOTerm> parents = myGO.get(go1).getAllParents();
    for ( GOTerm curPar : parents ) {
         if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
           pars.add(curPar.getId());
    }
    
    HashSet<String> parsString=rgt.translateToString(pars);
    gos.addAll(parsString);
        
        GOParents.put(functions.get(i), gos);
      }
        
         String path="";
         
         File input=new File(path+"BaselineOrganism1.arff");
        
        Classify cl=new Classify();
        FileDeleter del=new FileDeleter();
        Random rand = new Random();
        int runInd=0;
        ArrayList<HashMap<Integer,StringBuffer>> results=new ArrayList<>();
        
        for(int i=0;i<functions.size();i++){
            HashMap<Integer,StringBuffer> tmp=new HashMap<>();
            results.add(tmp);
        }            
        
        File FdrPvalsGZeroLOR=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pValsLORZeroFDR.txt");
         
          HashMap<Integer,ArrayList<Double>> fdrPvalsGZeroLOR=new HashMap<>();   
       
        for(int iter1=0;iter1<cgmap.GOtoIndex.keySet().size();iter1++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            fdrPvalsGZeroLOR.put(iter1, tmp);
        }
        
        BufferedReader bufRdr=null;
        
          try 
        {
            bufRdr = new BufferedReader(new FileReader(FdrPvalsGZeroLOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrpvalsGZeroLOR=Double.parseDouble(line.trim());
                fdrPvalsGZeroLOR.get(count).set(count1, fdrpvalsGZeroLOR);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
            bufRdr.close();
          }
           catch(Exception e){
               e.printStackTrace();
           }
        
        
          LoadOR lor=new LoadOR();
          
         int seed=rand.nextInt((100 - 1) + 1) + 1;
          System.out.println("Iterations... "+runInd);
      for(int iter=0;iter<functions.size();iter++){  
         
           ARFF f=new ARFF();
          String GO=functions.get(iter);
         OddsRatio or=new OddsRatio();
         
         File ORInp=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContAll/"+"ContingencyOut"+GO+".txt");
         lor.loadContValue(ORInp);
                            
         f.loadARFF(input);
         f.assignLabels(GO);
         ArrayList<Double> maxs=new ArrayList<>();
         ArrayList<Double> mins=new ArrayList<>();
         f.attrStats(maxs, mins);
    
         f.writeARFFPO(path+"SingleFunction"+GO+"BPP.arff",categories,true);
         
         f.writeARFFPOSF(path+"SingleFunction"+GO+"S.arff",functions.get(iter), categories, true, true, maxs, mins);
       }
      }
     }      

