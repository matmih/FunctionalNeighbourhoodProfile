/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import org.javatuples.Pair;
/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute single class FNP dataset obtained from functional neighbourhoods on original and randomized data
 */
public class SingleClassPerformanceTestNumGenes {
     
    public static void main(String[] args) {
        
         final HashMap<String,HashSet<Pair<String,String>>> geneOGMap=new HashMap<String,HashSet<Pair<String,String>>>();
              
        String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa

         BufferedReader reader=null;
            
             try {
                     Path path =Paths.get(new File(geneOGMapFile).getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOGMap.containsKey(gene)){
                                geneOGMap.put(gene, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> ogs=geneOGMap.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og[]=tmp[i].split(":");
                               Pair<String, String> np = new Pair(og[0].trim(),og[1].trim());
                                
                               ogs.add(np);
                            }
                            geneOGMap.put(gene, ogs);
                    }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}

              final ArrayList<HashMap<String,HashSet<Pair<String,String>>>> ogMaps=new ArrayList<>();
              
               File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
        OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);
              
    OntologyTools.GeneOntology myGO=null;

          ArrayList<HashSet<String>> categories=new ArrayList<>();
        GODevider gd=new GODevider();
        try{
            categories=gd.ExtractCategories(new File("go_daily-termdb.obo-xml.gz"),cgmap);
        }
        catch(Exception e){
            e.printStackTrace();
        }  
        
        ArrayList<String> functions=new ArrayList<>();
        

        functions.add("GO0006022");
        functions.add("GO0006030");
        functions.add("GO0046189");
        functions.add("GO0042438");
        functions.add("GO0044550");
        functions.add("GO1901617");
                           
         String path="";
         
         File input=new File(path+"Metazoa_3_10_30k=10All.arff");
        
        Random rand = new Random();
        int runInd=0;
        ArrayList<HashMap<Integer,StringBuffer>> results=new ArrayList<>();
        
        for(int i=0;i<functions.size();i++){
            HashMap<Integer,StringBuffer> tmp=new HashMap<>();
            results.add(tmp);
        }            


         int seed=rand.nextInt((100 - 1) + 1) + 1;
          System.out.println("Iterations... "+runInd);
      for(int iter=0;iter<functions.size();iter++){  
         
           ARFF f=new ARFF();
          String GO=functions.get(iter);
                   
         f.loadARFF(input);
         f.assignLabels(GO);
         ArrayList<Double> maxs=new ArrayList<>();
         ArrayList<Double> mins=new ArrayList<>();
         f.attrStats(maxs, mins);
    
         f.writeARFFPO(path+"SingleFunction"+GO+"BPP.arff",categories,true);
         System.out.println("Stats...");
         System.out.println(functions.get(iter));
         System.out.println(f.attributes.size());
         f.writeARFFPOSF(path+"SingleFunction"+GO+"S.arff",functions.get(iter), categories, true, true, maxs, mins);
       }
      }
     }      

