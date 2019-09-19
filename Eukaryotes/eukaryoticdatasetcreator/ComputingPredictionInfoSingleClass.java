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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to count the number of genes that were assigned a function with a precision larger or equal to the precision threshold
  * input parameters
 * 0 - prediction input
 * 1 - GO function
 * 2 - precision threshold
 */
 
public class ComputingPredictionInfoSingleClass {

    static public void main(String[] args){
        
        final ClusPredictionInfo predInfo=new ClusPredictionInfo();
     
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

        File predInput=null;//new File("DistanceLocLogD.train.1.pred.arff");//BaselineOGsk=4.train.1.pred.arff

        predInput=new File(args[0].trim());
        
        String GO=args[1].trim();
        predInfo.readPredictionInfoSCXval(cgmap, GO ,predInput);
        
      Iterator<String> gos=cgmap.GOtoIndex.keySet().iterator();

       ArrayList<Integer> global=new ArrayList<>();
       global.add(0); global.add(0); global.add(0); global.add(0);
       //params: parallel (0-no, 1-yes), method (0-GFP, 1-LocLogDist), onlySpec (0-allFunctions, 1-IC>4 only)
        
        ArrayList<Double> labs=predInfo.originalLabels.get(GO);
        ArrayList<Double> classScores=predInfo.classifierScores.get(GO);
        
        System.out.println("Original labels: ");
        for(int i=0;i<labs.size();i++){
            System.out.print(labs.get(i)+" ");
        }
        System.out.println();
        
        System.out.println("Class labels: ");
        for(int i=0;i<classScores.size();i++){
            System.out.print(classScores.get(i)+" ");
        }
        System.out.println();
        
        System.out.println("Velicina: "+classScores.size());
        ArrayList<Double> sortedScores=predInfo.sortedScores(classScores);
         
         double maxRec=0.0, score=0.0, scorePrev=-1, precTreshold=0.5,precision=0.0;
         
         precTreshold=Double.parseDouble(args[2]);
        
          for(int i=0;i<sortedScores.size();i++){
               if(scorePrev==sortedScores.get(i)){
                   scorePrev=sortedScores.get(i);
                   continue;
               }
               else{
                ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, sortedScores.get(i));
                ArrayList<Double> preRec=predInfo.precRec(predictedClass, labs);
                System.out.println("PR: "+preRec.get(0));
                if(preRec.get(0)>=precTreshold){
                    if(preRec.get(1)> maxRec){
                        maxRec=preRec.get(1);
                        precision=preRec.get(0);
                        score=sortedScores.get(i);
                        scorePrev=sortedScores.get(i);
                    }
                }
          }
          }
          
          System.out.println("score: "+score);
          System.out.println("maxRec: "+maxRec);
          System.out.println("precision: "+precision); 
          
          if(precision<precTreshold)
              return;
          
         ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, score);
         
         for(int i=0;i<predictedClass.size();i++)
             System.out.print(predictedClass.get(i)+" ");
         
         System.out.println();
         
         ArrayList<Integer> difference=predInfo.numberOfNewGenes(predictedClass, labs, geneOGMap);

         System.out.println("difference: "+difference.get(0)+" "+difference.get(1)+" "+difference.get(2)+" "+difference.get(3));
         
         global.set(0, global.get(0)+difference.get(0));
         global.set(1, global.get(1)+difference.get(1));
         global.set(2, global.get(2)+difference.get(2));
         global.set(3, global.get(3)+difference.get(3));
         
      System.out.println("Global difference: ");
      System.out.println("Original: "+global.get(0)+" "+"Predicted: "+global.get(1)+" "+"Difference: "+global.get(2)+" New: "+global.get(3));
     }
    
}