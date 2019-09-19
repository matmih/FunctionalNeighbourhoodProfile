/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.*;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the number of genes obtaining a prediction at a given precision level 
 */

/**
 *
 * @author matej
 * input parameters
 * 0 - use parallel computation
 * 1 - number of parallel threads
 * 2 - method to be used (0 - GFP, 1 - LocLog)
 * 3 - use all functions (0) or only functions with IC>4 (1)
 * 4 - precisionTreshold (numeric [0,1])
 * 5- ICTreshold (numeric)
 */
public class ComputingPredictionInformation {
    
    static public void main(String[] args){
        
        final ClusPredictionInfo predInfo=new ClusPredictionInfo();
        File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
        File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        String redCog="COGRed150.txt"; String redGO="GORed5.txt";

           int test=1;
        
        System.out.println("Paramters: "+"Parallel: "+Integer.parseInt(args[0])+" Number of threads: "+Integer.parseInt(args[1]));
        
        final GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);
        final ArrayList<GeneOGMapping> ogMaps=new ArrayList<>();
        
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
        
        
        for(int i=0;i<Integer.parseInt(args[1]);i++){
            GeneOGMapping tmp=new GeneOGMapping(geneOGMap);
            ogMaps.add(tmp); 
            System.out.println("Copied map "+(i+1));
        }
       
        File predInput=null;
        
        int method=Integer.parseInt(args[2]);
        
         if(method==0){
            predInput=new File("BaselineOGsk=4.train.1.pred.arff");
        }
         else if(method==1){
            predInput=new File("DistanceLocLogD.train.1.pred.arff");
        }
         else if(method==2){
            predInput=new File("genome5.train.1.pred.arff");
         }
         else if(method==3){
             predInput=new File("GFP.train.pred.arff");
         }
        
        predInfo.readPredictionInfo(cgmap, predInput);
        
        semanticsimilarity.GOMap mapMod=new semanticsimilarity.GOMap();
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
        
      Iterator<String> gos=cgmap.GOtoIndex.keySet().iterator();       
       ArrayList<Integer> global=new ArrayList<>();
       global.add(0); global.add(0); global.add(0); global.add(0);
        int count=0, parallel=Integer.parseInt(args[0]), onlySpec=Integer.parseInt(args[3]);
       //params: parallel (0-no, 1-yes), method (0-GFP, 1-LocLogDist), onlySpec (0-allFunctions, 1-IC>4 only)
        
        
     if(parallel==0){  
          double ICThr=Double.parseDouble(args[5]);
          double positive=0.0,negative=0.0;
          
          if(ICThr>=0.0){
              positive=1.0;
              negative=0.0;
          }
          
          if(ICThr<0.0){
              positive=0.0;
              negative=1.0;
          }
          
      while(gos.hasNext()){
          String go=gos.next();
          
           if(onlySpec==1){
          double ICTest=mapMod.frequency.get(mapMod.GOmap.get(go));
             
        if(positive==1.0){  
        if(Math.abs(Math.log10(ICTest)/Math.log10(2))<=ICThr){
            ++count;
            continue;
        }
       }
        else if(negative==1.0){
             if(Math.abs(Math.log10(ICTest)/Math.log10(2))>=ICThr){
            ++count;
            continue;
        }
        }
      }
          
          System.out.println("Function number: "+(++count));
          go=go.replace("GO", "GO:");
          System.out.println("go: "+go);
        ArrayList<Double> labs=predInfo.originalLabels.get(go);
        ArrayList<Double> classScores=predInfo.classifierScores.get(go);
        System.out.println("Velicina: "+classScores.size());
        ArrayList<Double> sortedScores=predInfo.sortedScores(classScores);
          
         double maxRec=0.0, score=0.0, scorePrev=-1, precTreshold=0.5,precision=0.0;
         
         precTreshold=Double.parseDouble(args[4]);
        
          for(int i=0;i<sortedScores.size();i++){
               if(scorePrev==sortedScores.get(i)){
                   scorePrev=sortedScores.get(i);
                   continue;
               }
               else{
                ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, sortedScores.get(i));
                ArrayList<Double> preRec=predInfo.precRec(predictedClass, labs);
         
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
          
         ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, score);
         ArrayList<Integer> difference=predInfo.numberOfNewGenes(predictedClass, labs, geneOGMap);

         System.out.println("difference: "+difference.get(0)+" "+difference.get(1)+" "+difference.get(2)+" "+difference.get(3));
         
         global.set(0, global.get(0)+difference.get(0));
         global.set(1, global.get(1)+difference.get(1));
         global.set(2, global.get(2)+difference.get(2));
         global.set(3, global.get(3)+difference.get(3));
         
      }
      
      System.out.println("Global difference: ");
      System.out.println("Original: "+global.get(0)+" "+"Predicted: "+global.get(1)+" "+"Difference: "+global.get(2)+" New: "+global.get(3));
     }
     else{
        double ICThr=Double.parseDouble(args[5]);
         ExecutorService exec = Executors.newFixedThreadPool(Integer.parseInt(args[1]));
                final Lock lock = new ReentrantLock();
                final Lock lock1 = new ReentrantLock();
                final ArrayList<Integer> globResults=new ArrayList<>();
                final int numTr=Integer.parseInt(args[1]);
      
    try {
        int k=0;
       while(gos.hasNext()){
          String go=gos.next();
     
          if(onlySpec==1){
          double ICTest=mapMod.frequency.get(mapMod.GOmap.get(go));
             
        if(Math.abs(Math.log10(ICTest)/Math.log10(2))<=ICThr){
            ++count;
            continue;
        }
      }
          
          go=go.replace("GO", "GO:");
          
          final String goFin=go;
          final int count1=(++count);
          final double prec= Double.parseDouble(args[4]);
          final GeneOGMapping ogs=ogMaps.get((k++)%numTr);
          
          exec.submit(new Runnable() {
                        @Override
                        public void run() {
                            System.out.println("Function number: "+(count1));
                            System.out.println("go: "+goFin);
                            ArrayList<Double> labs=predInfo.originalLabels.get(goFin);
                                          
                            ArrayList<Double> classScores=predInfo.classifierScores.get(goFin);
                            System.out.println("Velicina class: "+classScores.size());
                            System.out.println("Velicina labs: "+labs.size());
                            ArrayList<Double> sortedScores=predInfo.sortedScores(classScores);
                            
                                double maxRec=0.0, score=0.0, scorePrev=-1, precTreshold=0.5,precision=0.0;
        
                                 precTreshold=prec;
                                
                            for(int i=0;i<sortedScores.size();i++){
                                  if(scorePrev==sortedScores.get(i)){
                                       scorePrev=sortedScores.get(i);
                                       continue;
                                     }
                                   else{
                                     ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, sortedScores.get(i));
                                     ArrayList<Double> preRec=predInfo.precRec(predictedClass, labs);
         
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
                             System.out.println("precTresh: "+precTreshold);
                        ArrayList<Double> predictedClass=new ArrayList<>(Collections.nCopies(labs.size(), 0.0));
                           if(precision>=precTreshold) 
                                 predictedClass=predInfo.predictedClass(classScores, score);
                           
                           ArrayList<Integer> difference=null;
                           
                                     difference=predInfo.numberOfNewGenes(predictedClass, labs, ogs/*geneOGMap*/);
                            
                            if(precision<precTreshold){
                                difference.set(1,0); difference.set(2, 0); difference.set(3, 0);
                            }
                             System.out.println("difference: "+difference.get(0)+" "+difference.get(1)+" "+difference.get(2)+" "+difference.get(3));
                             
                              try{
                                  lock.lock();
                                  try{
                                    globResults.add(difference.get(0));
                                    globResults.add(difference.get(1));
                                    globResults.add(difference.get(2));
                                    globResults.add(difference.get(3));
                                  }finally{
                                            lock.unlock();
                                        }
                              }catch( Exception e ) {
                                    e.printStackTrace();
                                    }
                        }
                          }); 
      }
      
      System.out.println("Global difference: ");
      System.out.println("Original: "+global.get(0)+" "+"Predicted: "+global.get(1)+" "+"Difference: "+global.get(2));  
         
         
     }
    catch(Exception e){
                        e.printStackTrace();
                    }
                finally {
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                
                for(int i=0;i<globResults.size();i++){
                    if(i%4==0){
                        global.set(0, global.get(0)+globResults.get(i));
                    }
                    else if(i%4==1){
                        global.set(1, global.get(1)+globResults.get(i));
                    }
                    else if(i%4==2){
                        global.set(2, global.get(2)+globResults.get(i));
                    }
                     else if(i%4==3){
                        global.set(3, global.get(3)+globResults.get(i));
                    }
                }
                
                System.out.println("Global difference: ");
                System.out.println("Original: "+global.get(0)+" "+"Predicted: "+global.get(1)+" "+"Difference: "+global.get(2)+" new:"+global.get(3));           
      }       
    }
    
    
}
