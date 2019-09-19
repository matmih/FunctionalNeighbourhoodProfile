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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import static recursivefilesearch.COG.ENCODING;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the number of OGGOPairs obtaining a prediction with a FDR score below a threshold
  *params: parallel (0-no, 1-yes), method (0-GFP, 1-LocLogDist), onlySpec (0-allFunctions, 1-IC>4 only), precision, ICTreshold, novelAnnotations (1-report only novel predictions, not known before, 0 -report all)
 */
public class ComputingPredictionInformationFDR {
    static public void main(String[] args){
        
        final ClusPredictionInfo predInfo=new ClusPredictionInfo();
        int method = 0;
        File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
        File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        String redCog="COGRed150.txt"; String redGO="GORed5.txt";
        String taxIDFilePath  = "/home/mmihelcic/GOdatasetGenerator/GOanotacije/TaxIDs.txt";
        String targetSpecies = "targetSpecies.txt";
        double precisionG = 0.0;
        HashMap<String,ArrayList<Double>> predictedOGGOFDR = new HashMap<>();
        
           File SpeciesFile = new File(targetSpecies);
         HashSet<Integer> speciesTax = new HashSet<>();
         BufferedReader reader;
         
         try {
             Path path =Paths.get(SpeciesFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          speciesTax.add(Integer.parseInt(line.trim()));
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
         
         
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
        
        final GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);
        
        geneOGMap.reduceMap(red);
         
         String methods=args[2].trim();
         String tmpS[] = methods.split(",");
         
         ArrayList<Integer> methodsInt = new ArrayList<>();
         
         for(String s:tmpS)
             methodsInt.add(Integer.parseInt(s));
         
         File root=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/tmpFiles/proba");
        String[] extensions = {"ptt"};
        Boolean recursive = true;
        
        String taxID = "";
            
        GenesAndOgsOrganism organismMappings = new GenesAndOgsOrganism();
        
           int test=1;
        
        System.out.println("Paramters: "+"Parallel: "+Integer.parseInt(args[0])+" Number of threads: "+Integer.parseInt(args[1]));
    
        File predInput=null;
      
     for(int M:methodsInt){
        
         method=M;
        
         if(method==0){
            predInput=new File("BaselineOGsk=4_oob.preds");
        }
         else if(method==2){
            predInput=new File("genome5.train.1.pred.arff");
         }
         else if(method==3){
             predInput=new File("GFP.train.pred.arff");
         }
         
         File inputData = new File("BaselineOrganism1k=4.arff");
        
         if(method==0){
             
             File predInputTmp = new File("genome5.train.1.pred.arff");
                     
              ARFF d = new ARFF();
            d.loadARFFTarget(inputData);
           predInfo.loadGOOrder(predInputTmp);
            predInfo.readPredictionInfoOOB(cgmap, d, predInput);
         }
         else
        predInfo.readPredictionInfo(cgmap, predInput);
        
        semanticsimilarity.GOMap mapMod=new semanticsimilarity.GOMap();
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
        
      Iterator<String> gos=cgmap.GOtoIndex.keySet().iterator();
    
       ArrayList<Integer> global=new ArrayList<>();
       global.add(0); global.add(0); global.add(0); global.add(0);
        int count=0, parallel=Integer.parseInt(args[0]), onlySpec=Integer.parseInt(args[3]);
       //params: parallel (0-no, 1-yes), numThreads, method (0-GN, 1-LocLogDist, 2-baseline, 3-GFP), onlySpec (0-allFunctions, 1-IC>4 only), precTreshold, ICTreshold, inputOrganismFile
             
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
         HashSet<String> allNewGenes = new HashSet<>();
         ArrayList<Integer> difference=predInfo.numberOfNewGenesOrganism(predictedClass, labs, organismMappings,allNewGenes);

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
                final HashSet<String> allNewlyPredictedGenes = new HashSet<>();
                final int numTr=Integer.parseInt(args[1]);
                final ArrayList<ReentrantLock> locks = new ArrayList<>();
                final HashMap<String,ArrayList<Double>> predictedOGGOFDRL = predictedOGGOFDR;
      
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
          precisionG = prec;
          
          exec.submit(new Runnable() {
                        @Override
                        public void run() {
                            HashMap<String,ArrayList<Double>> predictionsTmp = new HashMap<>();
                            
                            System.out.println("Function number: "+(count1));
                            System.out.println("go: "+goFin);
                            ArrayList<Double> labs=predInfo.originalLabels.get(goFin);
                            predictionsTmp.put(goFin, new ArrayList<Double>(Collections.nCopies(labs.size(), 1.0)));
                            ArrayList<Double> classScores=predInfo.classifierScores.get(goFin);
                            System.out.println("Velicina class: "+classScores.size());
                            System.out.println("Velicina labs: "+labs.size());
                            ArrayList<Double> sortedScores=predInfo.sortedScores(classScores);
                            
                                double maxRec=0.0, score=0.0, scorePrev=-1, precTreshold=0.5,precision=0.0;
        
                                 precTreshold=prec;
                                
                             ArrayList<Double> tmp = predictionsTmp.get(goFin);    
                            for(int i=0;i<sortedScores.size();i++){
                                  if(scorePrev==sortedScores.get(i)){
                                       scorePrev=sortedScores.get(i);
                                       continue;
                                     }
                                   else{
                                     ArrayList<Double> predictedClass=predInfo.predictedClass(classScores, sortedScores.get(i));
                                     ArrayList<Double> preRec=predInfo.precRec(predictedClass, labs);
         
                                     double prec = preRec.get(0);
                                     
                                     for(int j=0;j<predictedClass.size();j++){
                                         if(tmp.get(j)>=(1.0-prec) && (predictedClass.get(j) == 1.0))
                                             tmp.set(j,1.0-prec);
                                     }

                                        }
                                 }
                             
                              try{
                                  lock.lock();
                                  try{
                                          if(!predictedOGGOFDRL.containsKey(goFin)){
                                              predictedOGGOFDRL.put(goFin, predictionsTmp.get(goFin));
                                              System.out.println("Function "+goFin+" set up");
                                          } 
                                  }finally{
                                            lock.unlock();
                                         }
                              }catch( Exception e ) {
                                    e.printStackTrace();
                                    }
                        }
                          }); 
      }
         
         
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
                
     }
     
     }
      
     System.out.println("NumFunctions: "+predictedOGGOFDR.keySet().size());
     
      String outputName = "OutputPrecFDR.txt";
      String info="";
      
      int novelAnnots = 0;
      
      if(args.length>5)
          novelAnnots = Integer.parseInt(args[6]);
      
        try{ 
          FileWriter fw=new FileWriter(outputName);
          
          fw.write("Method\tGOId\tGOId\tFDR"+"\n");
          fw.write("_________________________________________________________________________________________\n\n");
         
      for(int met=0;met<predInfo.gos.size();met++){
         
          String go = predInfo.gos.get(met);
          System.out.println("pred go: "+go);
          
          if(!predictedOGGOFDR.containsKey(go))
              continue;
          
          for(int j=0;j<predInfo.cogToIndex.keySet().size();j++){
              if(predictedOGGOFDR.get(predInfo.gos.get(met)).get(j)>(1-precisionG))
                  continue;
               info="";
             if(method==0){
                  info+="GFN\t";
              }
              else if(method==2)
                  info+="Baseline\t";
              else if(method==3)
                  info+="GFP\t";
          
         info+=predInfo.indexToCog.get(j)+"\t";
         info+=predInfo.gos.get(met)+"\t";
         info+=predictedOGGOFDR.get(predInfo.gos.get(met)).get(j);
         
         ArrayList<Double> labs = predInfo.originalLabels.get(predInfo.gos.get(met));
         
         if(novelAnnots == 1 && labs.get(j) == 1){
             continue;
         }
         
         fw.write(info+"\n");
          }
          
          fw.write("\n\n");
      }
          
          fw.close();
          
      }
         catch(Exception e){
             e.printStackTrace();
         }
    }
}
