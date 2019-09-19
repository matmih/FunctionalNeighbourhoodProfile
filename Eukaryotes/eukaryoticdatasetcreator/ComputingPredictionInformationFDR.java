/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantLock;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to compute the FDRs for predictions (eukaryotic organisms)
  * input parameters
 * 0 - use parallel computation: 0-no, 1-yes (1 should be used)
 * 1 - number of threads to be used for parallel computation
 * 2 - method to be used
 * 3 - use only specific GO functions
 * 4 -  precision threshold
 * 5 - Information content (IC) threhsold (lower bound)
 * 6 - Search type (0 - exhaustiv, 1 - heuristic)
 * 7 - Data type (0 - fungy, 1 - metazoa)
 * 8 - Only novel predictions (0 - false, 1 -true)
 */
public class ComputingPredictionInformationFDR {
    static public void main(String[] args){
        
        int numGOs = 0;
        
        final ClusPredictionInfo predInfo=new ClusPredictionInfo();
        final HashMap<String,HashSet<Pair<String,String>>> geneOGMap=new HashMap<String,HashSet<Pair<String,String>>>();
        HashMap<String,ArrayList<Double>> predictedOGGOFDRG = new HashMap<>(); 
        double precisionG = 0.0;
        
        int SearchType = Integer.parseInt(args[6].trim());//0-exhaustiv, 1-heuristic
        int dataType = Integer.parseInt(args[7].trim()); //0-fungy, 1-metazoa
       // String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa
        
       // String winMPE= "id_conversionReducedEns.tsv";//fungy
         String winMPE= "id_conversionReducedEnsMz.tsv";//metazoa
         
         if(dataType == 0)
             winMPE= "id_conversionReducedEns.tsv";
         else winMPE= "id_conversionReducedEnsMz.tsv"; 
         
       // String geneOGMapFile = "geneOGMappingFinalOKNewnonCH.txt";//fungy
            String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa
            
            if(dataType == 0)
                geneOGMapFile = "geneOGMappingFinalOKNewnonCH.txt";
            else geneOGMapFile = "geneOGMappingFinalOKNew.txt";
                
            
        String taxIDFilePath = "taxID.txt";
         Mappings map=new Mappings();
         if(dataType == 0)
          map.loadMappings(new File(winMPE)); 
         else
           map.loadMappingsMetazoa(new File(winMPE)); 
            
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
             
         HashMap<Integer,Integer> taxIDTranslationMap= new HashMap<>();
        
        //BufferedReader reader;
            
        String taxIDTranslation = "TaxIDTranslationFile.txt";
        
        File translationFile=new File(taxIDTranslation);
        
             try {
                     Path path =Paths.get(translationFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            int speciesTax = Integer.parseInt(tmp[1]);
                            int strainTax = Integer.parseInt(tmp[2]);
                            taxIDTranslationMap.put(strainTax, speciesTax);
                            }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}     
                   
        HashMap<Integer,HashMap<Integer, ArrayList<Integer>>> outputResults=new HashMap<>();     
        String targetSpecies = "targetSpecies.txt";
        File SpeciesFile = new File(targetSpecies);     
        HashSet<Integer> speciesTax = new HashSet<>();
         
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

       // geneOGMap.loadGOGMapping(geneOgs);
        final ArrayList<HashMap<String,HashSet<Pair<String,String>>>> ogMaps=new ArrayList<>();
        
        //geneOGMap.printGeneOGMapping();

        File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
        //File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        
        if(dataType==0)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");
        else ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
                   
        OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);

        for(int i=0;i<Integer.parseInt(args[1]);i++){
            HashMap<String,HashSet<Pair<String,String>>> tmp=new HashMap<String,HashSet<Pair<String,String>>>(geneOGMap);
            ogMaps.add(tmp); 
            System.out.println("Copied map "+(i+1));
        }
       
        File predInput=null;//new File("DistanceLocLogD.train.1.pred.arff");//BaselineOGsk=4.train.1.pred.arff
        //File predInput=new File("BaselineOGsk=4.train.1.pred.arff");
        
        int method=0;//Integer.parseInt(args[2]);
        
        String methods=args[2].trim();
         String tmpS[] = methods.split(",");
         
         ArrayList<Integer> methodsInt = new ArrayList<>();
         
         for(String s:tmpS)
             methodsInt.add(Integer.parseInt(s));
       
          CreateReducedMappingsFile rmf = new CreateReducedMappingsFile();
          //rmf.loadNogMembersMappings(new File("fuNOG.members.tsv")); //fungi
        // rmf.loadNogMembersMappingsMetazoa(new File("meNOG.members.tsv")); //metazoa
         
         if(dataType==0){
             rmf.loadNogMembersMappings(new File("fuNOG.members.tsv"));
         }
         else{
             rmf.loadNogMembersMappingsMetazoa(new File("meNOG.members.tsv"));
         }
          
        //File inputOrganism = new File(args[6]);//args[6] organism input file path
        String taxID = "";//args[6].trim();
            
        GenesAndOgsOrganism organismMappings = new GenesAndOgsOrganism();
        
           int test=1;
        
        System.out.println("Paramters: "+"Parallel: "+Integer.parseInt(args[0])+" Number of threads: "+Integer.parseInt(args[1]));
        
        //File root=new File("../FungiExtracted"); //fungi
         File root=new File("/home/mmihelcic/MatejQnap/MetazoaExtracted"); //metazoa
         
         if(dataType == 0)
             root=new File("../FungiExtracted");
         else root=new File("/home/mmihelcic/MatejQnap/MetazoaExtracted");
             
         String extensions[] = {"dat"};
        Boolean recursive=true;
        
        final ArrayList<GenesAndOgsOrganism> ogOrganismMaps = new ArrayList<>();
        final ArrayList<ClusPredictionInfo> predInfos = new ArrayList<>();
        //organismMappings.getGenesAndOgs(inputOrganism.getAbsolutePath(), geneOGMap);
         
        for(int i=0;i<Integer.parseInt(args[1]);i++){
            //HashMap<String,HashSet<Pair<String,String>>> tmp=new HashMap<String,HashSet<Pair<String,String>>>(geneOGMap);
            //ogMaps.add(tmp); 
            
            GenesAndOgsOrganism tmp1 = new GenesAndOgsOrganism(organismMappings);
            
           ogOrganismMaps.add(tmp1);
            
            System.out.println("Copied map "+(i+1));
        }
                
     for(int M:methodsInt){
        predInfos.clear();
          method=M;
         
         if(method==0){
          //  predInput=new File("BaselineOGsFungyNew.train.1.pred.arff");//BaselineOGsk=4.train.1.pred.arff //fungi
            // predInput=new File("BaselineOGsMetazoaNew1.train.1.pred.arff"); //metazoaTrain
             if(dataType==0)
                 predInput=new File("BaselineOGsFungyBB_oob.preds"); 
             else
                 predInput=new File("BaselineOGsMetazoaBB_oob.preds"); //metazoaOOB
        
         }
         else if(method==1){
            predInput=new File("DistanceLocLogD.train.1.pred.arff");
        }
         else if(method==2){
           // predInput=new File("genomeFungy.train.1.pred.arff"); //fungi
             //predInput=new File("genome5.train.1.pred.arff");//metazoa
             
             if(dataType == 0){
                 predInput=new File("genomeFungy.train.1.pred.arff");
             }
             else  predInput=new File("genome5.train.1.pred.arff");
             
         }
         else if(method==3){
             predInput=new File("GFP.train.pred.arff");
         }
        
        File inputData = null;
        if(method==0){
          // inputData = new File("BaselineOGsFungyBB_oob.preds");//fungi
            if(dataType == 0)
                inputData = new File("FungyFA_3_10_30k=10BB.arff");
            else inputData = new File("MetazoaFA_3_10_30k=10New.arff");
        }
      
        if(method==0){
            ARFF d = new ARFF();
            d.loadARFFTarget(inputData);
            
            File predInputTmp = null;
            
            if(dataType == 0){
                 predInputTmp=new File("genomeFungy.train.1.pred.arff");
             }
             else  predInputTmp=new File("genome5.train.1.pred.arff");
            
            predInfo.loadGOOrder(predInputTmp);//nesto se ovdje ne slaze
            predInfo.readPredictionInfoOOB(cgmap, d, predInput);
            numGOs = predInfo.gos.size();
            System.out.println("Num loaded gos: "+predInfo.gos.size());
        }
        else{
            //predInfo.loadGOOrder(predInput);
            predInfo.readPredictionInfo(cgmap, predInput);
            numGOs = predInfo.gos.size();
        }
        
        /*test=1;//remove after tests done
        if(test==1)
            continue;
        
        if(method==2)
            return;*/
        
        for(int i=0;i<Integer.parseInt(args[1]);i++){
            
            ClusPredictionInfo tmp1 = new ClusPredictionInfo(predInfo);
            
            predInfos.add(tmp1);
            
            System.out.println("Copied info "+(i+1));
        }
        
        GOMap mapMod=new GOMap();
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("GOFrequency.txt");//new File("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");//change
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
        precisionG = precTreshold;
         
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
               // final PriorityQueue qD=new PriorityQueue($trainingData.getNbRows()-1);
                final Lock lock = new ReentrantLock();
                final Lock lock1 = new ReentrantLock();
                final ArrayList<Integer> globResults=new ArrayList<>();
                final int numTr=Integer.parseInt(args[1]);
                final HashSet<String> allNewlyPredictedGenes = new HashSet<>();
                final HashMap<String,ArrayList<Double>> predictedOGGOFDR = predictedOGGOFDRG;
      
    try {
        int k=0;
       System.out.println("Num GOs: "+cgmap.GOtoIndex.size());
       
       while(gos.hasNext()){
          String go=gos.next();
     
          if(onlySpec==1){
              if(!mapMod.GOmap.containsKey(go)){
                  System.out.println(go+" not contained in mapMod.GoMap!");
                  System.out.println();
                  continue;
              }
              if(!mapMod.frequency.containsKey(mapMod.GOmap.get(go))){
                  System.out.println(go);
                  System.out.println("Index: "+mapMod.GOmap.get(go));
                  System.out.println("Not contained in mapMod.freq");
                  continue;
              }
          double ICTest=mapMod.frequency.get(mapMod.GOmap.get(go));
         // System.out.println("IC: "+ICTest); System.out.println("go: "+go);
             
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
          int tmp=k++;
          final GenesAndOgsOrganism ogsT=ogOrganismMaps.get((tmp)%numTr);
          int pnum=-1, countCur=0;
          
         // System.out.println("Semaphore MT: "+access.availablePermits());
         // access.acquire();

          
          final GenesAndOgsOrganism ogs=ogsT;
          final ClusPredictionInfo infpT=predInfos.get((tmp)%numTr);
          if(!infpT.originalLabels.containsKey(go))
              continue;
          final int pnumF = (tmp)%numTr;
          final Mappings mapF = map;
         final int SearchT = SearchType;
        //  if(k>24)
        //      break;//remove later
          
          System.out.println("k: "+tmp);
          
          countCur++;
          
          exec.submit(new Runnable() {
              int procNum = pnumF;
              GenesAndOgsOrganism ogs = ogsT;
              ClusPredictionInfo infP = infpT;
                        @Override
                        public void run() {
                          //  System.out.println("Computing in parallel threads...");
                            System.out.println("Processor number: "+procNum);
                            
                            HashMap<String,ArrayList<Double>> predictionsTmp = new HashMap<>();
                          

                            System.out.println("Function number: "+(count1));
                            System.out.println("go: "+goFin);
                            ArrayList<Double> labs=infP.originalLabels.get(goFin);
                            
                            ArrayList<Double> classScores=infP.classifierScores.get(goFin);
                              predictionsTmp.put(goFin, new ArrayList<Double>(Collections.nCopies(labs.size(), 1.0)));
                            System.out.println("Velicina class: "+classScores.size());
                            System.out.println("Velicina labs: "+labs.size());
                            ArrayList<Double> sortedScores=infP.sortedScores(classScores);
                            System.out.println("Sorted scores computed!");
                            System.out.println("Sorted scores size: "+sortedScores.size());//moguce rusenje
                                double maxRec=0.0, score=0.0, scorePrev=-1, precTreshold=0.5,precision=0.0;
        
                                 precTreshold=prec;
                                 
                                 
                                 ArrayList<Double> tmpAr = new ArrayList<>(Collections.nCopies(labs.size(), 0.0));
                                 ArrayList<Double> recAr = new ArrayList<>(Collections.nCopies(2, 0.0));
                            
                            long startTime = System.currentTimeMillis();
                            int located = 0, leftIndex=0, rightIndex=sortedScores.size()-1, successfulIndex=0, leftSucc=0;
                            //for(int i=0;i<sortedScores.size();i++){
                            int i=(leftIndex+rightIndex)/2;
                            System.out.println("SearchT: "+SearchT);
                            if(SearchT==1){
                            while(true){
                                        if(i<0 || i==sortedScores.size())
                                            System.out.println("i: "+i);
                                       //ArrayList<Double> predictedClass=infP.predictedClass1(classScores, tmpAr, precTreshold);//infP.predictedClass(classScores, sortedScores.get(i));                    
                                      infP.predictedClass1(classScores, tmpAr, sortedScores.get(i));//infP.predictedClass(classScores, sortedScores.get(i));
                                     ArrayList<Double> predictedClass = tmpAr;
                                     System.out.println("predicted class computed!");
                                     System.out.println("Predicted class: "+predictedClass.size());
                                    // ArrayList<Double> preRec=infP.precRec(predictedClass, labs);
                                    infP.precRec1(predictedClass, labs,recAr);
                                    ArrayList<Double> preRec = recAr;
                                     System.out.println(preRec.get(0)+" "+preRec.get(1));

                                         if(preRec.get(0)>=precTreshold){
                                               if(preRec.get(1)> maxRec){
                                                     maxRec=preRec.get(1);
                                                     precision=preRec.get(0);
                                                     score=sortedScores.get(i);
                                                     scorePrev=sortedScores.get(i);
                                              }
                                              // break;
                                              successfulIndex=i;
                                              leftSucc = leftIndex;
                                              rightIndex = i-1;
                                              i=(rightIndex+leftIndex)/2;
                                             }
                                         else{
                                             leftIndex=i+1;
                                             if(leftIndex<successfulIndex)
                                                 leftSucc = leftIndex;
                                             i=(rightIndex+leftIndex)/2;
                                         }                                        
                                  if(leftIndex>rightIndex)
                                      break;
                                 }
                            
                            System.out.println("successfulIndex: "+successfulIndex);
                            
                            for(i=leftSucc;i<=successfulIndex;i++){

                                       //ArrayList<Double> predictedClass=infP.predictedClass1(classScores, tmpAr, precTreshold);//infP.predictedClass(classScores, sortedScores.get(i));                    
                                      infP.predictedClass1(classScores, tmpAr, sortedScores.get(i));//infP.predictedClass(classScores, sortedScores.get(i));
                                     ArrayList<Double> predictedClass = tmpAr;
                                    // System.out.println("predicted class computed!");
                                    // System.out.println("Predicted class: "+predictedClass.size());
                                    // ArrayList<Double> preRec=infP.precRec(predictedClass, labs);
                                    infP.precRec1(predictedClass, labs,recAr);
                                    ArrayList<Double> preRec = recAr;
                                    // System.out.println(preRec.get(0)+" "+preRec.get(1));
                                        
                                         if(preRec.get(0)>=precTreshold){
                                               //if(preRec.get(1)> maxRec){
                                                     maxRec=preRec.get(1);
                                                     precision=preRec.get(0);
                                                     score=sortedScores.get(i);
                                                     scorePrev=sortedScores.get(i);
                                             // }
                                                    break;
                                             }
                                        }
                            }
                            else{
                                ArrayList<Double> tmp = predictionsTmp.get(goFin);
                                for(i=0;i<sortedScores.size();i++){
                                     infP.predictedClass1(classScores, tmpAr, sortedScores.get(i));//infP.predictedClass(classScores, sortedScores.get(i));
                                     ArrayList<Double> predictedClass = tmpAr;
                                    infP.precRec1(predictedClass, labs,recAr);
                                    ArrayList<Double> preRec = recAr;
                                                                  
                                 double prec = preRec.get(0);
                                 /*System.out.println("precision: "+prec);
                                 System.out.println("tmp size: "+tmp.size());
                                 System.out.println("predictedClass size: "+predictedClass.size());*/
                                 
                                /* if(predictedClass.size()>tmp.size())
                                     System.out.println("Wrong size: predicted"+predictedClass.size()+" tmp size: "+tmp.size());
                                 else System.out.println("size: "+predictedClass.size()+" "+tmp.size());*/
                                 
                                 for(int j=0;j<predictedClass.size();j++){
                                     if((tmp.get(j)>=(1.0-prec)) && (predictedClass.get(j)==1.0))
                                         tmp.set(j, (1.0-prec));
                                     //napuni lokalnu strukturu
                                     if(prec>=0.8 && predictedClass.get(j) == 1.0){
                                        System.out.println("P: "+infP.indexToCog.get(j)+" "+goFin+" "+tmp.get(j));
                                     }
                                 }
                                }
                                System.out.println("Local computation complete...");
                            }
                            
                            long endTime = System.currentTimeMillis();
                            System.out.println("Searching score cuttof: ");
                            System.out.println(((endTime-startTime)/1000));

                              try{
                                  lock.lock();
                                  try{
                                           ArrayList<Double> tmpDeep = new ArrayList<>();
                                           
                                           for(int in=0;in<predictionsTmp.get(goFin).size();in++)
                                               tmpDeep.add(predictionsTmp.get(goFin).get(in));
                                           
                                      //copy to global structure
                                      if(!predictedOGGOFDR.containsKey(goFin)){
                                            predictedOGGOFDR.put(goFin, tmpDeep);
                                            System.out.println("Function "+goFin+" set up");
                                            int count=0;
                                            for(int i1=0;i1<tmpDeep.size();i1++)
                                                    if(tmpDeep.get(i1)<0.2)
                                                        count++;
                                            System.out.println("Num OK: "+count+" "+goFin);
                                                }
                                      
                                  }finally{
                                            
                                           lock.unlock();
                                            //access.release();
                                           // System.out.println("Semaphore: "+access.availablePermits());
                                            
                                        }
                              }catch( Exception e ){
                                    e.printStackTrace();
                                    }
                        }
                          }); 
      }
       System.out.println("Shutting down");
       exec.shutdown();
       exec.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
    
     }
    catch(Exception e){
                        //System.out.println("Something went wrong...");
                        e.printStackTrace();
                    }
                finally {
                    /*System.out.println("Shutting down...");
                        exec.shutdown();*/

               /* try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                        System.out.println("Exception first shutdown");
                    }*/
    }
    }
  }
     
     System.out.println("Num functions: "+predictedOGGOFDRG.keySet().size());
     
      
      String info="";
      String outputName = "";
      
      int onlyNovel = 0;
      
      if(args.length>8)
          onlyNovel=Integer.parseInt(args[8]);
      
      if(args.length<9)
          outputName = "OutputPrecFDR.txt";
      else  outputName = "OutputPrecFDRNovel.txt";
          
      
        try{ 
          FileWriter fw=new FileWriter(outputName);
         
          fw.write("Method\tOGId\tGOId\tFDR"+"\n");
          fw.write("________________________________________________________________________________\n\n");
         
          System.out.println("gos size: "+predInfo.gos.size());
          
          System.out.println("struct: "+predictedOGGOFDRG.keySet().size());
         
          for(int met=0;met<predInfo.gos.size();met++){
              String go = predInfo.gos.get(met);
              System.out.println("pred go: "+go);
             /* go = go.replace("GO", "");
              int g = Integer.parseInt(go);
              go = "GO:"+g;*/
              System.out.println("pred go after: "+go);
              if(!predictedOGGOFDRG.containsKey(go))
                  continue;
          
              for(int j=0;j<predInfo.cogToIndex.keySet().size();j++){
                  if(predictedOGGOFDRG.get(predInfo.gos.get(met)).get(j)>(1-precisionG)){
                     // System.out.println("Final value: "+predictedOGGOFDRG.get(predInfo.gos.get(met)).get(j));
                      continue;
                  }
                  
               info="";
                if(method==0){
                  info+="GFN\t";
              }
              else if(method==1)
                      info+="LocDist\t";
              else if(method==2)
                  info+="Baseline\t";
              else if(method==3)
                  info+="GFP\t";
               info+=predInfo.indexToCog.get(j)+"\t";
               info+=predInfo.gos.get(met)+"\t";
               info+=predictedOGGOFDRG.get(predInfo.gos.get(met)).get(j);
              // predInfo.originalLabels
               
               ArrayList<Double> labs = predInfo.originalLabels.get(predInfo.gos.get(met));
               System.out.println("lab size: "+labs.size());
               System.out.println("only novel: "+onlyNovel);
               System.out.println("labs[j]: "+labs.get(j));
               int s=0;
               for(int ct = 0;ct<labs.size();ct++)
                   if(labs.get(ct) == 1.0)
                       s++;
               
               System.out.println("Num existing annotations: "+s);
               
               if(labs.get(j) == 1.0 && onlyNovel == 1)
                   continue;
               
               System.out.println("met: "+met);
               System.out.println("info: "+info);
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
