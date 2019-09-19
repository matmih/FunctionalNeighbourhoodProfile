/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;

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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import org.apache.commons.io.FileUtils;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import static recursivefilesearch.COG.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to create FNP features, features used to train the k-NN classifier and features used to train the GFP classifier
  *this class also contains functions to compute contingency tables, to count COG/NOG occurrences and GO occurrence and datasets to compute accretion
 */
public class MultiGenomeSet {
  File root;

  MultiGenomeSet(String rootPath){
      root=new File(rootPath);
  }
  
  public void countCogNogsGos(String taxIDFilePath, String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap, String output,String output1,int isTrain, boolean randomize){
       
       
       BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashMap<Integer,HashSet<String>> taxCogNogs=new HashMap<>();
         HashMap<Integer,HashSet<String>> taxGOs = new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         HashMap<Integer,Integer> NumPlasmids=new HashMap<>();
         
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
       
        try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGbaselinefeatures cbf=new COGbaselinefeatures();
         CircularDistance dist=new CircularDistance();
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
         Random rand=new Random();

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {

                 File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);
                  String fileName=input.getName().substring(0,input.getName().length()-4);
                  int tid=taxIDMap.get(fileName);
                  
                   if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }
                  
                  for(int z=0;z<cg.cogs.size();z++){
                     String gene= cg.cogs.get(z).getValue2();
                     if(!geneOGMap.geneOGsMap.containsKey(gene))
                         continue;
                    HashSet<String> cnogs =  geneOGMap.geneOGsMap.get(gene);
                    
                    Iterator<String> cit=cnogs.iterator();
                    
                    while(cit.hasNext()){
                        String cog=cit.next();
                        if(!taxCogNogs.containsKey(tid))
                            taxCogNogs.put(tid, new HashSet<String>());
                        
                        taxCogNogs.get(tid).add(cog);
                        
                        if(!taxGOs.containsKey(tid))
                            taxGOs.put(tid, new HashSet<String>());
                        
                        ArrayList<String> gosa=cgmap.CogGOmap.get(cog);
                       // System.out.println("cog: "+cog);
                        
                        if(!cgmap.CogGOmap.keySet().contains(cog)){
                                 continue;
                          }
                        
                        for(int z1=0;z1<gosa.size();z1++)
                            taxGOs.get(tid).add(gosa.get(z1));
                    }
                    
                  }
                  
                  if(randomize==true)
                      cg.randomizeLocation(rand,100);
                  
                  numGenomes++;
                 System.out.println("Broj obradenih dokumenata: "+numGenomes);
            }

            File outputCounts=new File(output);
            
            
            FileWriter fw = new FileWriter(outputCounts); 
        
            Iterator<Integer> it=taxCogNogs.keySet().iterator();
            HashMap<String,Integer> countsCogs=new HashMap<>();
            HashMap<String,Integer> countsGOs=new HashMap<>();
            
            while(it.hasNext()){
                int tidS=it.next();
                
                HashSet<String> tcn=taxCogNogs.get(tidS);
                
                Iterator<String> itn=tcn.iterator();
                
                while(itn.hasNext()){
                String c=itn.next();
                
                if(!countsCogs.containsKey(c))
                    countsCogs.put(c, 1);
                else{
                    int cc=countsCogs.get(c);
                    cc=cc+1;
                    countsCogs.put(c, cc);
                }
            }
                
              tcn=taxGOs.get(tidS);
              itn=tcn.iterator();
              
                while(itn.hasNext()){
                String c=itn.next();
                
                if(!countsGOs.containsKey(c))
                    countsGOs.put(c, 1);
                else{
                    int cc=countsGOs.get(c);
                    cc=cc+1;
                    countsGOs.put(c, cc);
                }
            }
          }
            
            Iterator<String> it1=countsCogs.keySet().iterator();
            
            while(it1.hasNext()){
                String c=it1.next();
                int count=countsCogs.get(c);
                fw.write(c+"\t"+count+"\n");
            }
            fw.close();
            
         outputCounts=new File(output1);
            
            
            fw = new FileWriter(outputCounts); 
            
            it1=countsGOs.keySet().iterator();
            
            while(it1.hasNext()){
                String c=it1.next();
                int count=countsGOs.get(c);
                fw.write(c+"\t"+count+"\n");
            }
            
        /*for(int i=0;i<cogNFC.CogIndex.keySet().size();i++){
            ArrayList<Double> tmp=finalStat.get(i);
            fw.write(cogNFC.IndexCog.get(i)+" \t ");
            for(int itt=0;itt<tmp.size();itt++){
                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+" \t ");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
            }
        }*/
        fw.close();
            
           //tt.createLocationTrainHeader(output, cgmap,sim, isTrain);
           //tt.appendLocationRowsToSet(output, cgmap, sim, isTrain);

   }catch (Exception e) {
            e.printStackTrace();
        }
 }
   
 
   void createAssociationGraph(String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,int useLogDist, String output,int isTrain){
       
        try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGbaselinefeatures cbf=new COGbaselinefeatures();
         CircularDistance dist=new CircularDistance();
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
         LocationSimilarity sim=new LocationSimilarity();
         sim.initializeSimilarity(cgmap);
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                 File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }
                  if(useLogDist==0)
                    sim.computeSimilarity(cg, cgmap, geneOGMap,dist); 
                  else if(useLogDist==1)
                      sim.computeLogSimilarity(cg, cgmap,geneOGMap ,dist);
                  else
                      sim.computeLogSimilarityLoc(cg, cgmap, dist);
                  numGenomes++;
                 System.out.println("Broj obradenih dokumenata: "+numGenomes);
            }
            sim.normalize1();
            sim.transformToWeights();
            
            if(useLogDist==1)
                output=(output.substring(0,output.length()-5)+"LogD.arff");
            else if(useLogDist>1)
                output=(output.substring(0,output.length()-5)+"LogL.arff");
                             
          // tt.createLocationTrainHeader(output, cgmap,sim, isTrain);//crete graph dataset
            tt.createAssociationGraphHeader(output, cgmap, sim, isTrain);
            tt.appendWeightsToGraph(output, cgmap, sim, isTrain);
          // tt.appendLocationRowsToSet(output, cgmap, sim, isTrain);
           sim.COGIndex.clear();
           sim.locationN.clear();
   }catch (Exception e) {
            e.printStackTrace();
        }
 }
   
    ArrayList<Double> numGenes(String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,int useLogDist, String output,int isTrain){
       
        ArrayList<Double> res = new ArrayList<>();
        res.add(0.0); res.add(0.0);
        
        try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         int numGenomes=0;

            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                 File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }
                  
                  res.set(0, res.get(0)+cg.cogs.size());
                  
                  HashSet<String> gos = new HashSet<>();
                  for(Triplet<Integer,Integer,String> t:cg.cogs){
                      if(geneOGMap.geneOGsMap.containsKey(t.getValue2())){
                         HashSet<String> ogs = geneOGMap.geneOGsMap.get(t.getValue2());
                          for(String s:ogs){
                              if(!cgmap.CogGOmap.containsKey(s))
                                  continue;
                              gos.addAll(cgmap.CogGOmap.get(s));
                          }
                           res.set(1, res.get(1)+gos.size());
                           gos.clear();
                      }
                  }
                 
            }     
   }catch (Exception e) {
            e.printStackTrace();
        }
        return res;
 }
   
   void createCoocurence(String[] extensions, boolean recursive, int k,COGGOMap cgmap,GeneOGMapping geneOGMap ,String output,int isTrain){
    
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGCoocurence cbf=new COGCoocurence();
         cbf.initializeSimilarity(cgmap);
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                cbf.createFeatures(sim, cgmap);
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            tt.createCoocurenceTrainHeader(output, cgmap,cbf, isTrain);
            tt.appendCoocurenceRowsToSet(output, cgmap, cbf, isTrain);
         //tt.createSparseTrainHeader(output, cgmap, isTrain, numGenomes);
        // tt.appendSparseRowsToSet(output, cgmap, csf, isTrain, numGenomes);
         //csf.COGsparsefeaturesmap.clear();
            cbf.coocurence.clear();
            cbf.COGIndex.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
   
   void createDistanceBaseline(String[] extensions, boolean recursive,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output,int isTrain, boolean randomize){
    
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGdistancefeatures cdf=new COGdistancefeatures();
         cdf.initializeDistances(cgmap);
         //COGbaselinefeatures cbf=new COGbaselinefeatures();
         //COGSparsefeatures csf=new COGSparsefeatures();
         Random rand=new Random();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);
                  
                  if(randomize==true)
                      cg.randomize(rand);
                
                cdf.ComputeDistances(cg, cgmap, geneOGMap);
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
            }
            tt.createDistanceTrainHeader(output,cdf ,isTrain);
            tt.appendDistanceRowsToSet(output, cgmap, cdf, isTrain);
           // tt.createTrainHeader(output, cgmap, isTrain);
            //tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
            cdf.distances.clear();
            cdf.maxCoord.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
  
     void computeAccritionDiffOnt(String taxIDFilePath,String accretionFilePath,String predictionFilePath,String[] extensions, boolean recursive, ArrayList<HashSet<String>> categories  ,int k,COGGOMap cgmap,GeneOGMapping geneOGMap, String output,int isTrain, boolean randomize){
    
        // String fileName=input.getName().substring(0,input.getName().length()-4);
        //        int taxid=taxIDMap.get(fileName);
         
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
          
         
         HashMap<String,Double> accretions = new HashMap<>();
          HashMap<String,HashSet<String>> predictions = new HashMap<>();
         
            try {
                File accretionFile = new File(accretionFilePath);
             Path path =Paths.get(accretionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(",");
          String og = tok[0].trim();
          //fix og
          int num = Integer.parseInt(og);
            
         int  nDigGo=0;
         
               while(num>0){
                  nDigGo++;
                  num/=10;
             }
         
         String GOName = "GO";
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=og;
          
          Double accr = Double.parseDouble(tok[1].trim());
          
          accretions.put(GOName, accr);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
            
            
       try {
                File predictionFile = new File(predictionFilePath);
             Path path =Paths.get(predictionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null){
          String[] tok=line.split("\t");
          String og = tok[0].trim();

          for(int j=1;j<tok.length;j++){
              if(!predictions.containsKey(og))
                  predictions.put(og, new HashSet<String>());
              
              predictions.get(og).add(tok[j].trim());
          }
      }
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}     
       
       System.out.println("predictions size: "+predictions.keySet().size());
       System.out.println("accretion size: "+accretions.keySet().size());
         
       double onlyKnownAnnotationAccretionBP=0.0, onlyKnownAnnotationAccretionMF=0.0, onlyKnownAnnotationAccretionCC=0.0, 
               newlyPredictedAnnotationAccretionBP=0.0, newlyPredictedAnnotationAccretionMF=0.0, newlyPredictedAnnotationAccretionCC=0.0, 
               knownAndPredictedAccretionBP = 0.0,  knownAndPredictedAccretionMF = 0.0,  knownAndPredictedAccretionCC = 0.0;
       double numGenes = 0.0; 
       
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         Random rand = new Random();
         
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap); //change to 0 to get COGs from file
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                   int taxid=taxIDMap.get(fileName);
                   
                 //  if(taxid!=511145)
                    //   continue;
                  
                  //create accretionCounting
                  
                  for(int g=0;g<cg.ancogs.size();g++){
                      String gen = cg.ancogs.get(g).getValue2();
                      HashSet<String> ogs = geneOGMap.geneOGsMap.get(gen);
                      
                       ArrayList<HashSet<String>> allgos = new ArrayList<>();
                      ArrayList<HashSet<String>> predictedgos = new ArrayList<>();
                      
                      for(int z=0;z<categories.size();z++){
                          allgos.add(new HashSet<String>());
                          predictedgos.add(new HashSet<String>());
                      }
                      
                      int useful = 0;
                      
                      for(String og:ogs){
                         ArrayList<String> go = cgmap.CogGOmap.get(og);
                         if(!cgmap.CogGOmap.containsKey(og))
                             continue;
                         
                         if(go.size()>0)
                             useful = 1;
                         
                         for(String gg:go){
                             if(categories.get(0).contains(gg))
                                     allgos.get(0).add(gg);
                             else if(categories.get(1).contains(gg))
                                     allgos.get(1).add(gg);
                             else if(categories.get(2).contains(gg))
                                     allgos.get(2).add(gg);
                             else 
                                 System.out.println("Can not find a function: "+gg);
                      }
                         
                         if(predictions.containsKey(og)){
                             HashSet<String> f = predictions.get(og);
                             for(String gg:f){
                              if(categories.get(0).contains(gg))
                                     predictedgos.get(0).add(gg);
                             else if(categories.get(1).contains(gg))
                                     predictedgos.get(1).add(gg);
                             else if(categories.get(2).contains(gg))
                                     predictedgos.get(2).add(gg);
                             }
                             
                         }
                         
                    }
  
                      if(useful == 1)
                            numGenes++;
                   
                   for(int z=0;z<3;z++){     
                      for(String go:allgos.get(z)){
                          if(!accretions.containsKey(go))
                              System.out.println("Error: "+go);
                          if(!predictedgos.get(z).contains(go)){
                               if(z==0)
                                    onlyKnownAnnotationAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   onlyKnownAnnotationAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   onlyKnownAnnotationAccretionCC+=accretions.get(go);
                          }
                          else if(predictedgos.get(z).contains(go)){
                                if(z==0)
                                    knownAndPredictedAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   knownAndPredictedAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   knownAndPredictedAccretionCC+=accretions.get(go);
                          }
                      }
                      
                      for(String go:predictedgos.get(z)){
                               if(!allgos.get(z).contains(go)){
                               if(z==0)
                              newlyPredictedAnnotationAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   newlyPredictedAnnotationAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   newlyPredictedAnnotationAccretionCC+=accretions.get(go);
                          }
                      }
                  }
                      
                  }
                 // numGenes+=cg.ancogs.size();

            }
            
            String outputTmp = output;
            
            for(int i=0;i<3;i++){  
              String o[] = outputTmp.split("\\.");
              String t = "";
              
              if(i==0)
                  t="BP";
              else if(i==1)
                  t="MF";
              else if(i==2)
                  t="CC";
              
              output = o[0]+t+"."+o[1];
            FileWriter fw = new FileWriter(output);
            
            double knownAA = 0.0, knownPredAA = 0.0, newlyPredAA = 0.0;
            
            if(i==0){
                  knownAA = onlyKnownAnnotationAccretionBP;
                  knownPredAA = knownAndPredictedAccretionBP;
                  newlyPredAA = newlyPredictedAnnotationAccretionBP;
            }
              else if(i==1){
                   knownAA = onlyKnownAnnotationAccretionMF;
                   knownPredAA = knownAndPredictedAccretionMF;
                   newlyPredAA = newlyPredictedAnnotationAccretionMF;
              }
              else if(i==2){
                   knownAA = onlyKnownAnnotationAccretionCC;
                   knownPredAA = knownAndPredictedAccretionCC;
                   newlyPredAA = newlyPredictedAnnotationAccretionCC;
              }
            
            fw.write("Only known: "+(knownAA/numGenes)+"\n");
            fw.write("Known and predicted: "+(knownPredAA/numGenes)+"\n");
            fw.write("Only predicted: "+(newlyPredAA/numGenes)+"\n");
            fw.close();
            
            System.out.println("Only known: "+(knownAA/numGenes));
            System.out.println("Known and predicted: "+(knownPredAA/numGenes));
            System.out.println("Only predicted: "+(newlyPredAA/numGenes));
          }
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
      
     void computeAccrition(String taxIDFilePath,String accretionFilePath,String predictionFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap,GeneOGMapping geneOGMap, String output,int isTrain, boolean randomize){
    
        // String fileName=input.getName().substring(0,input.getName().length()-4);
        //        int taxid=taxIDMap.get(fileName);
         
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
          
         
         HashMap<String,Double> accretions = new HashMap<>();
          HashMap<String,HashSet<String>> predictions = new HashMap<>();
         
            try {
                File accretionFile = new File(accretionFilePath);
             Path path =Paths.get(accretionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(",");
          String og = tok[0].trim();
          //fix og
          int num = Integer.parseInt(og);
            
         int  nDigGo=0;
         
               while(num>0){
                  nDigGo++;
                  num/=10;
             }
         
         String GOName = "GO";
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=og;
          
          Double accr = Double.parseDouble(tok[1].trim());
          
          accretions.put(GOName, accr);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
            
            
       try {
                File predictionFile = new File(predictionFilePath);
             Path path =Paths.get(predictionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null){
          String[] tok=line.split("\t");
          String og = tok[0].trim();

          for(int j=1;j<tok.length;j++){
              if(!predictions.containsKey(og))
                  predictions.put(og, new HashSet<String>());
              
              predictions.get(og).add(tok[j].trim());
          }
      }
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}     
       
       System.out.println("predictions size: "+predictions.keySet().size());
       System.out.println("accretion size: "+accretions.keySet().size());
         
       double onlyKnownAnnotationAccretion=0.0, newlyPredictedAnnotationAccretion=0.0, knownAndPredictedAccretion = 0.0;
       double numGenes = 0.0; 
       
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         Random rand = new Random();
         
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap); //change to 0 to get COGs from file
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                   int taxid=taxIDMap.get(fileName);
                   
                   //if(taxid!=511145)
                   //    continue;
                  
                  //create accretionCounting
                  
                  for(int g=0;g<cg.ancogs.size();g++){
                      String gen = cg.ancogs.get(g).getValue2();
                      HashSet<String> ogs = geneOGMap.geneOGsMap.get(gen);
                      
                      HashSet<String> allgos = new HashSet<>();
                      HashSet<String> predictedgos = new HashSet<>();
                      
                      int useful = 0;
                      
                      for(String og:ogs){
                         ArrayList<String> go = cgmap.CogGOmap.get(og);
                         if(!cgmap.CogGOmap.containsKey(og))
                             continue;
                         
                         if(go.size()>0)
                             useful = 1;
                         
                         allgos.addAll(go);
                         if(predictions.containsKey(og))
                             predictedgos.addAll(predictions.get(og));
                      }
                      
                      if(useful==1)
                        numGenes++;
                      
                      for(String go:allgos){
                          if(!accretions.containsKey(go))
                              System.out.println("Error: "+go);
                          if(!predictedgos.contains(go)){
                              onlyKnownAnnotationAccretion+=accretions.get(go);
                          }
                          else if(predictedgos.contains(go)){
                              knownAndPredictedAccretion+=accretions.get(go);
                          }
                      }
                      
                      for(String go:predictedgos){
                          if(!allgos.contains(go)){
                              newlyPredictedAnnotationAccretion+=accretions.get(go);
                          }
                      }
                      
                  }
                 // numGenes+=cg.ancogs.size();

            }
            
            
            FileWriter fw = new FileWriter(output);
            
            fw.write("Only known: "+(onlyKnownAnnotationAccretion/numGenes)+"\n");
            fw.write("Known and predicted: "+(knownAndPredictedAccretion/numGenes)+"\n");
            fw.write("Only predicted: "+(newlyPredictedAnnotationAccretion/numGenes)+"\n");
            fw.close();
            
            System.out.println("Only known: "+(onlyKnownAnnotationAccretion/numGenes));
            System.out.println("Known and predicted: "+(knownAndPredictedAccretion/numGenes));
            System.out.println("Only predicted: "+(newlyPredictedAnnotationAccretion/numGenes));

        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
      
     void createPid2GeneMapping(String taxIDFilePath, String[] extensions, boolean recursive, GeneOGMapping pid2OGMap, HashMap<String, HashSet<String>> pid2Gene, int taxID){
         
          HashMap<Integer,Integer> taxIDTranslationMap= new HashMap<>();
         String taxIDTranslation = "TaxIDTranslationFile.txt";
        
        File translationFile=new File(taxIDTranslation);
        BufferedReader reader = null;
        
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
         
         
          try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);

         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                 File input = (File) iterator.next();
                 System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                 int taxid=taxIDMap.get(fileName);
                 
                 if(!taxIDTranslationMap.containsKey(taxid))
                     continue;
                 
                 if(taxIDTranslationMap.get(taxid) != taxID)
                     continue;
                 
                  COG cg=new COG();
                  cg.createPid2GeneMapping(input, pid2OGMap, pid2Gene);
                 
                  
            }
          }
            catch(Exception e){
                    e.printStackTrace();
                    }
         
          File output = new File("pid2geneMapping");
          
          try{
                FileWriter fw = new FileWriter(output);
                
                Iterator<String> it = pid2Gene.keySet().iterator();
                
                while(it.hasNext()){
                    String pid = it.next();
                    fw.write(pid+": ");
                    HashSet<String> s = pid2Gene.get(pid);
                    
                    for(String g:s)
                        fw.write(g+" ");
                    
                    fw.write("\n");
                }
                
                fw.close();
          }
          catch(IOException e){
              e.printStackTrace();
          }
          
     }

      void createBaselineOrganisms(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap,GeneOGMapping geneOGMap, String output,int isTrain, boolean randomize){
    
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
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         COGbaselinefeatures cbf=new COGbaselinefeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap); //change to 0 to get COGs from file
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                if(usedTIDs.contains(taxid)){
                    //csf.imputeSparseFeatures(cf,cgmap, tidOccurence.get(taxid));
                    cbf.createFeaturesOrganism(sim, cgmap,geneOGMap);
                }
                else{
                cbf.createFeatures(sim, cgmap,geneOGMap);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            cbf.normalize();
            tt.createTrainHeader(output, cgmap, isTrain);
            tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
         //tt.createSparseTrainHeader(output, cgmap, isTrain, numGenomes);
        // tt.appendSparseRowsToSet(output, cgmap, csf, isTrain, numGenomes);
         //csf.COGsparsefeaturesmap.clear();
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
                
       void createContingency(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output,int isTrain){
    
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
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);//nebitno
         cot.initializeContingency(cgmap);
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
         int batchSize=30, batchCount=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
         ArrayList<COG> COGlist=new ArrayList<>();
         final ArrayList<COGSimilarity1> simList=new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                
                /*if(batchCount==0){
                    for(int i=0;i<batchSize;i++){
                        COG cg=new COG();
                        COGlist.add(cg);
                        COGSimilarity1 sim=new COGSimilarity1();
                        simList.add(sim);
                    }
                        
                }*/
                
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);

                 /* if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }*/

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                simList.add(sim);
                batchCount++;
               // termination--;
                //if(termination==0)
                 //   break;
                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
               // final PriorityQueue qD=new PriorityQueue($trainingData.getNbRows()-1);
                final Lock lock = new ReentrantLock();
                
               // System.out.println("Starting parallel threads...");
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_num = index;
                        @Override
                        public void run() {
                          //  System.out.println("Computing in parallel threads...");
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                                //  System.out.println("dist1: "+dist1);
                                  //try{
                                   //     lock.lock();
                                    //    try{
                                          //q.addElement(o, dist1);  
                                        //  System.out.println("queue size: "+q.getSize());
                                     //   }finally{
                                     //       lock.unlock();
                                      //  }
                               //   }catch( Exception e ) {
                                 //   e.printStackTrace();
                                 //   }
                                        
                                       // lock.unlock();
                             
                        }
                          });
                    }
                    }catch(Exception e){
                        //System.out.println("Something went wrong...");
                        e.printStackTrace();
                    }
                finally {
                  //  System.out.println("Shutting down...");
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                    
                    //cot.addContingency(sim, cgmap);
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);

            }
            
            if(batchCount!=0){
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
                        //System.out.println("Something went wrong...");
                        e.printStackTrace();
                    }
                finally {
                  //  System.out.println("Shutting down...");
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
            }
            
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
             
         void createContingencySelection(String taxIDFilePath, String selectiontaxIdFilePath ,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output,int isTrain){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
        File selTaxIDFile=new File(selectiontaxIdFilePath);
       HashSet<Integer> selectedTaxa = new HashSet<>();
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
        path =Paths.get(selTaxIDFile.getAbsolutePath());
        reader = Files.newBufferedReader(path,ENCODING);
        
      while ((line = reader.readLine()) != null){
         selectedTaxa.add(Integer.parseInt(line.trim()));
      }
      
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);//nebitno
         cot.initializeContingency(cgmap);
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
         int batchSize=30, batchCount=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
         ArrayList<COG> COGlist=new ArrayList<>();
         final ArrayList<COGSimilarity1> simList=new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                
                /*if(batchCount==0){
                    for(int i=0;i<batchSize;i++){
                        COG cg=new COG();
                        COGlist.add(cg);
                        COGSimilarity1 sim=new COGSimilarity1();
                        simList.add(sim);
                    }
                        
                }*/
                
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                 int tid = taxIDMap.get(fileName);
                 
                 if(!selectedTaxa.contains(tid))
                     continue;

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);

                 /* if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }*/

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                simList.add(sim);
                batchCount++;
               // termination--;
                //if(termination==0)
                 //   break;
                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
               // final PriorityQueue qD=new PriorityQueue($trainingData.getNbRows()-1);
                final Lock lock = new ReentrantLock();
                
               // System.out.println("Starting parallel threads...");
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_num = index;
                        @Override
                        public void run() {
                          //  System.out.println("Computing in parallel threads...");
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                                //  System.out.println("dist1: "+dist1);
                                  //try{
                                   //     lock.lock();
                                    //    try{
                                          //q.addElement(o, dist1);  
                                        //  System.out.println("queue size: "+q.getSize());
                                     //   }finally{
                                     //       lock.unlock();
                                      //  }
                               //   }catch( Exception e ) {
                                 //   e.printStackTrace();
                                 //   }
                                        
                                       // lock.unlock();
                             
                        }
                          });
                    }
                    }catch(Exception e){
                        //System.out.println("Something went wrong...");
                        e.printStackTrace();
                    }
                finally {
                  //  System.out.println("Shutting down...");
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                    
                    //cot.addContingency(sim, cgmap);
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);

            }
            
            if(batchCount!=0){
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
                        //System.out.println("Something went wrong...");
                        e.printStackTrace();
                    }
                finally {
                  //  System.out.println("Shutting down...");
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
            }
            
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
             
     void COGneighFuncComp(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap,String output, String GO, int isTrain){

         
         BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         HashMap<Integer,Integer> NumPlasmids=new HashMap<>();
         
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
         
         try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         //TrainTest tt=new TrainTest();
         COGNeighFC cogNFC=new COGNeighFC();
         cogNFC.initCog(cgmap);
         cogNFC.initialize(cgmap, GO); 
         //COGbaselinefeatures cbf=new COGbaselinefeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
         int proba=1;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
               
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);

                 /* if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }*/

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                
               /* if(!usedTIDs.contains(taxid)){
                    numGenomes++;
                    System.out.println("Num genomes: "+numGenomes);
                }
                //use to count the number of usable genomes...
                if(proba==1)
                    continue;*/
                
                
                
                if(usedTIDs.contains(taxid)){
                    cogNFC.addComp(sim, cgmap,geneOGMap,GO, k,tidOccurence.get(taxid));
                    //csf.imputeSparseFeatures(cf,cgmap, tidOccurence.get(taxid));
                    //cbf.createFeaturesOrganism(sim, cgmap);
                    int v=NumPlasmids.get(tidOccurence.get(taxid));
                    v=v+1;
                    NumPlasmids.put(tidOccurence.get(taxid), v);
                }
                else{
                    
               // cbf.createFeatures(sim, cgmap);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                 tidOccurence.put(taxid, numGenomes);
                 cogNFC.addComp(sim, cgmap, geneOGMap, GO, k,tidOccurence.get(taxid));
                 NumPlasmids.put(numGenomes, 1);
                 numGenomes++;
                
                }
                
                //numGenomes++;
                System.out.println("Obradeno bakterija: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            
            
            for(int i=0;i<cogNFC.cogNFc.get(0).get(cgmap.GOtoIndex.get(GO)).size();i++){
                for(int j=0;j<cogNFC.CogIndex.keySet().size();j++){
                    double totalPerc=cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(i);
                    totalPerc/=NumPlasmids.get(i);
                    cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).set(i, totalPerc);
                            }
            }
            
            HashMap<Integer,ArrayList<Double>> finalStat=new HashMap<>();
            
           for(int j=0;j<cogNFC.CogIndex.keySet().size();j++){ 
               
               int c=0;
                for(int k1=0;k1<cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).size();k1++){
                    if(cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1)>0)
                         c++;
                }
                    
               
               double tmp[]=new double[c/*cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).size()*/];
               
               
               for(int k1=0;k1<tmp.length;k1++)
                   if(cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1)>0)
                   tmp[k1]=cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1);
               
               if(c==0){
                   tmp=new double[1];
                   tmp[0]=0.0;
               }
               
            Statistics stat=new Statistics(tmp);
            double std=stat.getStdDev();
            double avg=stat.getMean();
            //System.out.println("avg: "+avg);
            ArrayList<Double> tmSt=new ArrayList<>();
            tmSt.add(avg); tmSt.add(std);
            finalStat.put(j, tmSt);
           }
            
           // cogNFC.normalize(cgmap, GO, numGenomes);
            
            FileWriter fw = null;
           String fileTmp="COGNeighCompPerc"+GO+".txt";
        fw = new FileWriter(fileTmp); 
        fw.write("COG \t "); fw.write("mean \t "); fw.write(" std \n");
        
        
        for(int i=0;i<cogNFC.CogIndex.keySet().size();i++){
            ArrayList<Double> tmp=finalStat.get(i);
            fw.write(cogNFC.IndexCog.get(i)+" \t ");
            for(int itt=0;itt<tmp.size();itt++){
                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+" \t ");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
            }
        }
        fw.close();
        
        //for(int i=0;i<cogNFC.cogNFc.get(0).size();i++)
           // if((i+1)<cogNFC.cogNFc.get(0).size())
           // fw.write("FreqORg"+i+" \t ");
        //else
               // fw.write("FreqORg"+i+" \n ");
       // fw.write("NumOrg: \n");
        
       /* for(int i=0;i<cogNFC.CogIndex.keySet().size();i++){
                          
                ArrayList<Double> tmp=cogNFC.cogNFc.get(i).get(cgmap.GOtoIndex.get(GO));
                fw.write(cogNFC.IndexCog.get(i)+" \t ");
                            for(int itt=0;itt<tmp.size();itt++){
                                //if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+" \t ");
                                //else
                                  //  fw.write(tmp.get(itt)+"\n"); 
                                //if(itt+1==tmp.size())
                                //fw.write("\n");
                            }
                            fw.write(numGenomes+"\n");
        }
       fw.close();*/
           // cbf.normalize();
           // tt.createTrainHeader(output, cgmap, isTrain);
            //tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
         //tt.createSparseTrainHeader(output, cgmap, isTrain, numGenomes);
        // tt.appendSparseRowsToSet(output, cgmap, csf, isTrain, numGenomes);
         //csf.COGsparsefeaturesmap.clear();
            //cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
       
   }
        
     void COGneighFuncCompRand(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output, String GO, int isTrain){

         
         BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         HashMap<Integer,Integer> NumPlasmids=new HashMap<>();
         
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
         
         try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         //TrainTest tt=new TrainTest();
         COGNeighFC cogNFC=new COGNeighFC();
         cogNFC.initCog(cgmap);
         cogNFC.initialize(cgmap, GO); 
         //COGbaselinefeatures cbf=new COGbaselinefeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0, numFiles=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
         Random rand=new Random();
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);
                  cg.randomize(rand);

                 /* if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }*/

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                if(usedTIDs.contains(taxid)){
                    cogNFC.addComp(sim, cgmap, geneOGMap ,GO,k,tidOccurence.get(taxid));
                    //csf.imputeSparseFeatures(cf,cgmap, tidOccurence.get(taxid));
                    //cbf.createFeaturesOrganism(sim, cgmap);
                    int v=NumPlasmids.get(tidOccurence.get(taxid));
                    v=v+1;
                    NumPlasmids.put(tidOccurence.get(taxid), v);
                }
                else{
                    
               // cbf.createFeatures(sim, cgmap);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                 tidOccurence.put(taxid, numGenomes);
                 cogNFC.addComp(sim, cgmap, geneOGMap ,GO,k ,tidOccurence.get(taxid));
                 NumPlasmids.put(numGenomes, 1);
                 numGenomes++;
                
                }
                
                //numGenomes++;
                System.out.println("Obradeno bakterija: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            
            
            for(int i=0;i<cogNFC.cogNFc.get(0).get(cgmap.GOtoIndex.get(GO)).size();i++){
                for(int j=0;j<cogNFC.CogIndex.keySet().size();j++){
                    double totalPerc=cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(i);
                    totalPerc/=NumPlasmids.get(i);
                    cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).set(i, totalPerc);
                            }
            }
            
            HashMap<Integer,ArrayList<Double>> finalStat=new HashMap<>();
            
           /*for(int j=0;j<cogNFC.CogIndex.keySet().size();j++){ 
               double tmp[]=new double[cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).size()];
               
               for(int k1=0;k1<tmp.length;k1++)
                   tmp[k1]=cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1);
               
            Statistics stat=new Statistics(tmp);
            double std=stat.getStdDev();
            double avg=stat.getMean();
            //System.out.println("avg: "+avg);
            ArrayList<Double> tmSt=new ArrayList<>();
            tmSt.add(avg); tmSt.add(std);
            finalStat.put(j, tmSt);
           }*/
            
            for(int j=0;j<cogNFC.CogIndex.keySet().size();j++){ 
               
               int c=0;
                for(int k1=0;k1<cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).size();k1++){
                    if(cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1)>0)
                         c++;
                }
                    
               
               double tmp[]=new double[c/*cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).size()*/];
               
               
               for(int k1=0;k1<tmp.length;k1++)
                   if(cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1)>0)
                   tmp[k1]=cogNFC.cogNFc.get(j).get(cgmap.GOtoIndex.get(GO)).get(k1);
               
               if(c==0){
                   tmp=new double[1];
                   tmp[0]=0.0;
               }
               
            Statistics stat=new Statistics(tmp);
            double std=stat.getStdDev();
            double avg=stat.getMean();
            //System.out.println("avg: "+avg);
            ArrayList<Double> tmSt=new ArrayList<>();
            tmSt.add(avg); tmSt.add(std);
            finalStat.put(j, tmSt);
           }
            
           // cogNFC.normalize(cgmap, GO, numGenomes);
            
            FileWriter fw = null;
           String fileTmp="COGNeighCompPercRand"+GO+".txt";
        fw = new FileWriter(fileTmp); 
        fw.write("COG \t "); fw.write("mean \t "); fw.write(" std \n");
        
        
        for(int i=0;i<cogNFC.CogIndex.keySet().size();i++){
            ArrayList<Double> tmp=finalStat.get(i);
            fw.write(cogNFC.IndexCog.get(i)+" \t ");
            for(int itt=0;itt<tmp.size();itt++){
                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+" \t ");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
            }
        }
        fw.close();
        
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
       
       void createContingencyBaseline(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output,int isTrain){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
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
          
         Random rand=new Random();
         
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         cot.initializeContingency(cgmap);
         int numGenomes=0;
         int batchSize=30, batchCount=0;
         ArrayList<COG> COGlist=new ArrayList<>();
         final ArrayList<COGSimilarity1> simList=new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);
                  cg.randomize(rand);

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);

                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                simList.add(sim);
                batchCount++;

                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                             int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c, cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
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
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);

            }
            
            if(batchCount!=0){
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                             int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
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
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
            }
            
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }    
   } 
             
    void createContingencySelectionBaseline(String taxIDFilePath, String selectiontaxIdFilePath, String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output,int isTrain){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       File selTaxIDFile = new File(selectiontaxIdFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> selectedTaxa = new HashSet<>();
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
      
      path =Paths.get(selTaxIDFile.getAbsolutePath());
        reader = Files.newBufferedReader(path,ENCODING);
        
      while ((line = reader.readLine()) != null){
         selectedTaxa.add(Integer.parseInt(line.trim()));
      }
      
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
         Random rand=new Random();
         
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         cot.initializeContingency(cgmap);
         int numGenomes=0;
         int batchSize=30, batchCount=0;
         ArrayList<COG> COGlist=new ArrayList<>();
         final ArrayList<COGSimilarity1> simList=new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                 int tid = taxIDMap.get(fileName);
                 
                 if(!selectedTaxa.contains(tid))
                     continue;      

                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap);
                  cg.randomize(rand);

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);

                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                simList.add(sim);
                batchCount++;

                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                             int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c, cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
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
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);

            }
            
            if(batchCount!=0){
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final COGGOMap c=cgmap;
                    final GeneOGMapping cGO=geneOGMap;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                             int proc_num = index;
                        @Override
                        public void run() {
                            
                                  contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO);
                             
                        }
                          });
                    }
                    }catch(Exception e){
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
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
            }
            
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }    
   } 
                
       void createContingency(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap,GeneOGMapping geneOGMap ,String output, String COG, int isTrain){
    
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
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         //cot.initializeContingency(cgmap);
         cot.initializeContingency(cgmap, COG);
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                cot.addContingency(sim, cgmap, COG);
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
       
    void createContingencyPO(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap ,String output, String COG, int isTrain,  ArrayList<HashSet<String>> categories){
    
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
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         //cot.initializeContingency(cgmap);
         cot.initializeContingency(cgmap, COG,categories);
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                cot.addContingency(sim, cgmap, COG,categories);
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }   
   
        void createContingencyRandom(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap, String output, ArrayList<Pair<String,String>> GOPairs, int isTrain){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
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
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         cot.initializeContingencyPairs(cgmap, GOPairs);
         int numGenomes=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                cot.addContingencyRandom(sim, cgmap, GOPairs);     
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
            }
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
      
        void createContingencyBaseline(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap,GeneOGMapping geneOGMap ,String output, ArrayList<String> GOs, ArrayList<HashSet<String>> categories, int isTrain){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
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
          Random rand=new Random();
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);
         cot.initializeContingency(cgmap, GOs,categories);
         int numGenomes=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);
                  cg.randomize(rand);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap,k);
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                cot.addContingency(sim, cgmap, GOs, categories);     
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
            }
            cot.writeContingency(cgmap, GOs,output);
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
        
       void createOccurence(String taxIDFilePath,String[] extensions, boolean recursive, int k,COGGOMap cgmap, GeneOGMapping geneOGMap, String output){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         COGOccurence occ=new COGOccurence();
         occ.initialize(cgmap);
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
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         TrainTest tt=new TrainTest();
        // COGbaselinefeatures cbf=new COGbaselinefeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         //COGSparsefeatures csf=new COGSparsefeatures();
         int numGenomes=0;
        // tt.createSparseTrainHeader(output, cgmap, isTrain,1);
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  COG cg=new COG();
                  cg.findCogs(input,0,geneOGMap);

                  if(cg.ancogs.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      continue;
                  }

                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(cg,cgmap,geneOGMap, k);
                //COGfeatures cf=new COGfeatures();
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                if(usedTIDs.contains(taxid)){
                    //csf.imputeSparseFeatures(cf,cgmap, tidOccurence.get(taxid));
                    //cbf.createFeaturesOrganism(sim, cgmap);
                }
                else{
                //cbf.createFeatures(sim, cgmap);
                 occ.createOccurence(cg, cgmap);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                //csf.createSparseFeatures(cf, cgmap, numGenomes);
                cg.cogs.clear();
                cg.ancogs.clear();
                sim.neighbours.clear();
               // cf.COGfeaturesmap.clear();
                // tt.appendRowsToSet(output, cgmap, cf, isTrain);
            }
            //cbf.normalize();
            //tt.createTrainHeader(output, cgmap, isTrain);
            //tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
            occ.writeToFile(output);
         //tt.createSparseTrainHeader(output, cgmap, isTrain, numGenomes);
        // tt.appendSparseRowsToSet(output, cgmap, csf, isTrain, numGenomes);
         //csf.COGsparsefeaturesmap.clear();
            //cbf.COGbaselinefeaturesmap.clear();
            occ.occurence.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
        
   }