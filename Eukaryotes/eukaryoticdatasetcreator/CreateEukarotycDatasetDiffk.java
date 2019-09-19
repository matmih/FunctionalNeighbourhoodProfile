/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create eukaryotic datasets with different neighbourhood size k
 */
public class CreateEukarotycDatasetDiffk {
     public static void main(String [] args){
         
         int OrganismGroup = 0;//0 fungy, 1 metazoa
         int Algorithm = 0;//0 GFA, 1 Baseline, 4-GFP
         int k=10;
         int stroj = 0;//flag used to choose different machine, can be deleted
         
         if(args.length>0)
            OrganismGroup=Integer.parseInt(args[0]);
         
         if(args.length>1)
             Algorithm = Integer.parseInt(args[1]);
         
         if(args.length>2)
              k = Integer.parseInt(args[2]);
         
         if(args.length>3)
             stroj = Integer.parseInt(args[3]);
         
        File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        if(OrganismGroup == 1)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        OGGOMapping oggomap=new OGGOMapping();
        
        oggomap.createCOGGOMapping(ogFuncFile);
        System.out.println("Num functions: "+oggomap.GOtoIndex.keySet().size());
        System.out.println("Num ogs: "+oggomap.CogGOmap.keySet().size());
        
        HashSet<String> ogSet = new HashSet();
        HashMap<String,HashSet<Integer>> ogOrgCount = new HashMap<>();
        HashMap<String,HashMap<Integer,Integer>> ogOccCount = new HashMap<>();
        
        Iterator<String> itog = oggomap.CogGOmap.keySet().iterator();
        
        while(itog.hasNext()){
            String og = itog.next();
            ogSet.add(og);
            ogOrgCount.put(og, new HashSet<Integer>());
            ogOccCount.put(og, new HashMap<Integer,Integer>());
        }
        
         String extensions[] = {"dat"};
           String winMPE= "id_conversionReducedEns.tsv";//fungy
           if(OrganismGroup==1)
               winMPE= "id_conversionReducedEnsMz.tsv";
          
          Mappings map=new Mappings();
          if(OrganismGroup==0)
           map.loadMappings(new File(winMPE)); //fungy
          if(OrganismGroup==1)
            map.loadMappingsMetazoa(new File(winMPE)); 
            
             File headerInput = new File("headerFungyGONew3_10_30NewnonCH.txt");//fungy
            String output = "FungyFA_3_10_30k="+k+"BB.arff";//fungy
            String outputTrain = "FungyFATrain_3_10_30k=10BB.arff";//fungy
            String outputTest = "FungyFATest_3_10_30k=10BB.arff";//fungy
            String outputSinteny = ""; //only for metazoa 
            
            if(OrganismGroup==1){
               headerInput = new File("headerFungyGONew3_10_30New.txt");
               output = "MetazoaFA_3_10_30k="+k+"New.arff";
               outputTrain = "MetazoaFATrain_3_10_30k=10New.arff";
               outputTest = "MetazoaFATest_3_10_30k=10New.arff";
               outputSinteny = "MetazoaFA_3_10_30k=10NewNoSB.arff";
            }

            ConstructionAlgorithms alg = null;
            
            if(stroj == 0 && OrganismGroup == 0){
                  alg=new ConstructionAlgorithms("/home/mmihelcic/GOdatasetGenerator/Eukariots/FungiExtracted");//fungy - abacus2
            }
            else if(stroj == 0 && OrganismGroup == 1){
                alg=new ConstructionAlgorithms("/home/mmihelcic/MatejQnap/MetazoaExtracted");//metazoa
            }
            else if(stroj == 1 && OrganismGroup == 0){
               alg=new ConstructionAlgorithms("FungiExtracted");
            }
            else if(stroj == 1 && OrganismGroup == 1){
                 alg=new ConstructionAlgorithms("/home/mmihelcic/MatejQnap/MetazoaExtracted");
            }

            HashMap<String,HashSet<Pair<String,String>>> geneOGMap = new HashMap<>();
            String geneOGMapFile = "geneOGMappingFinalOKNewnonCH.txt";//fungy
            if(OrganismGroup==1)
                geneOGMapFile = "geneOGMappingFinalOKNew.txt";
            String taxIDFilePath = "taxID.txt";
            
            if(OrganismGroup == 0)
                taxIDFilePath = "taxIDFungi.txt";
            else taxIDFilePath = "taxIDMetazoa.txt";
            
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
             
             System.out.println("Num genes: "+geneOGMap.keySet().size());
             
             int isTrain=1;
            
              HashMap<Integer,Integer> taxIDTranslationMap= new HashMap<>();
            
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
             
             File ensembleTIDsFile = new File("ensembleTIDs.txt");
             CreateReducedMappingsFile rmf = new CreateReducedMappingsFile();

             if(OrganismGroup==0)
                 rmf.loadNogMembersMappings(new File("fuNOG.members.tsv"));//fungy
              if(OrganismGroup==1)
                  rmf.loadNogMembersMappingsMetazoa(new File("meNOG.members.tsv"));
             
              int useLeaveOut = 0;
              
              HashSet<String> leaveoutOGs = new HashSet<>();
              
              File inputLeaveOutOGs = new File("fungiInfo.txt");
              
              if(OrganismGroup==1)
                  inputLeaveOutOGs = new File("metazoaInfo.txt");
             
              if(useLeaveOut == 1){
              try{
                  Path p = Paths.get(inputLeaveOutOGs.getAbsolutePath());
                  BufferedReader buf = Files.newBufferedReader(p,StandardCharsets.UTF_8);
                  
                  String line = "";
                  
                  while((line=buf.readLine())!=null){
                        String tmp[] = line.split("\t");
                        String og = tmp[2].trim();
                        leaveoutOGs.add(og);
                  }
                  
              }
              catch(IOException e){
                  e.printStackTrace();
              }
             }
       
    if(Algorithm==0 || Algorithm ==3)
         alg.createBaselineOrganismsNew(taxIDFilePath, ensembleTIDsFile, rmf, extensions, true, k, oggomap, geneOGMap, ogOrgCount, ogOccCount, map, headerInput, output, isTrain, false, taxIDTranslationMap, 1, OrganismGroup);
   
              String outputCount = "OGOccurencenonCH.txt";
              String outputLocation = "FungyDist_3_10_30BBk="+k+".arff";
              String outputGFP = "FungyGFPBBk="+k+".arff";
              String outputTrainLocation = "FungyDistTrain_3_10_30BB.arff";
              String outputTestLocation="FungyDistTest_3_10_30BB.arff";
              if(OrganismGroup==1){
                  outputLocation = "MetazoaDist_3_10_30BBk="+k+".arff";
                  outputGFP = "MetazoaGFPBBk="+k+".arff";
              }

      if(Algorithm ==1 || Algorithm ==3)
              alg.createDistanceForest(taxIDFilePath, ensembleTIDsFile, rmf,extensions, true, k,1 ,oggomap ,geneOGMap, ogOrgCount,ogOccCount, map, headerInput, outputLocation, isTrain, false,taxIDTranslationMap,1, OrganismGroup);
       
      if(Algorithm == 4){
          alg.createAssociationGraph(taxIDFilePath, ensembleTIDsFile, rmf, extensions, true, k, 1, oggomap, geneOGMap, ogOrgCount, ogOccCount, map, headerInput, outputGFP, isTrain, false, taxIDTranslationMap, 1, OrganismGroup);
      }
      
      if(Algorithm == 5){
          String outputContingency = "ContingencyOut";
          boolean randomize = false;
           
          if(randomize == true)
              outputContingency = "ContingecyOutRand";
                  
          alg.createContingencyNoSinc(taxIDFilePath, extensions, true, k, oggomap, geneOGMap, map, outputContingency, isTrain, randomize, OrganismGroup, taxIDTranslationMap);
      }
      
      if(Algorithm == 6){
    
      }
      if(Algorithm == 7){
      
      alg.countGOFrequencyNew(taxIDFilePath, ensembleTIDsFile, rmf, extensions, true, oggomap, geneOGMap, map, false, taxIDTranslationMap, k, 1); 
      }
    }
}
