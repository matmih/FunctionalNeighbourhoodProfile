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
  *@description class to create final datasets to train classification models (eukaryotic organisms), it also contains functions to create datasets to compute accrition
 */
public class CreateFungyDatasetNew {
    
     public static void main(String [] args){
         
         int OrganismGroup = 0;//0 fungy, 1 metazoa
         int Algorithm = 0;//0 GFA, 1 Baseline, 4-GFP
         int accrAlgorithm = 0; //0-FNP, 1-Baseline, 2-GFP
         double precTresh = 0.0;
         
         if(args.length>0)
            OrganismGroup=Integer.parseInt(args[0]);
         
         if(args.length>1)
             Algorithm = Integer.parseInt(args[1]);
         
         if(args.length>2)
             accrAlgorithm = Integer.parseInt(args[2]);
         
         if(args.length>3)
             precTresh = Double.parseDouble(args[3]);
         
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
            String output = "FungyFA_3_10_30k=10BB.arff";//fungy
            String outputTrain = "FungyFATrain_3_10_30k=10BB.arff";//fungy
            String outputTest = "FungyFATest_3_10_30k=10BB.arff";//fungy
            String outputSinteny = ""; //only for metazoa 
            
            if(OrganismGroup==1){
               headerInput = new File("headerFungyGONew3_10_30New.txt");
               output = "MetazoaFA_3_10_30k=10New.arff";
               outputTrain = "MetazoaFATrain_3_10_30k=10New.arff";
               outputTest = "MetazoaFATest_3_10_30k=10New.arff";
               outputSinteny = "MetazoaFA_3_10_30k=10NewNoSB.arff";
            }
			
            ConstructionAlgorithms alg=new ConstructionAlgorithms("/home/mmihelcic/GOdatasetGenerator/Eukariots/FungiExtracted");//fungy - abacus2
           
            if(OrganismGroup==1){
                alg=new ConstructionAlgorithms("/home/mmihelcic/MatejQnap/MetazoaExtracted");// metazoa - ada
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
             
             int k=10,isTrain=1;
            
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
              String outputLocation = "FungyDist_3_10_30BB.arff";
              String outputGFP = "FungyGFPBB.arff";
              String outputTrainLocation = "FungyDistTrain_3_10_30BB.arff";
              String outputTestLocation="FungyDistTest_3_10_30BB.arff";
              if(OrganismGroup==1){
                  outputLocation = "MetazoaDist_3_10_30BB.arff";
                  outputGFP = "MetazoaGFPBB.arff";
              }

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
           //create metazoa dataset without synthenic blocks
          HashMap<String,Integer> taxIDMapSinteny = new HashMap<>();
          File inputTaxS = new File("taxIDsSB.txt");
          
          try{
              Path p = Paths.get(inputTaxS.getAbsolutePath());
              BufferedReader rt = Files.newBufferedReader(p,StandardCharsets.UTF_8);
              
              String line = "";
              
              while((line = rt.readLine())!=null){
                  
                  String tmp[] = line.split("\t");
                  taxIDMapSinteny.put(tmp[1].trim(), Integer.parseInt(tmp[0].trim()));
                  
              }
              
              rt.close();
              
          }
          catch(IOException e){
              e.printStackTrace();
          }
          
          System.out.println("TID map size: "+taxIDMapSinteny.keySet().size());
          
          //call the function to create sintenic block structure and create the function to eliminate the required blocks of genome
          File inputTaxFolder = new File("/home/mmihelcic/GOdatasetGenerator/Eukariots/MetazoaFiles/AllFiles");
          String ext1[] = {"txt"};
          SynthenicBlocks blocks = new SynthenicBlocks();
          blocks.readBlocks(inputTaxFolder, ext1, taxIDMapSinteny);
          System.out.println("Block size: "+blocks.blocks.keySet().size()+" "+blocks.blocks.size());
          Iterator<Integer> it = blocks.blocks.keySet().iterator();
          
          while(it.hasNext()){
              int tid = it.next();
              System.out.println("taxid: "+tid+" Num Blocks: "+blocks.blocks.get(tid).size());
          }
          //eliminate blocks from genome and create dataset
         alg.createBaselineOrganismsNewNoSintheny(taxIDFilePath, ensembleTIDsFile, blocks , rmf, extensions, true, k, oggomap, geneOGMap, ogOrgCount, ogOccCount, map, headerInput, outputSinteny, isTrain, false, taxIDTranslationMap, 1, OrganismGroup);
    
      }
      if(Algorithm == 7){
      
      alg.countGOFrequencyNew(taxIDFilePath, ensembleTIDsFile, rmf, extensions, true, oggomap, geneOGMap, map, false, taxIDTranslationMap, k, 1); 
      }
      if(Algorithm == 8){
          String accretionFilePath = "";
          output = "accretion";
          
          if(OrganismGroup == 0){
              accretionFilePath= "IACFungi.txt";
              output+="Fungi";
          }
          else{ 
              accretionFilePath = "IACMetazoa.txt";
              output+="Metazoa";
          }
              
          String predictionFilePath = "";//prediction at precision file
          
          if(OrganismGroup == 0){
              predictionFilePath = "";
              if(accrAlgorithm == 0){
                  predictionFilePath = "PredictedCOGGOPairs0_"+precTresh+"_FNPFungi.txt";
                  output+="FNP.txt";
              }
              else if(accrAlgorithm == 1){
                  predictionFilePath = "PredictedCOGGOPairs2_"+precTresh+"_10NNFungi.txt";
                  output+="10NN.txt";
              }
              else if(accrAlgorithm == 2){
                  predictionFilePath = "PredictedCOGGOPairs3_"+precTresh+"_GFPFungi.txt";
                  output+="GFP.txt";
              }
          }
          else{
              predictionFilePath = "";
              output+="Metazoa";
              
              if(accrAlgorithm == 0){
                  predictionFilePath = "PredictedCOGGOPairs0_"+precTresh+"_FNPMetazoa.txt";
                  output+="FNP.txt";
              }
              else if(accrAlgorithm == 1){
                  predictionFilePath = "PredictedCOGGOPairs2_"+precTresh+"_10NNMetazoa.txt";
                  output+="10NN.txt";
              }
              else if(accrAlgorithm == 2){
                  predictionFilePath = "PredictedCOGGOPairs3_"+precTresh+"_GFPMetazoa.txt";
                  output+="GFP.txt";
              }
              
          }
          
          int diffOnt =1;
          if(diffOnt == 0)         
            alg.computeAccrition(taxIDFilePath, accretionFilePath, predictionFilePath, extensions, true, k, OrganismGroup, map, oggomap, geneOGMap, taxIDTranslationMap, output, isTrain, false);
          else{
               OntologyTools.GeneOntology myGO=null;
 OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);
          ArrayList<HashSet<String>> categories=new ArrayList<>();
        GODevider gd=new GODevider();
        try{
            categories=gd.ExtractCategories(new File("go_daily-termdb.obo-xml.gz"),cgmap);
        }
        catch(Exception e){
            e.printStackTrace();
        }  
        System.out.println("Categories size: "+categories.size());
              alg.computeAccritionDiffOnt(taxIDFilePath, accretionFilePath, predictionFilePath, extensions, true, categories ,k, OrganismGroup, map, oggomap, geneOGMap, taxIDTranslationMap, output, isTrain, true);
          }
      }
    }
}
