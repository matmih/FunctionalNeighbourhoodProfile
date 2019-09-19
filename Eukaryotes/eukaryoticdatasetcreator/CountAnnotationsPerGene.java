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
  *@description class to count the number of annotations per gene in eukaryotic datasets
 */
public class CountAnnotationsPerGene {
    
     public static void main(String [] args){
         
         int OrganismGroup = 0;//0 fungy, 1 metazoa
         int Algorithm = 0;//0 GFA, 1 Baseline, 4-GFP
         
         if(args.length>0)
            OrganismGroup=Integer.parseInt(args[0]);
         
         if(args.length>1)
             Algorithm = Integer.parseInt(args[1]);
          
       // File ogFuncFile=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_1_30F.txt");
        File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        if(OrganismGroup == 1)
            ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        //File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
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
          //String winMPE= "C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversionReducedEns.tsv";
           String winMPE= "id_conversionReducedEns.tsv";//fungy
           if(OrganismGroup==1)
               winMPE= "id_conversionReducedEnsMz.tsv";
             //String winMPE= "id_conversionReducedEnsMz.tsv";//metazoa
          
          Mappings map=new Mappings();
          if(OrganismGroup==0)
           map.loadMappings(new File(winMPE)); //fungy
          if(OrganismGroup==1)
            map.loadMappingsMetazoa(new File(winMPE)); 
            
            /*File headerInput = new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\headerFungyGONew3_1_30.txt");
            String output = "C:\\Users\\matej\\Downloads\\Eukaryot data\\FungyFA_3_1_30k=10.arff";*/
            
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
            
           // File headerInput = new File("headerFungyGONew3_10_30New.txt");//metazoa
           // String output = "FungyFA_3_10_30k=10NewMetazoa.arff";//metazoa
            
            /*ConstructionAlgorithms alg=new ConstructionAlgorithms("C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiExtracted");
            HashMap<String,HashSet<String>> geneOGMap = new HashMap<>();
            String geneOGMapFile = "C:\\Users\\matej\\Downloads\\Eukaryot data\\geneOGMappingFinalOK.txt";
            String taxIDFilePath = "C:\\Users\\matej\\Downloads\\Eukaryot data\\taxID.txt";*/
            
           // ConstructionAlgorithms alg=new ConstructionAlgorithms("FungiExtracted");
            ConstructionAlgorithms alg=new ConstructionAlgorithms("../FungiExtracted");//fungy
            
            if(OrganismGroup==1){
                alg=new ConstructionAlgorithms("/home/mmihelcic/MatejQnap/MetazoaExtracted");
            }
            
            //ConstructionAlgorithms alg=new ConstructionAlgorithms("/home/mmihelcic/MatejQnap/MetazoaExtracted");//metazoa
            
            HashMap<String,HashSet<Pair<String,String>>> geneOGMap = new HashMap<>();
            String geneOGMapFile = "geneOGMappingFinalOKNewnonCH.txt";//fungy
            if(OrganismGroup==1)
                geneOGMapFile = "geneOGMappingFinalOKNew.txt";
           // String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa
            String taxIDFilePath = "taxID.txt";
            
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
             
             File ensembleTIDsFile = new File("ensembleTIDs.txt");
             CreateReducedMappingsFile rmf = new CreateReducedMappingsFile();
            // rmf.loadEggnogProteins(new File("fuNOG.members.tsv"));
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

              String outputLocation = "FungyDist_3_10_30BB.arff";
              String outputGFP = "FungyGFPBB.arff";

              if(OrganismGroup==1){
                  outputLocation = "MetazoaDist_3_10_30BB.arff";
                  outputGFP = "MetazoaGFPBB.arff";
              }
              
              
              
              
          ArrayList<Double> numGenes = alg.numGenes(taxIDFilePath, ensembleTIDsFile, rmf, extensions, true, k, 1, oggomap, geneOGMap, ogOrgCount, ogOccCount, map, headerInput, outputGFP, isTrain, false, taxIDTranslationMap, 1, OrganismGroup);
          double numAnnots = 0.0;
          
          Iterator<String> it = oggomap.CogGOmap.keySet().iterator();
              while(it.hasNext()){
                  String og = it.next();
                  numAnnots+=oggomap.CogGOmap.get(og).size();               
              }
              
              
              System.out.println("NumGenes: "+numGenes.get(0));
              System.out.println("NumAnnots: "+numGenes.get(1));
              System.out.println("Annots per gene: "+(numGenes.get(1)/numGenes.get(0)));   
    }
}
