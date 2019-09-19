/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing a set of functions to compute different mappings/statistics
  *corresponding functions need to be uncommented to be used
 */
public class EukaryoticDatasetCreatorNew {
    
     public static void main(String[] args) {        
        CreateReducedMappingsFile rmf=new CreateReducedMappingsFile();
        
        String winEgMap="C:\\Users\\matej\\Downloads\\Eukaryot data\\fuNOG.members.tsv\\fuNOG.members.tsv";
        String idConv="C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversion.tsv";
        String winMp= "C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversionReduced.tsv";
        String winMPE= "C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversionReducedEns.tsv";
        String fungiGenomes="C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiExtracted";
        String geneIdFG= "C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiGenes.txt";
        String geneFunctionMapping = "C:\\Users\\matej\\Downloads\\Eukaryot data\\GeneFunctionsNew.txt";
        String uniprotMapping = "C:\\Users\\matej\\Downloads\\Eukaryot data\\goa_uniprot_all.gaf\\goa_uniprot_allReduced.gaf";
        String nogFunctionMapping = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctions0.3New.txt";
        String nogFunctionMappingFiltered = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsFiltered50_50nonCH.txt";
        String nogFunctionMappingProp = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_10_30New.txt";
        String nogFunctionMappingFinal = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsFinal3_30NewnonCH.txt";
        String nogFunctionMappingGORed = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsGORed3_10_30New.txt";
        String nogAtLeast20 = "C:\\Users\\matej\\Downloads\\Eukaryot data\\NogOccAtLeast3New.txt";//created in R
        String goCount = "C:\\Users\\matej\\Downloads\\Eukaryot data\\goCount.txt";
        String goCountRed = "C:\\Users\\matej\\Downloads\\Eukaryot data\\goReduced3_10_30New.txt";
        String nogCount = "C:\\Users\\matej\\Downloads\\Eukaryot data\\nogCountNew.txt";
        String functionOntology = "C:\\Users\\matej\\Downloads\\Eukaryot data\\go_daily-termdb.obo-xml\\go_daily-termdb.obo-xml";
        String geneOGMap = "C:\\Users\\matej\\Downloads\\Eukaryot data\\geneOGMappingFinalOKNew.txt";
        String taxIDTranslation = "C:\\Users\\matej\\Downloads\\Eukaryot data\\TaxIDTranslationFile.txt";
        String taxIDFilePath = "C:\\Users\\matej\\Downloads\\Eukaryot data\\taxID.txt";
        File ensembleTIDsFile = new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\ensembleTIDs.txt");
         String ensTaxs = "C:\\Users\\matej\\Downloads\\Eukaryot data\\ensembleTIDs.txt";
        
        
        /*String winEgMapMet="C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\meNOG.members.tsv\\meNOG.members.tsv";//windows
         String metazoaGenomes="C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\MetazoaExtracted";
         String geneIdMZ= "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\MetazoaGenes.txt";
         String winMPEMZ= "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\id_conversionReducedEnsMz.tsv";
         String geneFunctionMappingMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\GeneFunctionsNew.txt";
         String geneOGMapMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\geneOGMappingFinalOKNew.txt";*/
         
         
         /*String winEgMapMet="meNOG.members.tsv";//linux
         String metazoaGenomes="/home/mmihelcic/MatejQnap/MetazoaExtracted";//linux
         String geneIdMZ= "MetazoaGenes.txt";//linux
         String winMPEMZ= "id_conversionReducedEnsMz.tsv";//linux
         String geneFunctionMappingMet = "GeneFunctionsNew.txt";//linux
         String geneOGMapMet = "geneOGMappingFinalOKNew.txt";//linux
         String taxIDTranslationMet = "TaxIDTranslationFile.txt";//linux
         String nogFunctionMappingMet = "NogFunctions0.3New.txt";//linux*/
         
         
         String winEgMapMet="C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\meNOG.members.tsv\\meNOG.members.tsv";//win
         String metazoaGenomes="C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\MetazoaExtracted";//win
         String geneIdMZ= "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\MetazoaGenes.txt";//win
         String winMPEMZ= "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\id_conversionReducedEnsMz.tsv";//win
         String geneFunctionMappingMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\GeneFunctionsNew.txt";//win
         String geneOGMapMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\geneOGMappingFinalOKNew.txt";//win
         String taxIDTranslationMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\TaxIDTranslationFile.txt";//win
         String nogFunctionMappingMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctions0.3New.txt";//win
         String nogCountMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\nogCountNew.txt";//win
         String taxIDFilePathMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\taxID.txt";
         String nogAtLeast3Met = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogOccAtLeast3New.txt";//created in R
         String nogFunctionMappingFilteredMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsFiltered50_50nonCH.txt";
         String nogFunctionMappingPropMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsProp3_10_30New.txt";
         String nogFunctionMappingFinalMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsFinal3_30NewnonCH.txt";
         String nogFunctionMappingGORedMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsGORed3_10_30New.txt";
         String goCountRedMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\goReduced3_10_30New.txt";
         String ensTaxsMet = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\ensembleTIDs.txt";
         
        int test=1, useNonChromosomal = 1;
        
        /*rmf.filterNOGs(new File(nogFunctionMappingMet), new File(nogFunctionMappingFilteredMet));
        rmf.finalNOGSelection(new File(nogFunctionMappingFilteredMet), new File(nogAtLeast3Met), new File(nogFunctionMappingFinalMet));
        rmf.countGOFunctions(new File(nogFunctionMappingFinalMet), new File(goCountRedMet));*/
        rmf.filterAndReduceGOs(new File(nogFunctionMappingFinalMet), new File(goCountRedMet), new File(nogFunctionMappingGORedMet));
       rmf.propagateFunctions(new File(nogFunctionMappingGORedMet), functionOntology, new File(nogFunctionMappingPropMet));
       if(test==1)
            return;
        
        HashMap<Integer,Integer> taxIDTranslationMap= new HashMap<>();
        
        BufferedReader reader;
            
        
        File translationFile=new File(taxIDTranslationMet);
        
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
 
       // rmf.loadEggnogProteins(new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\fuNOG.members.tsv\\fuNOG.members.tsv"));//fungy
        //rmf.loadEggnogProteinsMetazoa(new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\meNOG.members.tsv\\meNOG.members.tsv"));//metazoa
       // rmf.countMappingsEnsemble(new File(idConv), new File(geneIdFG));
      //  rmf.writeGenesEnsemble(new File(fungiGenomes), new File(geneIdFG)); fungy
       // rmf.writeGenesEnsembleMetazoa(new File(metazoaGenomes), new File(geneIdMZ)); //metazoa
        //rmf.createReduced_mappingsEns(new File(idConv), new File(geneIdMZ), new File(winMPEMZ));//metazoa
        int t=1;
        /*        if(t==1)
                    return;*/

       // rmf.createReduced_mappings(new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversion.tsv"), new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversionReducedEx.tsv"));
       // rmf.createReduced_uniprot(new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\goa_uniprot_all.gaf\\goa_uniprot_all.gaf"), new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\goa_uniprot_all.gaf\\goa_uniprot_allReduced.gaf"));
        //rmf.writeGenesEnsemble(new File(fungiGenomes), new File(geneIdFG));
        //rmf.countMappingsSelected(new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversion.tsv"), new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiExtracted"), new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\eggnog4.protein_id_conversion.tsv\\id_conversionReducedEns.tsv"));
        
        HashSet<String> OKNOgs=new HashSet<>();
        
        //OKNOgs=rmf.loadOKNOgs(new File(nogFunctionMappingProp));
        //  rmf.loadNogMembersMappingsReduced(new File(winEgMap),OKNOgs);
       // rmf.loadNogMembersMappings(new File(winEgMap));// fungy
         rmf.loadNogMembersMappingsMetazoa(new File(winEgMapMet));//metazoa
         
        // rmf.writeGeneFunctionMappingNew(new File(fungiGenomes), new File(winMPE), new File(uniprotMapping), useNonChromosomal ,new File(geneFunctionMapping));//fungy
        // rmf.writeGeneFunctionMappingNewMetazoa(new File(metazoaGenomes), new File(winMPEMZ), new File(uniprotMapping), useNonChromosomal ,new File(geneFunctionMappingMet)); // metazoa
//        rmf.nonChromosomalTaxIDs(new File(metazoaGenomes), taxIDTranslationMap, new File(ensTaxsMet)); //metazoa
       
         
        //rmf.createGeneOGMappingNew(new File(fungiGenomes), new File(winMPE), taxIDTranslationMap, useNonChromosomal ,new File(geneOGMap)); //fungy
        //rmf.createGeneOGMappingNewMetazoa(new File(metazoaGenomes), new File(winMPEMZ), taxIDTranslationMap, useNonChromosomal ,new File(geneOGMapMet)); //metazoa
        /*rmf.AssignFunctionToNogsNewNew(new File(geneFunctionMappingMet),new File(geneOGMapMet), 0.3,new File(nogFunctionMappingMet));//metazoa
          if(t==1)
            return;*/
         
         HashMap<String,HashSet<Pair<String,String>>> geneOGMapM1 = new HashMap<>();
            //String geneOGMapFile = "geneOGMappingFinalOKNew.txt";
            
             try {
                     Path path =Paths.get(new File(geneOGMapMet).getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOGMapM1.containsKey(gene)){
                                geneOGMapM1.put(gene, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> ogs=geneOGMapM1.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og[]=tmp[i].split(":");
                               Pair<String, String> np = new Pair(og[0].trim(),og[1].trim());
                                
                               ogs.add(np);
                            }
                            geneOGMapM1.put(gene, ogs);
                    }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
         
         
         
          rmf.countNOGOccurenceNew(new File(metazoaGenomes), taxIDFilePathMet ,new File(geneOGMapMet), ensembleTIDsFile, taxIDTranslationMap, rmf, geneOGMapM1, useNonChromosomal, new File(nogCountMet));
       
        if(t==1)
            return;
        
        //rmf.countGOOccurence(new File(geneFunctionMapping), 0.3, new File(goCount));
       // rmf.countNOGOccurence(new File(fungiGenomes), new File(geneFunctionMapping), 0.3, new File(winMPE), new File(nogCount));
        
        
        /* File ogFuncFile=new File(nogFunctionMapping);
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
          
          Mappings map=new Mappings();
          */
        
                     if(useNonChromosomal == 1){
                 geneOGMap = geneOGMap.substring(0, geneOGMap.length()-4);
                 geneOGMap+="nonCH.txt";
                 nogCount = nogCount.substring(0, nogCount.length()-4);
                 nogCount+="nonCH.txt";
                 geneFunctionMapping = geneFunctionMapping.substring(0, geneFunctionMapping.length()-4);
                 geneFunctionMapping+="nonCH.txt";
                 nogFunctionMapping = nogFunctionMapping.substring(0,nogFunctionMapping.length()-4);
                 nogFunctionMapping+="nonCH.txt";
                 
             }
          
          
          HashMap<String,HashSet<Pair<String,String>>> geneOGMapM = new HashMap<>();
            //String geneOGMapFile = "geneOGMappingFinalOKNew.txt";
            
             try {
                     Path path =Paths.get(new File(geneOGMap).getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOGMapM.containsKey(gene)){
                                geneOGMapM.put(gene, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> ogs=geneOGMapM.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og[]=tmp[i].split(":");
                               Pair<String, String> np = new Pair(og[0].trim(),og[1].trim());
                                
                               ogs.add(np);
                            }
                            geneOGMapM.put(gene, ogs);
                    }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
          
        File ogOcc1=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\nogCountNewnonCHSC.txt");
             
       //  rmf.countGOOccWeighted(ogOcc1,new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_10_30NewF.txt"));
       //rmf.countNOGOccurenceNew(new File(fungiGenomes), taxIDFilePath ,new File(geneOGMap), ensembleTIDsFile ,taxIDTranslationMap, rmf, geneOGMapM, useNonChromosomal ,new File(nogCount));
       
       // File inputFolder, File geneFunctionFile, double perc, File mappingFile ,File outputFile
        
        //rmf.createGeneOGMappingNew(new File(fungiGenomes), new File(winMPE), taxIDTranslationMap, useNonChromosomal ,new File(geneOGMap));
       //rmf.writeGeneFunctionMappingNew(new File(fungiGenomes), new File(winMPE), new File(uniprotMapping), useNonChromosomal ,new File(geneFunctionMapping));
       //rmf.AssignFunctionToNogsNewNew(new File(geneFunctionMapping),new File(geneOGMap), 0.3,new File(nogFunctionMapping));
         //rmf.AssignFunctionToNogs(new File(geneFunctionMapping), 0.3, new File(nogFunctionMapping)); //old function
       // HashSet<Integer> res=rmf.countEggnogTax(new File(winEgMap), taxIDTranslationMap);
      //  rmf.countTaxSpecies(new File(fungiGenomes), taxIDTranslationMap, res);
        
        //rmf.countCoveredLocationsNew(new File(fungiGenomes), new File(geneOGMap), taxIDTranslationMap ,new File(winMPE));
        //rmf.countNOGOccurence(new File(fungiGenomes), new File(geneOGMap), new File(nogCount));
        
       
             

       // rmf.nonChromosomalTaxIDs(new File(fungiGenomes), taxIDTranslationMap, new File(ensTaxs));  
    }  
    
}
