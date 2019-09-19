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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to count annotated genes and OGs
 */
public class CountAnotatedGenesAndOgs {
    public static void main(String args[]){
        
        HashMap<Integer,Integer> geneCount = new HashMap<>();//tax - num An genes
        HashMap<Integer,Integer> OgCount = new HashMap<>();//tax - num An ogs
        
         final HashMap<String,HashSet<Pair<String,String>>> geneOGMap=new HashMap<String,HashSet<Pair<String,String>>>();
              
        String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa
        
        //String winMPE= "id_conversionReducedEns.tsv";//fungy
         String winMPE= "id_conversionReducedEnsMz.tsv";//metazoa
       // String geneOGMapFile = "geneOGMappingFinalOKNewnonCH.txt";//fungy
           // String geneOGMapFile = "geneOGMappingFinalOKNew.txt";//metazoa
        String taxIDFilePath = "taxID.txt";
         Mappings map=new Mappings();
         // map.loadMappings(new File(winMPE)); 
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
       // File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
        OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);
       
        File predInput=null;//new File("DistanceLocLogD.train.1.pred.arff");//BaselineOGsk=4.train.1.pred.arff
        //File predInput=new File("BaselineOGsk=4.train.1.pred.arff");
       
          CreateReducedMappingsFile rmf = new CreateReducedMappingsFile();
          //rmf.loadNogMembersMappings(new File("fuNOG.members.tsv"));
          rmf.loadNogMembersMappingsMetazoa(new File("meNOG.members.tsv"));
          
        GenesAndOgsOrganism organismMappings = new GenesAndOgsOrganism();
        
           int test=1;
        
       
        //File root=new File("../FungiExtracted");//fungi
        File root= new File("/home/mmihelcic/MatejQnap/MetazoaExtracted"); //metazoa
         String extensions[] = {"dat"};
        Boolean recursive=true;
       
        for(int taxIdT: speciesTax){
        organismMappings.traverseGenomesAndCount(root, taxIdT, extensions, recursive, taxIDFilePath , geneOGMap, map, rmf,geneCount,OgCount);
        organismMappings.clearMappins();
        }
        
         String outputName = "CountsAnotOGsGenes.txt";
      String info="";
      
        try{ 
          FileWriter fw=new FileWriter(outputName);
          
          fw.write("TaxID\t#OGs\t#genes\t\n");
          fw.write("____________________\n\n");
         
      for(int txIdT:speciesTax){         
          fw.write(txIdT+"\t"+OgCount.get(txIdT)+"\t"+geneCount.get(txIdT)+"\n");
          
          
          fw.write("\n\n");
      }
          
          fw.close();
          
      }
         catch(Exception e){
             e.printStackTrace();
         }        
    }
}
