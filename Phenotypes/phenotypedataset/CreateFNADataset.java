/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import static phenotypedataset.LoadGeneData.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create FNP phenotype dataset
 */
public class CreateFNADataset{
    public static void main(String [] args){
        Mappings map = new Mappings();
        
        File inputogFunc = new File("C:\\Users\\matej\\Documents\\Phenotypes\\ReducedCogOGMappingPhenotype.txt");
        File inputgeneOg = new File("C:\\Users\\matej\\Documents\\Phenotypes\\geneOgMapping.txt");
        
        String extension[] = {"gff"};
        File root = new File("C:\\Users\\matej\\Documents\\Phenotypes\\Bacteria");
        
        String pathStrains = "C:\\Users\\matej\\Documents\\Phenotypes\\BacteriaIntersection.txt";
        String pathPhenotypes = "C:\\Users\\matej\\Documents\\Phenotypes\\PhenotypeIntersection.txt";
        
        map.loadGeneFunctionMapping(inputogFunc);
        map.loadGeneOGMapping(inputgeneOg);    
        
        String output = "FunctionalNeighbourhoodFull.arff";
        
        //load targets
        
        HashMap<String,ArrayList<Integer>> bacterialTargets = new HashMap<>();
        
        File inputTargets = new File("C:\\Users\\matej\\Documents\\Phenotypes\\Targets.txt");
        
          BufferedReader reader;
         try {
      Path path =Paths.get(inputTargets.getAbsolutePath());
      System.out.println("Path: "+inputTargets.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      while ((line = reader.readLine()) != null) {
        String tmp[] = line.split(",");
        String  bacteria = tmp[0].trim();
        bacterialTargets.put(bacteria, new ArrayList<>());
        
        for(int i = 1; i<tmp.length;i++)
            bacterialTargets.get(bacteria).add(Integer.parseInt(tmp[i]));
         
      }
         }
         catch(Exception e){
             e.printStackTrace();
         }
        
        
        int k=4;
        
               try{
         Collection files = FileUtils.listFiles(root, extension, true);
         TrainTest tt=new TrainTest(output);
       
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         
         PhenotypesAndStrains selection = new PhenotypesAndStrains();
         selection.loadPhenotypesAndStrains(pathStrains, pathPhenotypes);
         
           

         int numGenomes=0, numFiles=0;

         tt.createTrainHeader(selection,map);
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                FunctionalFeatures cbf=new FunctionalFeatures(map);
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 String tmp[] = input.getName().split("_");
                 String bactName = tmp[0].trim();
                 
                 System.out.println(bactName);
                 
                 if(!selection.strains.contains(bactName))
                     continue;
                 
                 System.out.println("Loading genes");
                 LoadGeneData1 genes = new LoadGeneData1();
                 genes.loadData(input);

                 System.out.println("NCont: "+genes.geneContigs.size());
                 System.out.println("FCont: "+genes.geneContigs.get(0).size());
                  COGSimilarity1 sim=new COGSimilarity1();
                sim.computeSpatialNeighboors(genes,map,k);

                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                String fileName=input.getName().substring(0,input.getName().length()-4);
                    
                    cbf.createFeatures(sim, map);
                    tt.writeBact(bactName);
                    tt.createTrainBody(cbf,bacterialTargets,map,bactName);

                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
            }
           
          tt.Close();
        } catch (Exception e) {
            e.printStackTrace();
        }   
    }
}
