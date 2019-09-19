/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description helper class to count the number of annotated genes per strain in the data
 */
public class CountAnnotatedGenesPerStrain {
    static public void main(String [] args){
        
         Mappings map = new Mappings();
        
        File inputogFunc = new File("C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\ReducedCogOGMappingPhenotype.txt");
        File inputgeneOg = new File("C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\geneOgMapping.txt");
        
        String extension[] = {"gff"};
        File root = new File("C:\\Users\\mmihelcic\\Downloads\\Bacterial strains genomes\\Bacteria");
        
        String pathStrains = "C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\BacteriaIntersection.txt";
        String pathPhenotypes = "C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\PhenotypeIntersection.txt";
        
        map.loadGeneFunctionMapping(inputogFunc);
        map.loadGeneOGMapping(inputgeneOg);    
        
         String output = "AnnotatedGenes.txt";
         
          int k=4;
        
               try{
         Collection files = FileUtils.listFiles(root, extension, true);
         TrainTest tt=new TrainTest(output);
       
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         
         PhenotypesAndStrains selection = new PhenotypesAndStrains();
         selection.loadPhenotypesAndStrains(pathStrains, pathPhenotypes);

         int numGenomes=0, numFiles=0;
         
         FileWriter  fw = new FileWriter(output);
         
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
                  sim.countAnnotatedGenes(genes, map);

                  fw.write(bactName+"\t"+sim.numgenes+"\t"+sim.annotatedGenes+"\n");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
            }
            fw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }   
        
    }
        
}
