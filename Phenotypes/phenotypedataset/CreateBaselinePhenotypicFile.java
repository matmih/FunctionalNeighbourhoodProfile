/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import static phenotypedataset.LoadGeneData.ENCODING;
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
import org.apache.commons.io.FileUtils;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create a conditional score file
 */
public class CreateBaselinePhenotypicFile {
 
    public static void main(String [] args){
        
        String inputCondScore = "C:\\Users\\mmihelcic\\Downloads\\elife-31035-supp4-v2\\conditional_scores.tsv";
        String inputGrowDefect = "C:\\Users\\mmihelcic\\Downloads\\elife-31035-supp5-v2\\phenotypic_data.tsv";
          String genomesPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains genomes\\Bacteria";
        String sequencesPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains sequences\\Sequences";
        String outputPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains sequences\\bactSeq.fa";
         String extension[] = {"gff"};
         String extension1[] = {"fasta"};
       
        File inputGenomes = new File(genomesPath);
        File inputSequences = new File(sequencesPath);
        String output = "BaselinePhenotype.arff";
        String outputPhenotypes = "PhenotypeIntersection.txt";
        String outputBacteria = "BacteriaIntersection.txt";
        
        File cScore = new File(inputCondScore);
        File cGD = new File(inputGrowDefect);
        
        HashMap<String,ArrayList<Double>> condS = new HashMap<>();
        HashMap<String,ArrayList<Integer>> growD = new HashMap<>();
        HashMap<String,Pair<String,Integer>> growDTemp = new HashMap<>();
        
        ArrayList<String> phenoDTemp = new ArrayList<>();
        ArrayList<String> bactDTemp = new ArrayList<>();
        ArrayList<Integer> deffDTemp = new ArrayList<>();
        
        ArrayList<String> bacteria = new ArrayList<>();
        ArrayList<String> phenotypes = new ArrayList<>();

        
        int count=0;
        
         BufferedReader reader;
         try {
      Path path =Paths.get(cGD.getAbsolutePath());
      System.out.println("Path: "+cGD.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null){
          
         String tmp[] = line.split("\t");
         
         if(count==0){
             count++;
             System.out.println("First: "+tmp[0]);
            // for(int i=1;i<tmp.length;i++)
              //   bacteria.add(tmp[i]);
             continue;          
         }
         else{
             count++;
             String p = tmp[0].trim();
            // System.out.println("PehnoReading: "+p);
             String b = tmp[1].trim();
             int def = Integer.parseInt(tmp[4]);
             
             phenoDTemp.add(p);
             bactDTemp.add(b);
             deffDTemp.add(def);
            //growDTemp.put(b, new Pair<String,Integer>(p,def));
             
         }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
         count=0;
         
         try {
      Path path =Paths.get(cScore.getAbsolutePath());
      System.out.println("Path: "+cScore.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null){
          
         String tmp[] = line.split("\t");
         
         if(count==0){
             count++;
             for(int j=0;j<tmp.length;j++){
                 if(j==0)
                     continue;
                 bacteria.add(tmp[j].trim());
             }
         }
         else{
             count++;
            
                tmp = line.split("\t");
                
                String p = tmp[0].trim();
                phenotypes.add(p);
                
                //System.out.println("tmp: "+tmp[0]);
                
                
                for(int i=0;i<tmp.length;i++){

                    if(i==0)
                        condS.put(p, new ArrayList<>()); 
                    else if(!tmp[i].isEmpty()) condS.get(p).add(Double.parseDouble(tmp[i]));
                    else condS.get(p).add(0.0);
                }
         }          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
         HashSet<String> d0 = new HashSet<>();
         for(String s:bactDTemp){
                d0.add(s);
         }
         
         System.out.println("numBact growthDeff: "+d0.size());
         
         HashSet<String> d = new HashSet<>();
         for(String s:phenoDTemp){
                d.add(s);
         }
         
         System.out.println("numPheno growthDeff: "+d.size());
         
         System.out.println("numPhenoCS: "+condS.keySet().size());
          Iterator<String> its = condS.keySet().iterator();
         System.out.println("numBactCS: "+condS.get(its.next()).size());
         
         for(String s:phenotypes)
             if(!d.contains(s))
                 System.out.println("phenotype "+s+" can not be used");
        
         HashSet<String> bactsGenomes = new HashSet<>();
         HashSet<String> bactsSequences = new HashSet<>();
         
          try{
         Collection files = FileUtils.listFiles(inputGenomes, extension, true);

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 //System.out.println("gen: "+input.getName());
                 String bactName = input.getName().replace(".gff", "");
                 String n[] = bactName.split("_");
                 bactName = n[0].trim();
                 
               //  System.out.println("bactName: "+bactName);
                 bactsGenomes.add(bactName);
                 
            }
        
                } catch (Exception e) {
            e.printStackTrace();
        }
          
             try{
         Collection files = FileUtils.listFiles(inputSequences, extension1, true);

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                
              //  System.out.println("seq: "+input.getName());
                 
                 String bactName = input.getName().replace(".fasta", "");
                 String n[] = bactName.split("_");
                 bactName = n[0].trim();
                 
                 System.out.println("bactName: "+bactName);
                 bactsSequences.add(bactName);
                 
            }
        
                } catch (Exception e) {
            e.printStackTrace();
        }
             
             
              HashSet<String> bactsTmpInter = new HashSet<>();
              HashSet<String> bactsInterSet = new HashSet<>();
             
             System.out.println("bacteria size: "+bacteria.size());
             
             for(String s: bacteria)
                 if(d0.contains(s))
                     bactsTmpInter.add(s);
             
             for(String s: bactsTmpInter){
                 if(bactsGenomes.contains(s) && bactsSequences.contains(s))
                     bactsInterSet.add(s);
             }

             
          try{
 
            int numIter = 1;//(int)Math.ceil((double)numSeq/numPerFile);
            

        
            FileWriter fw = new FileWriter(new File(outputPhenotypes)); 
            
            for(String s:phenotypes)
                fw.write(s+"\n");
            
        fw.close();
          
   }catch (Exception e) {
            e.printStackTrace();
        }
          
          
            try{
 
            int numIter = 1;//(int)Math.ceil((double)numSeq/numPerFile);

        
            FileWriter fw = new FileWriter(new File(outputBacteria)); 
            for(String s:bactsInterSet)
                fw.write(s+"\n");
                       
        fw.close();
          
   }catch (Exception e) {
            e.printStackTrace();
        }
            
            
             try{
 
            int numIter = 1;//(int)Math.ceil((double)numSeq/numPerFile);

        
            FileWriter fw = new FileWriter(new File(output)); 
            
            fw.write("@RELATION   phenData\n\n");
            
            fw.write("@ATTRIBUTE   ID"+"    "+"string"+"\n");
            
            for(int i=0;i<phenotypes.size();i++){
                fw.write("@ATTRIBUTE   "+phenotypes.get(i)+"    "+"numeric"+"\n");
            }
            
            
             for(int i=0;i<phenotypes.size();i++){
                 if(i+1<phenotypes.size())
                      fw.write("@ATTRIBUTE   pheno"+(i+1)+"   {1,0}\n");
                 else fw.write("@ATTRIBUTE   pheno"+(i+1)+"   {1,0}\n\n");
             }
             
            fw.write("@DATA\n");
            
            
            
            for(int i=0;i<bacteria.size();i++){
                String line = bacteria.get(i)+", ";
                if(!bactsInterSet.contains(bacteria.get(i)))
                    continue;
                for(int j=0;j<phenotypes.size();j++){
                    ArrayList<Double> tmp = condS.get(phenotypes.get(j));
                    if(tmp.size()<695){
                        System.out.println("Wrong size: "+j);
                        System.out.println("tmp size: "+tmp.size());
                        System.out.println(phenotypes.get(j));
                    }
                        
                    line+=tmp.get(i)+", ";
                   // fw.write(tmp.get(i)+", ");
                }
                
               for(int j=0;j<phenotypes.size();j++){ 
                   int found=0;
                for(int k=0;k<phenoDTemp.size();k++){
                    if(phenoDTemp.get(k).equals(phenotypes.get(j)) && bactDTemp.get(k).equals(bacteria.get(i))){
                        found=1;
                        if(j+1<phenotypes.size())
                            line+=deffDTemp.get(k)+", ";
                            //fw.write(deffDTemp.get(k)+", ");
                        else
                            line+=deffDTemp.get(k)+"\n";
                            //fw.write(deffDTemp.get(k)+"\n");
                        break;
                    }       
                }
                if(found==0){
                    if(j+1<phenotypes.size())
                            line+="0, ";
                            //fw.write(deffDTemp.get(k)+", ");
                        else
                            line+="0\n";
                   
                }
               }
               fw.write(line);
            }
                       
        fw.close();
          
   }catch (Exception e) {
            e.printStackTrace();
        }
         
        
    }
    
}
