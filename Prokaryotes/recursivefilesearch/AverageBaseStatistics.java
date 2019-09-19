/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import org.javatuples.Pair;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create comparative model performance file (contains comparative AUC/AUPRC scores)
 */
public class AverageBaseStatistics {
     public static void main(String[] args) {
        HashMap<String,Pair<Double,Double>> average=new HashMap<>();
        HashMap<String,Pair<Double,Double>> baseline=new HashMap<>();
       // File functionFile= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\BaselineOGsk=4_600_.oob");//OK prokaryot
        //File baselineFile = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genome5LogK=10.out");//OK prokaryot
       // File functionFile= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\BaselineOGsFungyBB_600_.oob");//OK metazoa
       // File baselineFile = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genomBBLogDistK=10.out");//OK metazoa
        File functionFile= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\BaselineOGsMetazoaBB_600_.oob");
        File baselineFile = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genomeMetazoaBB1.out");
        //File functionFile= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\BacterialSettingDifOnt_600_J=1.0.oob"); //different ontology setting
        //File baselineFile = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\BacterialSettingDifOnt_600_BPOnly.oob");
        HashMap<String,Double> realFreq=new HashMap<>();
        int knntest=1;
        int useAUC=1;
        
        //String frequencyFile = "Z:/databases/NCBI_genomes/auxilliary_files_withIDs/Uniprot-freqs-from idmapping-2014-12-01.txt";//OK prokaryot
       // String frequencyFile = "C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\GOFrequencyFungyN.txt";//OK fungi
        String frequencyFile = "C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\GOFrequencyMetazoaN.txt";//OK metazoa
        
         
        try (BufferedReader bufRdr = new BufferedReader(new FileReader(frequencyFile)))
        {
            String line;
            
            while ((line = bufRdr.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");

                int go = Integer.parseInt(parts[0]);
                String gString=go+"";
                while(gString.length()<7){
                    gString="0"+gString;
                }
                gString="GO"+gString;
                double generality = Double.parseDouble(parts[1]);
                realFreq.put(gString, generality);
                System.out.println("GO: "+gString);
            }
        }
        catch(Exception e){}

        BufferedReader reader;
         try {
      Path path =Paths.get(functionFile.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              String AUC = tmp[1];
              AUPRC=tmp[2].split(": ")[1];
               AUC=tmp[1].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              if(useAUC==0)
                 average.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
              else
                 average.put(first, new Pair(Double.parseDouble(AUC),realFreq.get(first))); 
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
      if(knntest==0){   
             try {
      Path path =Paths.get(baselineFile.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }
       
          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              baseline.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)));
          }        
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            } 
      }
         
       if(knntest==1){
         try {
      Path path =Paths.get(baselineFile.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int summary=1, dataSection=0, test=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Summary")){
              summary=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(summary==1){
              if(line.contains("Testing error")){
                  test=1;
                  summary=0;
                  continue;
              }
          }
          
          if(test==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  test=0;
                  continue;
              }
          }
        
          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("firstBase: "+first);
              String AUPRC=tmp[2];
              String AUC = tmp[1];
              AUPRC=tmp[2].split(": ")[1];
              AUC=tmp[1].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
               if(useAUC==0)
                baseline.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
               else
                  baseline.put(first, new Pair(Double.parseDouble(AUC),realFreq.get(first))); 
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            } 
       }
       
        // String out="C:\\Users\\matej\\Desktop\\ComparisonUniprot.txt";//OK prokaryot
         // String out="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\ComparisonUniprotFungi.txt";//OK fungi
           String out="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\ComparisonUniprotMetazoa.txt";//OK metazoa
          
          
         try{
         FileWriter fw = new FileWriter(out);
         
         for(String s:average.keySet()){
             if(!baseline.containsKey(s))
                 System.out.println(s+" not contained in baseline");
             fw.write(baseline.get(s).getValue0()+" "+average.get(s).getValue0()+" "+average.get(s).getValue1()+"\n");
         }
         
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}         
     }
}
