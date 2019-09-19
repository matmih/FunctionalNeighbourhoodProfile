/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import javafx.util.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to load bacterial gene data.
 */
public class LoadGeneData {
    final static Charset ENCODING = StandardCharsets.UTF_8;
    public HashMap<String,ArrayList<Integer>> contigs;
    public HashMap<String,Pair<ArrayList<Integer>,String>> genes; //array contains minInd, maxInd, strain
   
    public LoadGeneData(){
        contigs = new HashMap<>();
        genes = new HashMap<>();
    }
    
    public void loadData(File genome){
        
        
        BufferedReader reader;
         try {
      Path path =Paths.get(genome.getAbsolutePath());
      System.out.println("Path: "+genome.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("##gff-version"))
              continue;
        
          if(line.contains("##sequence-region")){
              String tmp[] = line.split(" ");
              String contid = tmp[1].trim();
              int min = Integer.parseInt(tmp[2].trim());
              int max = Integer.parseInt(tmp[3].trim());
              ArrayList<Integer> tmpL = new ArrayList<>();
              tmpL.add(min); tmpL.add(max);
              contigs.put(contid, tmpL);
          }
          else{
              String tmp[] = line.split("\t");
              if(!(tmp[2].trim()).equals("gene"))
                  continue;
              
              int min = Integer.parseInt(tmp[3].trim());
              int max =  Integer.parseInt(tmp[4].trim());
              String straind = tmp[6].trim();
              int strInt;
              
              if(straind.equals("-"))
                  strInt =0;
              else strInt = 1;
              
              String desc = tmp[8];
              String descA[] = desc.split(";");
              String geneID="";
              
              for(int i=0;i<descA.length;i++)
                  if(descA[i].contains("ID") && !descA[i].contains("GeneID")){
                      String t[] = descA[i].trim().split("=");
                      geneID = t[1].trim();
                  }
                  else if(descA[i].contains("GeneID")){
                      String t[] = descA[i].trim().split(",");
                      String t1[] = t[1].trim().split(":");
                      geneID = t1[1].trim();
                  }
                  else if(!descA[i].contains("ID") && !descA[i].contains("GeneID")){
                      if(descA[i].contains("locus_tag")){
                      String t[] = descA[i].trim().split("=");
                      geneID = t[1].trim()+"_gene";
                      }
                  } 
              
              if(geneID.equals("")){
                  System.out.println("Warning.............!!!!");
                  System.out.println(genome.getAbsolutePath());
                  System.out.println("Empty gene name.....!!!!");
              }
              
              String contigID = tmp[0].trim();
              
              ArrayList<Integer> tmpArr = new ArrayList<Integer>();
              tmpArr.add(min); tmpArr.add(max); tmpArr.add(strInt);
              
              Pair<ArrayList<Integer>,String> p = new Pair<>(tmpArr,contigID);
              genes.put(geneID, p);
          }      
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }      
        
    }
    
    public void shortOutput(){
           System.out.println("Num Contigs: "+contigs.keySet().size());
           
           Iterator<String> it = genes.keySet().iterator();
           int count=0;
           
           while(it.hasNext()){
               String gene = it.next();
               Pair<ArrayList<Integer>,String> p = genes.get(gene);
               System.out.println(gene+" "+p.getKey().get(0)+" "+p.getKey().get(1)+" "+p.getValue());
               if(count>8)
                   return;
               count++;
           }
    }
    
}
