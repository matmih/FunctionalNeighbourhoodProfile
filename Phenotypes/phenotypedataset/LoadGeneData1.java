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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import javafx.util.Pair;
import org.javatuples.Triplet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to load bacterial gene data.
 */
public class LoadGeneData1 {
     final static Charset ENCODING = StandardCharsets.UTF_8;
    public ArrayList<ArrayList<Triplet<Integer,Integer,String>>> geneContigs; 
    public HashMap<String,Integer> contigs;
    
   
    public LoadGeneData1(){
      geneContigs = new ArrayList<>();
      contigs = new HashMap<>();
    }
    
    public void loadData(File genome){
        
        int count=0;
        
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
               continue;

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
              
              Triplet<Integer,Integer,String> t = new Triplet(min,max,geneID);
              
              ArrayList<Integer> tmpArr = new ArrayList<Integer>();
              tmpArr.add(min); tmpArr.add(max); tmpArr.add(strInt);
              
              Pair<ArrayList<Integer>,String> p = new Pair<>(tmpArr,contigID);
              
              if(!contigs.keySet().contains(contigID)){
                  geneContigs.add(new ArrayList<>());
                  contigs.put(contigID,count++);
                  geneContigs.get(contigs.get(contigID)).add(t);
              }
              else{
                  geneContigs.get(contigs.get(contigID)).add(t);
              }             
             // genes.put(geneID, p);
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
           
          /* Iterator<String> it = genes.keySet().iterator();
           int count=0;
           
           while(it.hasNext()){
               String gene = it.next();
               Pair<ArrayList<Integer>,String> p = genes.get(gene);
               System.out.println(gene+" "+p.getKey().get(0)+" "+p.getKey().get(1)+" "+p.getValue());
               if(count>8)
                   return;
               count++;
           }*/
    }
}
