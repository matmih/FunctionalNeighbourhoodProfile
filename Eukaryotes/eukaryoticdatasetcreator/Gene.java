/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import org.apache.commons.lang3.math.NumberUtils;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about genes and their locations
 */
public class Gene {
    final static Charset ENCODING = StandardCharsets.UTF_8;
    ArrayList<Triplet<Integer,Integer,String>> genes=new ArrayList<Triplet<Integer,Integer,String>>();//all genes
    int minCoordinate=0,maxCoordinate=0;//we save the maximum coordinate and take maxCoordinate=0 to compute min distance in both direction

    
    ArrayList<ArrayList<Triplet<Integer,Integer,String>>> genesContig=new ArrayList<ArrayList<Triplet<Integer,Integer,String>>>();
    
    HashMap<String,Integer> cogToIndex=new HashMap<String,Integer>();//contains COG name COG index mapping
    HashSet<String> anotGen=new HashSet<>();
    
    void findCogs(File input,  HashMap<String,HashSet<String>> geneOgs, Mappings map){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 geneSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    
    void findCogsNew(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null){//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                              
                              
             if(gen.equals("FOXG_05859") || (gen.equals("FVEG_03736T0")) || (gen.equals("FVEG_03736")) || (gen.equals("FOXG_05859P0")) || (gen.equals("NechaG86793")) || (gen.equals("NechaG82508")) || (gen.equals("NechaP86793")) || (gen.equals("NechaP82508"))){
                 System.out.println(gen+" gg");
            }
                                 geneSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    
     void findCogsNewNoDuplicates(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null){//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen)){
                                    int contained  =0;
                                    
                                   if(map.mappings.containsKey(gen)){ //eliminate duplicate proteins
                                     
                         HashSet<Pair<String,String>> cand = map.mappings.get(gen);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gen)){
                            if(anotGen.contains(p.getValue1()))
                                contained=1;
                            break;
                        }      
                    }
                  }
                                if(contained  == 0)  
                                 anotGen.add(gen);
                              }
                              
                              
             if(gen.equals("FOXG_05859") || (gen.equals("FVEG_03736T0")) || (gen.equals("FVEG_03736")) || (gen.equals("FOXG_05859P0")) || (gen.equals("NechaG86793")) || (gen.equals("NechaG82508")) || (gen.equals("NechaP86793")) || (gen.equals("NechaP82508"))){
                 System.out.println(gen+" gg");
            }
                                 geneSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    
     void findCogsNewNoSyntheni(File input, SynthenicBlocks blocks, HashMap<Integer,Integer> taxIDTranslationMap, int taxid , HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      ArrayList<Pair<Integer,Integer>> p = new ArrayList<>();
      
      int translated = taxIDTranslationMap.get(taxid);
                              
                              Iterator<Integer> it = blocks.blocks.keySet().iterator();
                              
                              while(it.hasNext()){
                                  int ti = it.next();
                                  if(!taxIDTranslationMap.containsKey(ti))
                                      continue;
                                  if(taxIDTranslationMap.get(ti) == translated){
                                      p=blocks.blocks.get(ti);
                                  }
                              }
  
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null){//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                              
                              
                                for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }
                                
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                              
                              
             if(gen.equals("FOXG_05859") || (gen.equals("FVEG_03736T0")) || (gen.equals("FVEG_03736")) || (gen.equals("FOXG_05859P0")) || (gen.equals("NechaG86793")) || (gen.equals("NechaG82508")) || (gen.equals("NechaP86793")) || (gen.equals("NechaP82508"))){
                 System.out.println(gen+" gg");
            }
                                 geneSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    void findCogsNewLeaveOut(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, HashSet<String> leaveoutOGs){ //used to create Cafa2 dataset
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null){//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              int leaveOut = 0;
                              
                              HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              
                              for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0)
                                 anotGen.add(gen);
                              
                              
             if(gen.equals("FOXG_05859") || (gen.equals("FVEG_03736T0")) || (gen.equals("FVEG_03736")) || (gen.equals("FOXG_05859P0")) || (gen.equals("NechaG86793")) || (gen.equals("NechaG82508")) || (gen.equals("NechaP86793")) || (gen.equals("NechaP82508"))){
                 System.out.println(gen+" gg");
            }
                                 geneSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    
    void findCogsNewMetazoaLeaveOut(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, HashSet<String> leaveoutOGs, int taxonId){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                               int leaveOut = 0;
                              
                              HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              
                              for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0)
                                 anotGen.add(gen);

                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;

                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                                 
                                  int leaveOut = 0;
                              
                              HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              
                              for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0)
                                 anotGen.add(gen);
                                 proteinSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }    
    
     void findCogsNewMetazoa(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, int taxonId){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                               geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                          
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 proteinSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
     void findCogsNewMetazoaNoDuplicates(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, int taxonId){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                               if(geneOgs.containsKey(gen)){
                                    int contained  =0;
                                    
                                   if(map.mappings.containsKey(gen)){ //eliminate duplicate proteins
                                     
                         HashSet<Pair<String,String>> cand = map.mappings.get(gen);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gen)){
                            if(anotGen.contains(p.getValue1()))
                                contained=1;
                            break;
                        }      
                    }
                  }
                                if(contained  == 0)  
                                 anotGen.add(gen);
                              }
                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;

                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 proteinSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
     
     
     void findCogsNewMetazoaNoSyntheni(File input, SynthenicBlocks blocks, HashMap<Integer,Integer> taxIDTranslationMap, int taxid, HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, int taxonId){ 
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
       ArrayList<Pair<Integer,Integer>> p = new ArrayList<>();
      
      int translated = taxIDTranslationMap.get(taxid);
                              
                              Iterator<Integer> it = blocks.blocks.keySet().iterator();
                              
                              while(it.hasNext()){
                                  int ti = it.next();
                                  if(!taxIDTranslationMap.containsKey(ti))
                                      continue;
                                  if(taxIDTranslationMap.get(ti) == translated){
                                      p=blocks.blocks.get(ti);
                                  }
                              }
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                    
                                    
                                for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }
                                    
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                          
                                
                               for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen))
                                 anotGen.add(gen);
                                 proteinSection=0;
                              }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
     
     
      int findCogsNewnonCHMetazoaLeaveOut(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, HashSet<String> leaveoutOGs ,int taxonId){ 
        int lineNum=0;
  int maxContigSize=0;
  HashSet<String> anotTmp = new HashSet<>();
  
        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];

                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);

                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              int leaveOut = 0;
                               HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1){
  
                              for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0)
                                 anotGen.add(gen);
                              }
                                                                 
                                 anotTmp.add(gen);
                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                          
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);

                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                                 int leaveOut = 0;
                                 HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1){
                                   for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0)
                                 anotGen.add(gen);
                              }
                                 anotTmp.add(gen);
                                 
                                 proteinSection=0;
                              }
                               
                                if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }                            
                         }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
    
    int findCogsNewnonCHMetazoa(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, int taxonId){ 
        int lineNum=0;
  int maxContigSize=0;
  HashSet<String> anotTmp = new HashSet<>();
  
        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1)
                                 anotGen.add(gen);
                                 anotTmp.add(gen);
                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;

                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1)
                                 anotGen.add(gen);
                                 anotTmp.add(gen);
                                 proteinSection=0;
                              }
                               
                                if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }                            
                         }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
    
     int findCogsNewnonCHMetazoaNoSyntheni(File input,  SynthenicBlocks blocks, HashMap<Integer,Integer> taxIDTranslationMap, int taxid  , HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, int taxonId){ 
        int lineNum=0;
  int maxContigSize=0;
  HashSet<String> anotTmp = new HashSet<>();
  
        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
       ArrayList<Pair<Integer,Integer>> p = new ArrayList<>();
      
      int translated = taxIDTranslationMap.get(taxid);
                              
                              Iterator<Integer> it = blocks.blocks.keySet().iterator();
                              
                              while(it.hasNext()){
                                  int ti = it.next();
                                  if(!taxIDTranslationMap.containsKey(ti))
                                      continue;
                                  if(taxIDTranslationMap.get(ti) == translated){
                                      p=blocks.blocks.get(ti);
                                  }
                              }
      
      int geneSection=0,proteinSection=0;
      int x=-1,y=-1;
      int lineCount=0;
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               if(NumberUtils.isParsable(coord[0]) && NumberUtils.isParsable(coord[1])){
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               }
                               else{x=-1; y=-1;}
                               continue;
                        }
                           if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183)){
                                    continue;
                                }
                                else{
                                    
                                    
                                for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }
                                    
                                     Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1)
                                 anotGen.add(gen);
                                 anotTmp.add(gen);

                                 geneSection=0;
                              }
                                }
                                
                               if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxonId!=6183 && taxonId!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                                                
                               
                                for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }                 

                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                               if(x!=-1 && y!=-1)
                                genes.add(t);

                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen) && x!=-1 && y!=-1)
                                 anotGen.add(gen);
                                 anotTmp.add(gen);

                                 proteinSection=0;
                              }
                               
                                if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }                            
                         }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
    
     
     int findCogsNewnonCHLeaveOut(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map, HashSet<String> leaveoutOGs){ 
        int lineNum=0;
        int maxContigSize=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
      HashSet<String> anotTmp = new HashSet<>();
              
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                 if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates NCH!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen); 
                                 
                                  int leaveOut = 0;
                              
                              HashSet<Pair<String,String>> ogsgene = geneOgs.get(gen);
                              
                              for(Pair<String,String> p:ogsgene){
                                  if(leaveoutOGs.contains(p.getValue0().trim())){
                                        leaveOut=1;
                                        break;
                                  }
                              }
                              
                              if(geneOgs.containsKey(gen) && leaveOut == 0){
                                 anotGen.add(gen);
                                 anotTmp.add(gen);
                              }
                                 geneSection=0;
                              }
                          
                          if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }                         
                       }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
    
    
     int findCogsNewnonCH(File input,  HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map){ 
        int lineNum=0;
        int maxContigSize=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
      HashSet<String> anotTmp = new HashSet<>();
              
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                 if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates NCH!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen)){
                                 anotGen.add(gen);
                                 anotTmp.add(gen);
                              }
                                 geneSection=0;
                              }
                          
                          if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }                         
                         }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
    
     int findCogsNewnonCHNoSyntheni(File input, SynthenicBlocks blocks, HashMap<Integer,Integer> taxIDTranslationMap, int taxid, HashMap<String,HashSet<Pair<String,String>>> geneOgs, Mappings map){ 
        int lineNum=0;
        int maxContigSize=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
       ArrayList<Pair<Integer,Integer>> p = new ArrayList<>();
      
      int translated = taxIDTranslationMap.get(taxid);
                              
                              Iterator<Integer> it = blocks.blocks.keySet().iterator();
                              
                              while(it.hasNext()){
                                  int ti = it.next();
                                  if(!taxIDTranslationMap.containsKey(ti))
                                      continue;
                                  if(taxIDTranslationMap.get(ti) == translated){
                                      p=blocks.blocks.get(ti);
                                  }
                              }
      
      int geneSection=0;
      int x=-1,y=-1;
      int lineCount=0;
      HashSet<String> anotTmp = new HashSet<>();
              
         while ((line = reader.readLine()) != null) {//parse a .dat file
            
                            if(line.contains("source") && line.contains("..") && !line.contains("[") && !line.contains("]")){
                                String tmp[]=line.split("source");
                                tmp[1]=tmp[1].trim();
                                String coC[]=tmp[1].split("\\.\\.");
                                minCoordinate=Integer.parseInt(coC[0]);
                                maxCoordinate=Integer.parseInt(coC[1]);
                                 if(Double.isNaN(minCoordinate) || Double.isNaN(maxCoordinate) || minCoordinate>=maxCoordinate){
                                    System.out.println("Wrong coordinates NCH!");
                                    System.out.println(minCoordinate+" "+maxCoordinate);
                                }
                            }
             
                          if(line.contains("db_xref=\"taxon:")){
                              continue;
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              String tmp[]=line.split("gene");
                              tmp[1]=tmp[1].trim();
                              if(tmp[1].contains("complement")){
                                  tmp[1]=tmp[1].replace("complement","");
                                  tmp[1]=tmp[1].replace("(", "");
                                  tmp[1]=tmp[1].replace(")", "");
                              }
                              String coord[]=tmp[1].split("\\.\\.");
                               x=Integer.parseInt(coord[0]); y=Integer.parseInt(coord[1]);
                               continue;
                        }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                               
                                 for(int i=0;i<p.size();i++){
                                    if(x>=p.get(i).getValue0() && x<= p.get(i).getValue1()){
                                        gen+="SynthenicRemove";
                                        break;
                                    }
                                }
                              
                               Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,gen);
                               
                                genes.add(t);
                                 HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                              
                              if(geneOgs.containsKey(gen)){
                                 anotGen.add(gen);
                                 anotTmp.add(gen);
                              }

                                 geneSection=0;
                              }
                          
                          if(line.trim().equals("//")){
                              ArrayList<Triplet<Integer,Integer,String>> ttA = new ArrayList<>();
                             
                              for(int c = 0; c<genes.size();c++)
                                  ttA.add(genes.get(c));
                              if(genes.size()>maxContigSize)
                                  maxContigSize=genes.size();
                              genesContig.add(ttA);
                              genes = new ArrayList<>();
                              
                              if(anotTmp.size()>maxContigSize)
                                  maxContigSize = anotTmp.size();
                              anotTmp.clear();
                          }
                          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         return maxContigSize;
    }
     
     
    void randomize(Random randGen){
       ArrayList<Triplet<Integer,Integer,String>> ShuffledCogs=new ArrayList<Triplet<Integer,Integer,String>>(); 
       HashSet<Integer> usedIndex=new HashSet<>();
       
       for(int i=0;i<genes.size();i++){
           int ind=randGen.nextInt(genes.size()-usedIndex.size());
           int count=0;
           
           for(int j=0;j<genes.size();j++)
               if(!usedIndex.contains(j)){
                   if(count==ind){
                       usedIndex.add(j);
                       ShuffledCogs.add(genes.get(j));
                       break;
                   }
                   count++;
               }           
       }  
       genes=ShuffledCogs;
    }
    
    void randomizeLocation(Random randGen, int numShufflingSteps){
       ArrayList<Triplet<Integer,Integer,String>> ShuffledCogs=new ArrayList<Triplet<Integer,Integer,String>>(); 
       HashSet<Integer> usedIndex=new HashSet<>();
       
    for(int k=0;k<numShufflingSteps;k++){
       for(int i=0;i<genes.size();i++){
           int ind=randGen.nextInt(genes.size()-i);
           int count=0;
       
                       int a, b;
                       String c;
                       
                       a=genes.get(ind+i).getValue0();
                       b=genes.get(ind+i).getValue1();
                       c=genes.get(ind+i).getValue2();
                       
                       Triplet tmp=genes.get(ind+i);
                       tmp=tmp.setAt0(genes.get(i).getValue0());
                       tmp=tmp.setAt1(genes.get(i).getValue1());
                       
                       Triplet tmp1=genes.get(i);
                       tmp1=tmp1.setAt0(a);
                       tmp1=tmp1.setAt1(b);
                       
                       genes.set(ind+i,tmp);
                       genes.set(i, tmp1);
                       
                       ShuffledCogs.add(tmp);
                   count++;
       }  
     }
    }
}
