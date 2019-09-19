/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create geneOG mappings
 */
 
public class GenesAndOgsOrganism {
    final static Charset ENCODING = StandardCharsets.UTF_8;
    public HashSet<String> organismOgs;
    public HashMap<String,HashSet<Pair<String,String>>> geneOGMappingOrganism;
    public HashMap<String,HashSet<String>> OGGeneMappingOrganism;
    
    public GenesAndOgsOrganism(){
        organismOgs = new HashSet<>();
        geneOGMappingOrganism = new HashMap<>();
        OGGeneMappingOrganism = new HashMap<>();
    }
    
    public GenesAndOgsOrganism(GenesAndOgsOrganism c){
        organismOgs = new HashSet<>();
        geneOGMappingOrganism = new HashMap<>();
        OGGeneMappingOrganism = new HashMap<>();
        
        organismOgs.addAll(c.organismOgs);
        
        Iterator<String> it = c.geneOGMappingOrganism.keySet().iterator();
        
        while(it.hasNext()){
            String k = it.next();
            HashSet<Pair<String,String>> th = c.geneOGMappingOrganism.get(k);
            HashSet<Pair<String,String>> thNew = new HashSet<>();
            
            for(Pair p: th){
                Pair<String,String> pn = new Pair(p.getValue0(),p.getValue1());
                thNew.add(pn);
            }
            
            geneOGMappingOrganism.put(k, thNew);          
        }
        
        it = c.OGGeneMappingOrganism.keySet().iterator();
        
        while(it.hasNext()){
            String k = it.next();
            HashSet<String> th = c.OGGeneMappingOrganism.get(k);
            HashSet<String> thNew = new HashSet<>();
            
            for(String s: th){
                thNew.add(s);
                
                OGGeneMappingOrganism.put(k, thNew);
            }
        }       
    }
    
    public void traverseGenomes(File folderPath, int taxIdT, String[] extensions, Boolean recursive, String taxIDFilePath, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, Mappings map, CreateReducedMappingsFile rmf, int mode){
        
        File ensembleTIDsFile = new File("ensembleTIDs.txt");
          BufferedReader reader=null;
       File taxIDFile=new File(taxIDFilePath);
       HashMap<Integer,Integer> taxTranslation= new HashMap<>();
       
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
                            taxTranslation.put(strainTax, speciesTax);
                            }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();} 
       
       
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
        
        try{
         
         Collection files = FileUtils.listFiles(folderPath, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         
          HashSet<Integer> ensembleTaxs = new HashSet();
          
           BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
          
              HashSet<Integer> eggnogTaxs = new HashSet();
         
         Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
               
                for (Iterator iterator = files.iterator(); iterator.hasNext();){ 
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0; 
                
                    System.out.println("Translated tax of"+taxid+" +: "+taxTranslation.get(taxid));
                    System.out.println("Translated tax of"+taxIdT+" +: "+taxTranslation.get(taxIdT));
                
                    if(taxid!=taxIdT)
                        found=0;
                    else found=1;
              
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;

                if(!input.getAbsolutePath().contains(".nonchromosomal.")){
                    if(mode==0)
                         getGenesAndOgs(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);
                    else
                         getGenesAndOgsMetazoa(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);//metazoa
                }
                else{
                     //make a nonchromosomal version of getGenesAndOgs...
                    if(mode==0)
                        getGenesAndOgsNonChromosomal(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);
                    else
                        getGenesAndOgsNonChromosomalMetazoa(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);//metazoa
                }
                    
                }      
           }
           catch(Exception e){
               e.printStackTrace();
           }
    }
    
    
    
     public void traverseGenomesAndCount(File folderPath, int taxIdT, String[] extensions, Boolean recursive, String taxIDFilePath, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, Mappings map, CreateReducedMappingsFile rmf, HashMap<Integer,Integer> geneCount, HashMap<Integer,Integer> ogCount){
              
        File ensembleTIDsFile = new File("ensembleTIDs.txt");
          BufferedReader reader=null;
       File taxIDFile=new File(taxIDFilePath);
       HashMap<Integer,Integer> taxTranslation= new HashMap<>();
       
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
                            taxTranslation.put(strainTax, speciesTax);
                            }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();} 
       
       
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
        
        try{
         
         Collection files = FileUtils.listFiles(folderPath, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         
          HashSet<Integer> ensembleTaxs = new HashSet();
          
           BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
          
              HashSet<Integer> eggnogTaxs = new HashSet();
         
         Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
               
                for (Iterator iterator = files.iterator(); iterator.hasNext();){ 
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0; 
                
                    System.out.println("Translated tax of"+taxid+" +: "+taxTranslation.get(taxid));
                    System.out.println("Translated tax of"+taxIdT+" +: "+taxTranslation.get(taxIdT));
                
                    if(taxid!=taxIdT)
                        continue;
                  
                System.out.println("taxSpecies jednaki!");
                     
                if(!input.getAbsolutePath().contains(".nonchromosomal."))
                      getGenesAndOgsMetazoa(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);
                else{
                     //make a nonchromosomal version of getGenesAndOgs...
                    getGenesAndOgsNonChromosomalMetazoa(input.getAbsolutePath(),geneOGMap, map,taxid,taxTranslation);
                }
                    
                }      
           }
           catch(Exception e){
               e.printStackTrace();
           }
        
        int c=0;
        
        for(String s: geneOGMappingOrganism.keySet()){
            HashSet<Pair<String,String>> st = geneOGMappingOrganism.get(s);
            if(st.size()>0)
                c++;
        }
        
        int c1=0;
        
        for(String s:OGGeneMappingOrganism.keySet()){
            HashSet<String> st = OGGeneMappingOrganism.get(s);
            if(st.size()>0)
                c1++;
        }
        
        geneCount.put(taxIdT, c);
        ogCount.put(taxIdT, c1);      
    }
     
     public void clearMappins(){
         OGGeneMappingOrganism.clear();
         geneOGMappingOrganism.clear();
     }
    
    public void getGenesAndOgs(String input, HashMap<String,HashSet<Pair<String,String>>> geneOGmap, Mappings map, int taxID,  HashMap<Integer,Integer> taxTranslation){//for fungy

        int lineNum=0, minCoordinate, maxCoordinate;
        HashSet<String> genes = new HashSet<>();

        File inF = new File(input);
        
        BufferedReader reader;
         try {
      Path path =Paths.get(inF.getAbsolutePath());
      System.out.println("Path: "+inF.getAbsolutePath());
      reader = new BufferedReader(new InputStreamReader(new FileInputStream(inF.getAbsolutePath()),"utf-8"));
      String line = null;
      int geneSection = 0;
      int x=-1,y=-1;
      while ((line = reader.readLine()) != null) {
          lineNum++;
                    
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
                              
                                if(gen.contains(".") && taxID!=284812){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                                }
                               
         
         genes.add(gen);
         
          HashSet<Pair<String,String>> ogs = new HashSet();
          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }
                                      }
      geneSection=0;
    }
   }
      System.out.println("GeneOG mapping size: "+geneOGmap.size());
      System.out.println("Chromosome New mapping size: "+geneOGMappingOrganism.size());//iz nekog razloga prazno...
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    public void getGenesAndOgsNonChromosomal(String input, HashMap<String,HashSet<Pair<String,String>>> geneOGmap, Mappings map,int taxID,  HashMap<Integer,Integer> taxTranslation){

        int lineNum=0, minCoordinate, maxCoordinate;
        HashSet<String> genes = new HashSet<>();

        File inF = new File(input);
        
        BufferedReader reader;
         try {
      Path path =Paths.get(inF.getAbsolutePath());
      System.out.println("Path: "+inF.getAbsolutePath());
      reader = new BufferedReader(new InputStreamReader(new FileInputStream(inF.getAbsolutePath()),"utf-8"));
      String line = null;
      int geneSection = 0;
      int x=-1,y=-1;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          
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
                              
                                if(gen.contains(".") && taxID!=284812){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                                }
                                     
         genes.add(gen);
 
         HashSet<Pair<String,String>> ogs = new HashSet();
          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }
          }
    }
   }
      System.out.println("Non-Chromosome New mapping size: "+geneOGMappingOrganism.size());
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
     public void getGenesAndOgsMetazoa(String input, HashMap<String,HashSet<Pair<String,String>>> geneOGmap, Mappings map, int taxID,  HashMap<Integer,Integer> taxTranslation){//for fungy

        int lineNum=0, minCoordinate, maxCoordinate;
        HashSet<String> genes = new HashSet<>();

        File inF = new File(input);
        
        BufferedReader reader;
         try {
      Path path =Paths.get(inF.getAbsolutePath());
      System.out.println("Path: "+inF.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int geneSection = 0, proteinSection=0;
      int x=-1,y=-1;
      while ((line = reader.readLine()) != null) {
          lineNum++;
                    
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
                          
                          if(line.contains("CDS") && taxID != 7165 && taxID!=7460 && taxID!=7425 && taxID!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxID!=6183)
                                  gen = gen.split("\\.")[0];

                              geneSection=0;
                              
                                if(gen.contains(".")){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                               }
                                
                                if(!(taxID==7165 || taxID == 7460 || taxID == 7425 || taxID == 6183)){
                                    continue;
                                }

                               
        if(x!=-1 && y!=-1) 
         genes.add(gen);
        
          HashSet<Pair<String,String>> ogs = new HashSet();
          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }
                                      }
      geneSection=0;
    }
                          
                 if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxID!=6183 && taxID!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                               
                               if(x!=-1 && y!=-1)
                                genes.add(gen);
                               
                                HashSet<Pair<String,String>> ogs = new HashSet();
          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }

                                 proteinSection=0;
                              }
                 }
                          
                          
   }
      System.out.println("GeneOG mapping size: "+geneOGmap.size());
      System.out.println("Chromosome New mapping size: "+geneOGMappingOrganism.size());//iz nekog razloga prazno...
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    public void getGenesAndOgsNonChromosomalMetazoa(String input, HashMap<String,HashSet<Pair<String,String>>> geneOGmap, Mappings map,int taxID,  HashMap<Integer,Integer> taxTranslation){

        int lineNum=0, minCoordinate, maxCoordinate;
        HashSet<String> genes = new HashSet<>();

        File inF = new File(input);
        
        BufferedReader reader;
         try {
      Path path =Paths.get(inF.getAbsolutePath());
      System.out.println("Path: "+inF.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int geneSection = 0, proteinSection=0;
      int x=-1,y=-1;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          
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
                          
                          if(line.contains("CDS") && taxID != 7165 && taxID!=7460 && taxID!=7425 && taxID!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
  
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              
                              if(taxID!=6183)
                                  gen = gen.split("\\.")[0];
                               geneSection=0;
                              
                                if(gen.contains(".") && taxID!=284812){
                                   String tmp[] = gen.split("\\.");
                                   if(tmp.length==2)
                                       gen=tmp[0].trim();
                                   else if(tmp.length==3)
                                       gen=tmp[1].trim();
                                }
                                
                                 if(!(taxID==7165 || taxID == 7460 || taxID == 7425 || taxID == 6183)){
                                    continue;
                                }
       
         if(x!=-1 && y!=-1)
            genes.add(gen);
 
         HashSet<Pair<String,String>> ogs = new HashSet();
          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }   
      }
    }
                          
                      if(line.contains("/protein_id") && proteinSection==1){
                                        String  gen = line.split("=")[1].trim();
                                         gen = gen.replaceAll("\"", "");
                                         gen = gen.replace("-PA", "");
                                         gen = gen.replace("-tr", "");
                                         if(taxID!=6183 && taxID!=6239)
                                                  gen = gen.split("\\.")[0];
                                                proteinSection=0;
                               
                               if(x!=-1 && y!=-1)
                                genes.add(gen);
                               
                                HashSet<Pair<String,String>> ogs = new HashSet();

          if(geneOGmap.containsKey(gen)){
             HashSet<Pair<String,String>> OGsT=new HashSet<>();
              OGsT=geneOGmap.get(gen);
               HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:OGsT){
                if(Integer.parseInt(p.getValue1())==taxID || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxID))){
                    OGs.add(p.getValue0());
                    ogs.add(p);
                }
            }
              
                organismOgs.addAll(OGs);
        
         if(!geneOGMappingOrganism.containsKey(gen))
             geneOGMappingOrganism.put(gen, ogs);
         else{
             geneOGMappingOrganism.get(gen).addAll(ogs);
         }
         
         for(String og:OGs){
             if(!OGGeneMappingOrganism.containsKey(og)){
                 OGGeneMappingOrganism.put(og, new HashSet<String>());
                 OGGeneMappingOrganism.get(og).add(gen);
             }
             else{
                 OGGeneMappingOrganism.get(og).add(gen);
             }
         }
                                 proteinSection=0;
                              }
                 }                                            
   }
      System.out.println("Non-Chromosome New mapping size: "+geneOGMappingOrganism.size());
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
}
