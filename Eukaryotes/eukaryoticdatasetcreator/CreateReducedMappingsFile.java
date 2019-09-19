/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.io.InputStreamReader;
import java.util.ArrayList;
import org.apache.commons.io.FileUtils;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create reduced mapping files
 */
public class CreateReducedMappingsFile {
  final static Charset ENCODING = StandardCharsets.UTF_8;
        HashSet<String> proteinsEggnog=new HashSet<>();
        HashSet<String> assignedPIDsToUID=new HashSet<>();
        HashSet<String> containedPIDsInMap=new HashSet<>();
        HashSet<String> uniprot_acs=new HashSet<>();
        HashMap<String, HashSet<Pair<String,String>>> NOGProtein = new HashMap<>();
        HashMap<String, HashSet<Pair<String,String>>> ProteinNOG = new HashMap<>();

        public void loadEggnogProteins(File eggnog){
            
            BufferedReader reader;
            
             try {
                     Path path =Paths.get(eggnog.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                            System.out.println("protIDds size: "+protIds.length);
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                if(pid1.length<4)
                                proteinsEggnog.add(pid1[1].trim()); 
                                else if(pid1[1].contains("SP"))
                                    proteinsEggnog.add(pid1[1].trim()+"."+pid1[2].trim());
                            }
                    }
         System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        }
        
        
        public void loadEggnogProteinsMetazoa(File eggnog){
            
            BufferedReader reader;
            
             try {
                     Path path =Paths.get(eggnog.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                            System.out.println("protIDds size: "+protIds.length);
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                    
                                int tid = Integer.parseInt(pid1[0]);
                                if( tid == 7165 || tid == 7425 || tid == 7460){
                                    pid1[1] = pid1[1].replace("-PA", "");
                                }
                                else if(tid == 6183){
                                    pid1[1] = pid1[1].replace("__mRNA", "");
                                }
                                else if(tid == 6239){
                                    pid1[1] = pid1[1]+"."+pid1[2];
                                }
                                
                                if(pid1.length<4)
                                proteinsEggnog.add(pid1[1].trim()); 

                            }
                    }
         System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        }
        
       public void countMappingsEnsemble(File mapping, File geneFile){
           HashSet<String> containedAll=new HashSet<>();
           HashSet<String> containedAllUniprot=new HashSet<>();
            
            EnsembleGenes eGene=new EnsembleGenes();
            eGene.loadGenes(geneFile);
            
            BufferedReader reader;
            
             try{
                     Path path =Paths.get(mapping.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");

                           /* String id1=tmp[1].trim();
                            String id2=tmp[2].trim();*/
                            
                              String id1=tmp[1].trim();
                            if(id1.contains(".")){
                                String tmpt[]=id1.split("\\.");
                                if(tmpt.length<2)
                                    continue;
                                id1=tmpt[0].trim();
                            }
                                
                            String id2=tmp[2].trim();
                             if(id2.contains(".")){
                                String tmpt[]=id2.split("\\.");
                                if(tmpt.length<2)
                                    continue;
                                id2=tmpt[0].trim();
                            }

                            if(eGene.genes.contains(id1) || eGene.genes.contains(id2)){
                                
                               // HashSet<String> containedPIDsInMap=new HashSet<>();
                                if(eGene.genes.contains(id1)){
                                    containedPIDsInMap.add(id1);
                                    if(proteinsEggnog.contains(id1))
                                        containedAll.add(id1);
                                }
                                else if(eGene.genes.contains(id2)){ containedPIDsInMap.add(id2);
                                       if(proteinsEggnog.contains(id2))
                                        containedAll.add(id2);
                                }
                                
                                // System.out.println(id1+" "+id2);

                            String mapType=tmp[3];
                            String types[] = mapType.split("_");
                           
                            if(types[0].toLowerCase().contains("blast") && types[1].toLowerCase().contains("uniprot")){
                              //  System.out.println("conditions met1");
                                String uni=tmp[2];
                                if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(eGene.genes.contains(id1)){
                                        assignedPIDsToUID.add(id1);
                                         if(proteinsEggnog.contains(id1))
                                         containedAllUniprot.add(id1);
                                }
                                else if(eGene.genes.contains(id2)){
                                    assignedPIDsToUID.add(id2);
                                     if(proteinsEggnog.contains(id2))
                                    containedAllUniprot.add(id2);
                                }
                            }
                            else if(types.length==2){
                                    if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast"))){
                                   
                              // System.out.println("conditions met2");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(eGene.genes.contains(id1)){
                                        assignedPIDsToUID.add(id1);
                                         if(proteinsEggnog.contains(id1))
                                        containedAllUniprot.add(id1);
                                }
                                else if(eGene.genes.contains(id2)){
                                    assignedPIDsToUID.add(id2);
                                     if(proteinsEggnog.contains(id2))
                                    containedAllUniprot.add(id2);
                                }
                            } 
                        }
                            else if(types.length==3){
                                 if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast") || types[2].toLowerCase().contains("blast"))){
                               //      System.out.println("conditions met3");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(eGene.genes.contains(id1)){
                                        assignedPIDsToUID.add(id1);
                                         if(proteinsEggnog.contains(id1))
                                        containedAllUniprot.add(id1);
                                }
                                else if(eGene.genes.contains(id2)){
                                    assignedPIDsToUID.add(id2);
                                     if(proteinsEggnog.contains(id2))
                                    containedAllUniprot.add(id2);
                                    }
                                 }
                            }
                      }
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();}    
                
                System.out.println("Fraction of ensemble pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/eGene.genes.size()));
                 System.out.println("Fraction of ensemble pids contained in mapping: "+((double)containedPIDsInMap.size()/eGene.genes.size()));
                 
                  //System.out.println("Fraction of eggnog pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/proteinsEggnog.size()));
                // System.out.println("Fraction of eggnog pids contained in mapping: "+((double)containedPIDsInMap.size()/proteinsEggnog.size()));
                 
                 System.out.println("Fraction of eggnog pids assigned to uniprotID and mapped to ensemble: "+((double)containedAllUniprot.size()/proteinsEggnog.size()));
                 System.out.println("Fraction of eggnog pids contained in mapping and mapped to ensemble: "+((double)containedAll.size()/proteinsEggnog.size()));
                 
                 System.out.println("Good uni size: "+containedAllUniprot.size());
                 System.out.println("All map size: "+containedAll.size());
                 
}
        
        public void countMappingsSelected(File mapping, File inputFolder, File mappingFile){
            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            
            Mappings map=new Mappings();
            map.loadMappings(mappingFile);
            
            BufferedReader reader;
            
             try{
                     Path path =Paths.get(mapping.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");

                            String id1=tmp[1].trim();
                            String id2=tmp[2].trim();
                            
                            if(id1.contains(".")){
                                String tmpt[]=id1.split("\\.");
                                if(tmpt.length<2)
                                    continue;
                                id1=tmpt[0].trim();
                            }
                                
                             if(id2.contains(".")){
                                String tmpt[]=id2.split("\\.");
                                if(tmpt.length<2)
                                    continue;
                                id2=tmpt[0].trim();
                            }

                            if(proteinsEggnog.contains(id1) || proteinsEggnog.contains(id2)){
                                
                               // HashSet<String> containedPIDsInMap=new HashSet<>();
                                if(proteinsEggnog.contains(id1))
                                    containedPIDsInMap.add(id1);
                                else containedPIDsInMap.add(id2);
                                
                                // System.out.println(id1+" "+id2);

                            String mapType=tmp[3];
                            String types[] = mapType.split("_");
                           
                            if(types[0].toLowerCase().contains("blast") && types[1].toLowerCase().contains("uniprot")){
                              //  System.out.println("conditions met1");
                                String uni=tmp[2];
                                if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                            }
                            else if(types.length==2){
                                    if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast"))){
                                   
                              // System.out.println("conditions met2");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                            } 
                        }
                            else if(types.length==3){
                                 if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast") || types[2].toLowerCase().contains("blast"))){
                               //      System.out.println("conditions met3");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                                 }
                            }
                      }
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();}    
                
                System.out.println("Fraction of eggnog pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/proteinsEggnog.size()));
                 System.out.println("Fraction of eggnog pids contained in mapping: "+((double)containedPIDsInMap.size()/proteinsEggnog.size()));
            
                 //add code here to traverse the genomes

         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                // System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                         
                          //System.out.println(line);
                          //System.out.println("geneSection: "+geneSection);
                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              /*if(taxonId!=559292){
                                 wrong=1;
                                  break;
                              }*/
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              //System.out.println("gene: "+gen);
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }
                     // System.out.println("num genes: "+genes.size());
                        reader1.close();  
         if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;
             System.out.println("genes size: "+genes.size());
             int numProcG=0;
             for(String s:genes){
                  numProcG++; 
                 HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s))
                     continue;
                 
                /*for(String s1:matches){ 
                 if(containedPIDsInMap.contains(s1))
                     count++;
                 if(proteinsEggnog.contains(s1))
                     count1++;
                 break;
                }*/
                  int f=0;
                  
                  for(Pair<String,String> s1:matches){ 
                 if(containedPIDsInMap.contains(s1.getValue0()))
                     count++;
                 
                  if(proteinsEggnog.contains(s1.getValue0())){
                     count1++;
                     f=1;
                     break;
                     }
                  }
                  
                  if(f==1)
                      continue;
                  
                   /*for(String s1:matches){ 
                 
                  int numProc=0;
                 // System.out.println("PEgg size: "+proteinsEggnog.size());
                 for(String segg:proteinsEggnog){
                     
                     if(!map.mappings.containsKey(segg))
                         continue;
                     
                     HashSet<String> matchesEgg=map.mappings.get(segg);
                     
                         if(matchesEgg.contains(s1)){
                             count1++;
                             f=1;
                             break;
                         }
                     
                     if(f==1)
                         break;  
                    // numProc++;
                     //System.out.println("Num proc PEggs: "+numProc);
                 }
                 if(f==1)
                     break;
                }*/
               
               if((numProcG)%100==0)
                   System.out.println("Proc genes: "+numProcG);
             }
            // System.out.println("taxon, count, count1, geneSize");
            // System.out.println(taxonId+" "+count+" "+count1+" "+genes.size());
             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);
             System.out.println("genes size: "+genes.size());
             int numProcG=0;
              for(String s:genes){
                   numProcG++; 
                  HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                  if(!map.mappings.containsKey(s))
                     continue;
                  
                  
                  int f=0;
                  
                  for(Pair<String,String> s1:matches){ 
                 if(containedPIDsInMap.contains(s1.getValue0()))
                     count++;
                 
                  if(proteinsEggnog.contains(s1.getValue0())){
                     count1++;
                     f=1;
                     break;
                     }
                  }
                  
                  if(f==1)
                      continue;
                 
                  /*for(String s1:matches){ 
                   int numProc=0;
                // System.out.println("PEgg size: "+proteinsEggnog.size());
                 for(String segg:proteinsEggnog){

                     if(!map.mappings.containsKey(segg))
                         continue;
                     
                     HashSet<String> matchesEgg=map.mappings.get(segg);
                     
                         if(matchesEgg.contains(s1)){
                             count1++;
                             f=1;
                             break;
                         }
                     
                     if(f==1)
                         break;  
                     
                    // numProc++;
                    // System.out.println("Num proc PEggs: "+numProc);
                     
                 }
                 if(f==1)
                     break;
                  
               if((numProcG)%100==0)
                   System.out.println("Proc genes: "+numProcG);
                 
                }*/
                  numProcG++; 
             }
               
            // System.out.println("taxon, count, count1, geneSize");
            // System.out.println(taxonId+" "+count+" "+count1+" "+genes.size());
             taxCountMap.put(taxonId, count1);
         }
         
         if(!taxGeneMap.containsKey(taxonId)){
             taxGeneMap.put(taxonId, genes);
         }
         else{
             HashSet<String> genesOld=taxGeneMap.get(taxonId);
             
             for(String s:genes)
                 genesOld.add(s);
             
             taxGeneMap.put(taxonId, genesOld);
         }
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
        
        Iterator<Integer> it=taxCountMap.keySet().iterator();
        
        while(it.hasNext()){
            int tax=it.next();
           // System.out.println("tax: "+tax);
            if(taxCountMap.get(tax)==0){
               // System.out.println("Size: "+taxGeneMap.get(tax).size());
                for(String s:taxGeneMap.get(tax)){
                    if(proteinsEggnog.contains(s))
                        //System.out.println("egg gene found!");
                        
                    proteinsEggnog.remove(s);
                }
            }
            else{ 
                System.out.println("tax: "+tax);
                System.out.println("Num contained genes: "+taxCountMap.get(tax));
                System.out.println("Num genes total: "+taxGeneMap.get(tax).size());
            }
        }
        
      System.out.println("Fraction of eggnog pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/proteinsEggnog.size()));
      System.out.println("Fraction of eggnog pids contained in mapping: "+((double)containedPIDsInMap.size()/proteinsEggnog.size()));  
        
   }
        
        
        public HashSet<String> loadOKNOgs(File NOGS){
             BufferedReader reader;
            HashSet<String> res=new HashSet<>();
             
             try {
                     Path path =Paths.get(NOGS.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String NOG = tmp[0].trim();
                            NOG=NOG.replace(":", "");
                            
                            res.add(NOG);     
                      }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        return res;
        }
        
        public void loadNogMembersMappings(File NogMembers){
             BufferedReader reader;
            
             try {
                     Path path =Paths.get(NogMembers.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                            String NOG = tmp[1].trim();
                            
                            if(!NOGProtein.containsKey(NOG)){
                                NOGProtein.put(NOG, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> nGenes=NOGProtein.get(NOG);
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                
                                String protein = pid1[1].trim();
                                
                                /*if(pid1.length==4){
                                    protein+="."+pid1[2].trim();
                                    System.out.println("protein: "+protein);
                                }
                                */
                                
                                if(Integer.parseInt(pid1[0].trim()) == 6239){
                                     protein+="."+pid1[2].trim();
                                }
                                
                                //proteinsEggnog.add(pid1[1].trim());
                                if(!ProteinNOG.containsKey(protein))
                                    ProteinNOG.put(protein, new HashSet<Pair<String,String>>());
                                HashSet<Pair<String,String>> nogs=ProteinNOG.get(protein);
                                nGenes.add(new Pair(protein,pid1[0].trim()));
                                nogs.add(new Pair(NOG,pid1[0].trim()));
                                NOGProtein.put(NOG, nGenes);
                                ProteinNOG.put(protein, nogs);
                            }      
                            //System.out.println("Mapping: "+" "+NOG+" "+NOGProtein.get(NOG).size());
                    }
         System.out.println("Mappings loaded");
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        }
        
        
        public void loadNogMembersMappingsMetazoa(File NogMembers){
             BufferedReader reader;
            
             try {
                     Path path =Paths.get(NogMembers.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                            String NOG = tmp[1].trim();
                            
                            if(!NOGProtein.containsKey(NOG)){
                                NOGProtein.put(NOG, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> nGenes=NOGProtein.get(NOG);
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                
                                String protein = pid1[1].trim();
                                
                                  protein = protein.replace("-PA", "");
                                  protein = protein.replace("-PB", "");
                                  protein = protein.replace("-tr", "");
                                
                                if(pid1.length==4){
                                    protein+="."+pid1[2].trim();
                                    System.out.println("protein: "+protein);
                                }
                                
                                if(pid1[0].trim().equals("7460"))
                                    System.out.println("PS: "+protein);
                                
                                //proteinsEggnog.add(pid1[1].trim());
                                if(!ProteinNOG.containsKey(protein))
                                    ProteinNOG.put(protein, new HashSet<Pair<String,String>>());
                                HashSet<Pair<String,String>> nogs=ProteinNOG.get(protein);
                                nGenes.add(new Pair(protein,pid1[0].trim()));
                                nogs.add(new Pair(NOG,pid1[0].trim()));
                                NOGProtein.put(NOG, nGenes);
                                ProteinNOG.put(protein, nogs);
                            }      
                            //System.out.println("Mapping: "+" "+NOG+" "+NOGProtein.get(NOG).size());
                    }
         System.out.println("Mappings loaded");
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        }
        
        public void loadNogMembersMappingsReduced(File NogMembers, HashSet<String> okNogs){
             BufferedReader reader;
            
             try {
                     Path path =Paths.get(NogMembers.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                            String NOG = tmp[1].trim();
                            
                            if(!okNogs.contains(NOG))
                                continue;
                            
                            if(!NOGProtein.containsKey(NOG)){
                                NOGProtein.put(NOG, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> nGenes=NOGProtein.get(NOG);
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                //proteinsEggnog.add(pid1[1].trim());
                                if(!ProteinNOG.containsKey(pid1[1].trim()))
                                    ProteinNOG.put(pid1[1].trim(), new HashSet<Pair<String,String>>());
                                HashSet<Pair<String,String>> nogs=ProteinNOG.get(pid1[1].trim());
                                nGenes.add(new Pair(pid1[1].trim(),pid1[0].trim()));
                                nogs.add(new Pair(NOG,pid1[0].trim()));
                                NOGProtein.put(NOG, nGenes);
                                ProteinNOG.put(pid1[1].trim(), nogs);
                            }      
                            //System.out.println("Mapping: "+" "+NOG+" "+NOGProtein.get(NOG).size());
                    }
         System.out.println("Mappings loaded");
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
        }
        
        public void AssignFunctionToNogs(File geneFunctionFile, double perc ,File outputFile){
            
              BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneFunctions=new HashMap<>();
              HashMap<String, HashSet<String>> NOGFunctions=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneFunctionFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneFunctions.containsKey(gene)){
                                geneFunctions.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> func=geneFunctions.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String go=tmp[i].trim();
                                
                               func.add(go);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneFunctions.put(gene, func);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
            
            Iterator<String> it=NOGProtein.keySet().iterator();
            
            int numOKNOGS=0;
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<Pair<String,String>> genes=NOGProtein.get(nog);
                
                HashSet<String> allGenesWithFunction= new HashSet<String>();
                HashSet<String> allFunctions=new HashSet<>();
                HashSet<String> assignedFunctions=new HashSet<>();
                
                for(Pair<String,String> g:genes){
                    HashSet<String> f=geneFunctions.get(g.getValue0());
                    if(!geneFunctions.containsKey(g.getValue0()))
                        continue;
                    else if(f.size()>0){ allGenesWithFunction.add(g.getValue0());
                        allFunctions.addAll(f);
                    }
                }
                
                //System.out.println(nog+" "+genes.size()+" "+allFunctions.size());
                
                for(String fun:allFunctions){
                    int funCount=0;
                    for(String gen:allGenesWithFunction){
                        HashSet<String> f=geneFunctions.get(gen);
                        if(f.contains(fun))
                            funCount++;
                    }
                    
                    if((funCount/(double)allGenesWithFunction.size())>=perc)
                        assignedFunctions.add(fun);  
                }
                
                //System.out.println(nog+" "+assignedFunctions.size());
                NOGFunctions.put(nog, assignedFunctions); 
                if(assignedFunctions.size()>0)
                    numOKNOGS++;
            } 
            
            System.out.println("Fraction NOGs with function: "+(numOKNOGS/(double)NOGFunctions.keySet().size()));
            System.out.println("Number NOGs with function: "+numOKNOGS);
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 it=NOGFunctions.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<String> f=NOGFunctions.get(nog);
                     
                     fw.write(nog+": ");
                     for(String s:f){
                         fw.write(s+" ");
                     }
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        
        
        
         public void AssignFunctionToNogsNew(File geneFunctionFile, File geneOGFile ,double perc ,File outputFile){
            
              BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneFunctions=new HashMap<>();
              HashMap<String, HashSet<String>> NOGFunctions=new HashMap<>();
              HashMap<String,HashSet<String>> geneOG = new HashMap<>();
              HashMap<String,HashSet<String>> OGGene = new HashMap<>();
              
             try {
                     Path path =Paths.get(geneFunctionFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneFunctions.containsKey(gene)){
                                geneFunctions.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> func=geneFunctions.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String go=tmp[i].trim();
                                
                               func.add(go);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneFunctions.put(gene, func);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
             
             
              try {
                     Path path =Paths.get(geneOGFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length>2)
                                continue;
                           // tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneOG.containsKey(gene)){
                                geneOG.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOG.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].trim();
                                
                                 if(!OGGene.containsKey(og)){
                                OGGene.put(og, new HashSet<String>());
                            }
                             
                                 HashSet<String> tg=OGGene.get(og);
                                 tg.add(gene);
                                 OGGene.put(og, tg);
                                 
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOG.put(gene, ogs);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
            
            Iterator<String> it=OGGene.keySet().iterator();
            
            int numOKNOGS=0;
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<String> genes=OGGene.get(nog);
                 HashSet<Pair<String,String>> genes1=NOGProtein.get(nog);
                
                 for(Pair<String,String> g1:genes1)
                     genes.add(g1.getValue0());
                 
                HashSet<String> allGenesWithFunction= new HashSet<String>();
                HashSet<String> allFunctions=new HashSet<>();
                HashSet<String> assignedFunctions=new HashSet<>();
                
                for(String g:genes){
                    HashSet<String> f=geneFunctions.get(g);
                    if(!geneFunctions.containsKey(g))
                        continue;
                    else if(f.size()>0){ allGenesWithFunction.add(g);
                        allFunctions.addAll(f);
                    }
                }
                
                //System.out.println(nog+" "+genes.size()+" "+allFunctions.size());
                
                for(String fun:allFunctions){
                    int funCount=0;
                    for(String gen:allGenesWithFunction){
                        HashSet<String> f=geneFunctions.get(gen);
                        if(f.contains(fun))
                            funCount++;
                    }
                    System.out.println("funCount: "+funCount);
                    System.out.println("allGenesWithFunction: "+allGenesWithFunction.size());
                    if((funCount/(double)allGenesWithFunction.size())>=perc)
                        assignedFunctions.add(fun);  
                }
                
                //System.out.println(nog+" "+assignedFunctions.size());
                NOGFunctions.put(nog, assignedFunctions); 
                if(assignedFunctions.size()>0)
                    numOKNOGS++;
            } 

            System.out.println("Fraction NOGs with function: "+(numOKNOGS/(double)NOGFunctions.keySet().size()));
            System.out.println("Number NOGs with function: "+numOKNOGS);
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 it=NOGFunctions.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<String> f=NOGFunctions.get(nog);
                     
                     fw.write(nog+": ");
                     for(String s:f){
                         fw.write(s+" ");
                     }
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
         
           public void AssignFunctionToNogsNewNew(File geneFunctionFile, File geneOGFile ,double perc ,File outputFile){
            
              BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneFunctions=new HashMap<>();
              HashMap<String, HashSet<String>> NOGFunctions=new HashMap<>();
              HashMap<String,HashSet<String>> geneOG = new HashMap<>();
              HashMap<String,HashSet<String>> OGGene = new HashMap<>();
              
             try {
                     Path path =Paths.get(geneFunctionFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                           // tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneFunctions.containsKey(gene)){
                                geneFunctions.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> func=geneFunctions.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String tmp1[] = tmp[i].split("::");
                                String go=tmp1[0].trim();//tmp[i].trim();
                                
                               func.add(go);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneFunctions.put(gene, func);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
             
             
              try {
                     Path path =Paths.get(geneOGFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length>2)
                                continue;
                           // tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneOG.containsKey(gene)){
                                geneOG.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOG.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].split(":")[0].trim();
                                
                                 if(!OGGene.containsKey(og)){
                                OGGene.put(og, new HashSet<String>());
                            }
                             
                                 HashSet<String> tg=OGGene.get(og);
                                 tg.add(gene);
                                 OGGene.put(og, tg);
                                 
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOG.put(gene, ogs);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
            
            Iterator<String> it=OGGene.keySet().iterator();
            
            int numOKNOGS=0;
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<String> genes=OGGene.get(nog);
                 HashSet<Pair<String,String>> genes1=NOGProtein.get(nog);
                
                 for(Pair<String,String> g1:genes1)
                     genes.add(g1.getValue0());
                 
                HashSet<String> allGenesWithFunction= new HashSet<String>();
                HashSet<String> allFunctions=new HashSet<>();
                HashSet<String> assignedFunctions=new HashSet<>();
                
                for(String g:genes){
                    HashSet<String> f=geneFunctions.get(g);
                    if(!geneFunctions.containsKey(g))
                        continue;
                    else if(f.size()>0){ allGenesWithFunction.add(g);
                        allFunctions.addAll(f);
                    }
                }
                
                //System.out.println(nog+" "+genes.size()+" "+allFunctions.size());
                
                for(String fun:allFunctions){
                    int funCount=0;
                    for(String gen:allGenesWithFunction){
                        HashSet<String> f=geneFunctions.get(gen);
                        if(f.contains(fun))
                            funCount++;
                    }
                    System.out.println("funCount: "+funCount);
                    System.out.println("allGenesWithFunction: "+allGenesWithFunction.size());
                    if((funCount/(double)allGenesWithFunction.size())>=perc)
                        assignedFunctions.add(fun);  
                }
                
                //System.out.println(nog+" "+assignedFunctions.size());
                NOGFunctions.put(nog, assignedFunctions); 
                if(assignedFunctions.size()>0)
                    numOKNOGS++;
            } 

            System.out.println("Fraction NOGs with function: "+(numOKNOGS/(double)NOGFunctions.keySet().size()));
            System.out.println("Number NOGs with function: "+numOKNOGS);
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 it=NOGFunctions.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<String> f=NOGFunctions.get(nog);
                     
                     fw.write(nog+": ");
                     for(String s:f){
                         fw.write(s+" ");
                     }
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        public void countGOOccurence(File geneFunctionFile, double perc ,File outputFile){
             BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneFunctions=new HashMap<>();
              HashMap<String, HashSet<String>> NOGFunctions=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneFunctionFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneFunctions.containsKey(gene)){
                                geneFunctions.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> func=geneFunctions.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String go=tmp[i].trim();
                                
                               func.add(go);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneFunctions.put(gene, func);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
            
            Iterator<String> it=NOGProtein.keySet().iterator();
            
            int numOKNOGS=0;
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<Pair<String,String>> genes=NOGProtein.get(nog);
                
                HashSet<String> allGenesWithFunction= new HashSet<String>();
                HashSet<String> allFunctions=new HashSet<>();
                HashSet<String> assignedFunctions=new HashSet<>();
                
                for(Pair<String,String> g:genes){
                    HashSet<String> f=geneFunctions.get(g.getValue0());
                    if(!geneFunctions.containsKey(g.getValue0()))
                        continue;
                    else if(f.size()>0){ allGenesWithFunction.add(g.getValue0());
                        allFunctions.addAll(f);
                    }
                }
                
                //System.out.println(nog+" "+genes.size()+" "+allFunctions.size());
                
                for(String fun:allFunctions){
                    int funCount=0;
                    for(String gen:allGenesWithFunction){
                        HashSet<String> f=geneFunctions.get(gen);
                        if(f.contains(fun))
                            funCount++;
                    }
                    
                    if((funCount/(double)allGenesWithFunction.size())>=perc)
                        assignedFunctions.add(fun);  
                }
                
                //System.out.println(nog+" "+assignedFunctions.size());
                NOGFunctions.put(nog, assignedFunctions); 
                if(assignedFunctions.size()>0)
                    numOKNOGS++;
            } 
            
            System.out.println("Fraction NOGs with function: "+(numOKNOGS/(double)NOGFunctions.keySet().size()));
            System.out.println("Number NOGs with function: "+numOKNOGS);
            
            HashMap<String,Integer> goCount=new HashMap<>();
            
            it=NOGFunctions.keySet().iterator();
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<String> s = NOGFunctions.get(nog);
                for(String g:s){
                    if(!goCount.containsKey(g))
                        goCount.put(g, 1);
                    else{
                        int c=goCount.get(g);
                        c++;
                        goCount.put(g, c);
                    }
                }
            }
            
            
            FileWriter fw = null;
            
              try{
                  int num=0;
                fw = new FileWriter(outputFile); 
                
                
                 it=goCount.keySet().iterator();
                 
                 while(it.hasNext()){
                     String go=it.next();
                     
                    int c=goCount.get(go);
                    
                    if(c>0)
                        num++;
                    
                     fw.write(go+" "+c);
                     fw.write("\n");
                 }
                System.out.println("Number of GOs "+num);
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        
        public void countGOOccWeighted(File OGOcc, File OGFuncFile){
                   
            OGGOMapping oggomap=new OGGOMapping();
                    oggomap.createCOGGOMapping(OGFuncFile);
                     BufferedReader reader;
                     
                     HashMap<String,Integer> goCount = new HashMap<>();
                     HashMap<String,Integer> ogCount = new HashMap<>();
      
             try {
                     Path path =Paths.get(OGOcc.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            ogCount.put(tmp[0].trim(), Integer.parseInt(tmp[1].trim()));
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
             
             Iterator<String> it = ogCount.keySet().iterator();
                     //oggomap.CogGOmap.keySet().iterator();
             
             while(it.hasNext()){
                 
                 String og = it.next();
                 if(!oggomap.CogGOmap.containsKey(og))
                     continue;
                 ArrayList<String> fun=oggomap.CogGOmap.get(og);
                 
                 for(int i=0;i<fun.size();i++){
                     if(!goCount.containsKey(fun.get(i)))
                         goCount.put(fun.get(i), 1);
                     else{
                         int c = goCount.get(fun.get(i));
                         c++;
                         goCount.put(fun.get(i), c);
                     }
                 }
                 
             }  
             
             
            FileWriter fw = null;
            
              try{
                  int num=0;
                fw = new FileWriter("C:\\Users\\matej\\Downloads\\Eukaryot data\\goCountNewnonCHSP.txt"); 
                
                
                 it=goCount.keySet().iterator();
                 
                 while(it.hasNext()){
                     String go=it.next();
                     
                    int c=goCount.get(go);
                    
                     fw.write(go+" "+c);
                     fw.write("\n");
                 }
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }   
             
        }
        
         public void countNOGOccurenceNew(File inputFolder, String taxIDFilePath ,File geneOgFile, File ensembleTIDsFile , HashMap<Integer,Integer> taxIDTranslationMap, CreateReducedMappingsFile rmf ,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, int useNC ,File outputFile){//add tax translation mapping
                     BufferedReader reader;
                     HashMap<String,HashSet<Pair<String,String>>> geneOG=new HashMap<>();
                     HashSet<Integer> usedTIDs=new HashSet<>();
                     Mappings map= new Mappings();
                      try {
                     Path path =Paths.get(geneOgFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length>2)
                                continue;
                            //tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneOG.containsKey(gene)){
                                geneOG.put(gene, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> ogs=geneOG.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].split(":")[0].trim();
                                Pair<String,String> p = new Pair(og,tmp[i].split(":")[1].trim());
                               ogs.add(p);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOG.put(gene, ogs);
                    }
                      }
                      catch(Exception e){e.printStackTrace();}  
                      
                 
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
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
              
         String extensions[] = {"dat"};
         boolean recursive = true; 
          HashSet<Integer> eggnogTaxs = new HashSet();
          HashSet<Integer> ensembleTaxs = new HashSet();
          
          BufferedReader reader1;
          
          try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
         
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
                      
          HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>(); 
            HashSet<Integer> usedTids = new HashSet();
            HashMap<String,Integer> taxFNMapping = new HashMap<>();
            
                       try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
       //  BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 
                   if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                     String fileName = input.getName();
                       String fileName1=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName1);
               
                    //int taxid=taxIDMap.get(fileName);  
                     
                    int found=0;
                      
                      if(taxIDTranslationMap.containsKey(taxid)){
                if(eggnogTaxs.contains(taxid)){found=1;}
                else if(eggnogTaxs.contains(taxIDTranslationMap.get(taxid)) && taxIDTranslationMap.get(taxid)==taxid && ensembleTaxs.contains(taxid)) found=1;
                else if(!ensembleTaxs.contains(taxIDTranslationMap.get(taxid)) && eggnogTaxs.contains(taxIDTranslationMap.get(taxid))) found=1;
                else if(!ensembleTaxs.contains(taxIDTranslationMap.get(taxid)) && !eggnogTaxs.contains(taxIDTranslationMap.get(taxid))){

                    for(int td:eggnogTaxs){
                        if(!taxIDTranslationMap.containsKey(td))
                            continue;
                        if(taxIDTranslationMap.get(td).equals(taxIDTranslationMap.get(taxid))){
                            found=1;
                            break;
                        }
                  }
                    
                }
            }
                      else continue;
                      
                      if(taxIDTranslationMap.get(taxid) == 4932 && taxid!=764097)
                          found=0;
                      else if(taxid==764097) found=1;
                      
                System.out.println("found: "+found);
                if(found==0 || (usedSpecies.contains(taxIDTranslationMap.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
                
               // if(taxIDTranslationMap.get(taxid)!=4896)//4932 S.Cerv, 4896//S.Pombe
                // if(taxIDTranslationMap.get(taxid)!=367110 /*taxid != 764097*/ /*taxid != 227321*/)
                    //        continue;
                
                 Gene cg=new Gene();
                 
                  cg.findCogsNewMetazoa(input,geneOGMap,map,taxid); //change to 0 to get COGs from file
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(cg.anotGen.size()<10+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
                           
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              int spec;
                              if(taxIDTranslationMap.containsKey(taxonId)){
                                  spec=taxIDTranslationMap.get(taxonId);
                                  /*if(usedTids.contains(spec))
                                      break;
                                  else usedTids.add(spec);*/
                              }
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
              System.out.println("taxid: "+taxid);
              System.out.println("taxonID: "+taxonId);
             usedTIDs.add(taxid);
             usedSpecies.add(taxIDTranslationMap.get(taxid));
             taxFNMapping.put(fileName, taxid);
          if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;
             //System.out.println("genes size: "+genes.size());
             int numProcG=0, np=0;
             for(String s:genes){
                 
                HashSet<Pair<String,String>> ogs=geneOG.get(s);
                
                if(ogs==null)
                    continue;
                
                for(Pair<String,String> n:ogs){   
                    if(!nogCount.containsKey(n.getValue0())){
                             nogCount.put(n.getValue0(), new HashSet<Integer>());
                             
                             if(!taxIDTranslationMap.containsKey(Integer.parseInt(n.getValue1())))
                                 continue;
                             
                             if(taxIDTranslationMap.get(Integer.parseInt(n.getValue1())).equals(taxIDTranslationMap.get(taxonId)))
                                         nogCount.get(n.getValue0()).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(n.getValue0());
                             
                              if(!taxIDTranslationMap.containsKey(Integer.parseInt(n.getValue1())))
                                 continue;
                         
                          if(taxIDTranslationMap.get(Integer.parseInt(n.getValue1())).equals(taxIDTranslationMap.get(taxonId)))
                             ti.add(taxonId);
                          
                             nogCount.put(n.getValue0(), ti);
                         }   
                }
                
              }  
                 np++;
                // System.out.println("genes remaining: "+(genes.size()-np));
             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);
         //    System.out.println("genes size: "+genes.size());
             int numProcG=0,np=0;
              for(String s:genes){
                 
                 HashSet<Pair<String,String>> ogs=geneOG.get(s);
                
                if(ogs==null)
                    continue;
                
                for(Pair<String,String> n:ogs){   
                    if(!nogCount.containsKey(n.getValue0())){
                             nogCount.put(n.getValue0(), new HashSet<Integer>());
                             
                             if(!taxIDTranslationMap.containsKey(Integer.parseInt(n.getValue1())))
                                 continue;
                             
                             if(taxIDTranslationMap.get(Integer.parseInt(n.getValue1())).equals(taxIDTranslationMap.get(taxonId)))
                                         nogCount.get(n.getValue0()).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(n.getValue0());
                             
                              if(!taxIDTranslationMap.containsKey(Integer.parseInt(n.getValue1())))
                                 continue;
                         
                          if(taxIDTranslationMap.get(Integer.parseInt(n.getValue1())).equals(taxIDTranslationMap.get(taxonId)))
                             ti.add(taxonId);
                          
                             nogCount.put(n.getValue0(), ti);
                         }   
                }
              }  

             taxCountMap.put(taxonId, count1);
         }

         System.out.println("Fraction genes with function: "+(geneswithFunction.size()/(double)genes.size()));
         System.out.println("Num genes: "+genes.size());
         System.out.println("Genes with function: "+geneswithFunction.size());
         
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
               /* Iterator<String>*/ it=nogCount.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<Integer> f=nogCount.get(nog);
                     
                     fw.write(nog+" "+f.size());
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }
              
             // usedTIDs.add(taxonId);
             //taxFNMapping.put(fileName, taxonId);
               
                  try{
                      
                      String utf= "usedTaxFile.txt";
                      
                     if(useNC==1)
                         utf="usedTaxFilenonCH";
                      
                fw = new FileWriter(utf); 
                
               /* Iterator<String>*/ it=taxFNMapping.keySet().iterator();
                 
                 while(it.hasNext()){
                     String fn=it.next();
                     
                     int f=taxFNMapping.get(fn);
                     
                     fw.write(fn+" "+f);
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }
                 
                 
                     try{
                         
                  String utf= "usedTaxes.txt";
                      
                     if(useNC==1)
                         utf="usedTaxesnonCH";        
                         
                fw = new FileWriter(utf); 
                
                Iterator<Integer> it1=usedTIDs.iterator();
                 
                 while(it1.hasNext()){
                     Integer tid=it1.next();
                     
                     fw.write(taxIDTranslationMap.get(tid)+"");
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }             
          }
        
          public void countNOGOccurence(File inputFolder, File geneOgFile ,File outputFile){
                     BufferedReader reader;
                     HashMap<String,HashSet<String>> geneOG=new HashMap<>();
                     
                      try {
                     Path path =Paths.get(geneOgFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length>2)
                                continue;
                            tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneOG.containsKey(gene)){
                                geneOG.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOG.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].trim();
                                
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOG.put(gene, ogs);
                    }
                      }
                      catch(Exception e){e.printStackTrace();}     
              
         String extensions[] = {"dat"};
         boolean recursive = true;             
                      
          HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>(); 
            
                       try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
             
                        
          if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;
             //System.out.println("genes size: "+genes.size());
             int numProcG=0, np=0;
             for(String s:genes){
                 
                HashSet<String> ogs=geneOG.get(s);
                
                if(ogs==null)
                    continue;
                
                for(String n:ogs){
                     if(!nogCount.containsKey(n)){
                             nogCount.put(n, new HashSet<Integer>());
                             nogCount.get(n).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(n);
                             ti.add(taxonId);
                             nogCount.put(n, ti);
                         }   
                }
                
              }  
                 np++;
                // System.out.println("genes remaining: "+(genes.size()-np));
             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);
         //    System.out.println("genes size: "+genes.size());
             int numProcG=0,np=0;
              for(String s:genes){
                 
                HashSet<String> ogs=geneOG.get(s);
                
                if(ogs==null)
                    continue;
                
                for(String n:ogs){
                     if(!nogCount.containsKey(n)){
                             nogCount.put(n, new HashSet<Integer>());
                             nogCount.get(n).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(n);
                             ti.add(taxonId);
                             nogCount.put(n, ti);
                         }   
                } 
              }  

             taxCountMap.put(taxonId, count1);
         }

         System.out.println("Fraction genes with function: "+(geneswithFunction.size()/(double)genes.size()));
         System.out.println("Num genes: "+genes.size());
         System.out.println("Genes with function: "+geneswithFunction.size());
         
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                Iterator<String> it=nogCount.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<Integer> f=nogCount.get(nog);
                     
                     fw.write(nog+" "+f.size());
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
                  
          }
        
        public void countNOGOccurence(File inputFolder, File geneFunctionFile, double perc, File mappingFile ,File outputFile){
             BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneFunctions=new HashMap<>();
              HashMap<String, HashSet<String>> NOGFunctions=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneFunctionFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            tmp[0]=tmp[0].replace(":", "");
                            String gene = tmp[0].trim();
                            
                            if(!geneFunctions.containsKey(gene)){
                                geneFunctions.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> func=geneFunctions.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String go=tmp[i].trim();
                                
                               func.add(go);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneFunctions.put(gene, func);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
            
            Iterator<String> it=NOGProtein.keySet().iterator();
            
            int numOKNOGS=0;
            
            while(it.hasNext()){
                String nog=it.next();
                HashSet<Pair<String,String>> genes=NOGProtein.get(nog);
                
                HashSet<String> allGenesWithFunction= new HashSet<String>();
                HashSet<String> allFunctions=new HashSet<>();
                HashSet<String> assignedFunctions=new HashSet<>();
                
                for(Pair<String,String> g:genes){
                    HashSet<String> f=geneFunctions.get(g.getValue0());
                    if(!geneFunctions.containsKey(g.getValue0()))
                        continue;
                    else if(f.size()>0){ allGenesWithFunction.add(g.getValue0());
                        allFunctions.addAll(f);
                    }
                }
                
                //System.out.println(nog+" "+genes.size()+" "+allFunctions.size());
                
                for(String fun:allFunctions){
                    int funCount=0;
                    for(String gen:allGenesWithFunction){
                        HashSet<String> f=geneFunctions.get(gen);
                        if(f.contains(fun))
                            funCount++;
                    }
                    
                    if((funCount/(double)allGenesWithFunction.size())>=perc)
                        assignedFunctions.add(fun);  
                }
                
                //System.out.println(nog+" "+assignedFunctions.size());
                NOGFunctions.put(nog, assignedFunctions); 
                if(assignedFunctions.size()>0)
                    numOKNOGS++;
            } 
            
            System.out.println("Fraction NOGs with function: "+(numOKNOGS/(double)NOGFunctions.keySet().size()));
            System.out.println("Number NOGs with function: "+numOKNOGS);
            
            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>();
            
            Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
         if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;
             //System.out.println("genes size: "+genes.size());
             int numProcG=0, np=0;
             for(String s:genes){
                 
                 HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s))
                     continue;
                 
                 for(Pair<String,String> s1:matches){
                 
                 if(geneFunctions.containsKey(s1.getValue0())){
                     if(geneFunctions.get(s1.getValue0()).size()>0)
                         geneswithFunction.add(s1.getValue0());
                     else
                         continue;
                 }
                 
                 if(!ProteinNOG.containsKey(s1.getValue0()))
                     continue;
                 
                 HashSet<Pair<String,String>> ns= ProteinNOG.get(s1.getValue0());
                 if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     if(!nogCount.containsKey(nss.getValue0())){
                             nogCount.put(nss.getValue0(), new HashSet<Integer>());
                             nogCount.get(nss.getValue0()).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(nss.getValue0());
                             ti.add(taxonId);
                             nogCount.put(nss.getValue0(), ti);
                         }   
                 }
              }  
                 np++;
                // System.out.println("genes remaining: "+(genes.size()-np));
             }

             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);
         //    System.out.println("genes size: "+genes.size());
             int numProcG=0,np=0;
              for(String s:genes){
                  
                  HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s))
                     continue;
                  
                 for(Pair<String,String> s1:matches){ 
                   if(geneFunctions.containsKey(s1.getValue0())){
                     if(geneFunctions.get(s1.getValue0()).size()>0)
                         geneswithFunction.add(s1.getValue0());
                     else 
                         continue;
                 }
 
                   
                   if(!ProteinNOG.containsKey(s1))
                     continue;
                   
                  HashSet<Pair<String,String>> ns= ProteinNOG.get(s1);
 
                 if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     if(!nogCount.containsKey(nss.getValue0())){
                             nogCount.put(nss.getValue0(), new HashSet<Integer>());
                             nogCount.get(nss.getValue0()).add(taxonId);
                         }
                          else{
                             HashSet<Integer> ti=nogCount.get(nss.getValue0());
                             ti.add(taxonId);
                             nogCount.put(nss.getValue0(), ti);
                         }   
                 }
               } 
                  np++;
                 //System.out.println("genes remaining: "+(genes.size()-np));   
             }

             taxCountMap.put(taxonId, count1);
         }

         System.out.println("Fraction genes with function: "+(geneswithFunction.size()/(double)genes.size()));
         System.out.println("Num genes: "+genes.size());
         System.out.println("Genes with function: "+geneswithFunction.size());
         
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 it=nogCount.keySet().iterator();
                 
                 while(it.hasNext()){
                     String nog=it.next();
                     
                     HashSet<Integer> f=nogCount.get(nog);
                     
                     fw.write(nog+" "+f.size());
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        
        
         public void nonChromosomalTaxIDs(File inputFolder, HashMap<Integer,Integer> taxTranslation ,File outputFile){
             BufferedReader reader;

         String extensions[] = {"dat"};
         boolean recursive = true;
         HashSet<Integer> OKTIDs = new HashSet();        
         
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {

                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                  //if(input.getName().contains(".nonchromosomal."))
                   //  continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              OKTIDs.add(taxonId);
                              if(taxTranslation.containsKey(taxonId))
                              taxonSpecies = taxTranslation.get(taxonId);
                              break;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 Iterator<Integer> it=OKTIDs.iterator();
                 
                 while(it.hasNext()){
                     fw.write(it.next()+"\n");
                     }
                   
                    // fw.write("\n");
     
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
         
         
          public void createGeneOGMappingNew(File inputFolder, File mappingFile, /*File ensembleTIDsFile,*/ HashMap<Integer,Integer> taxTranslation, int useNC ,File outputFile){
              //add checks for matching eggnog and ensemble organisms, only one with a given species id
              //for each ensemble taxID -> check if found in eggnog
              //get all taxIDs with the same species id and check if found in eggnog
              //check if species taxID found in eggnog
              //add taxID information to geneOGmapping.
              
              /*HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader;
             
             
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
                  
              }*/
              
              HashMap<String,HashSet<Pair<String,String>>> geneOGs=new HashMap<>();

            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>();
            
            Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                System.out.println("geneOGMap size: ");
                System.out.println(geneOGs.size());
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 
               if(input.getName().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              if(taxTranslation.containsKey(taxonId))
                              taxonSpecies = taxTranslation.get(taxonId);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(gen.contains(".")){
                                  String tmp[] = gen.split("\\.");
                                  
                                  if(tmp.length==2 && gen.trim().matches("SP[^.]+\\.[^.]+")){ 
                                      gen = gen.trim();
                                      System.out.println("gen spec: "+gen);
                                  }
                                  else if(tmp.length==2)
                                      gen=tmp[0].trim();
                                  else if(tmp.length==3)
                                      gen=tmp[1].trim();
                              }
                                  
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
         if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;

             int numProcG=0, np=0;
             for(String s:genes){
                 
                 HashSet<Pair<String,String>> matches = map.mappings.get(s);


                 if(!map.mappings.containsKey(s) && !ProteinNOG.containsKey(s))
                     continue;
                 else if(!map.mappings.containsKey(s) && ProteinNOG.containsKey(s)){
                         if(ProteinNOG.containsKey(s)){
                      
                  if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                           
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s);
                      HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     int eggSp1=-1;
                     
                      if(taxTranslation.containsKey(taxonId))
                              eggSp1 = taxTranslation.get(taxonId);
                      
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp1==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                 }
                         continue;
                 }
           
                 System.out.println("Gene: "+s);//add changed check
                 for(Pair<String,String> s1:matches){
                 System.out.println(" "+s1.getValue0());
                 if(ProteinNOG.containsKey(s1.getValue0())){
                     
                     
                      int eggSp=-1;
                     if(taxTranslation.containsKey(Integer.parseInt(s1.getValue1())))
                              eggSp = taxTranslation.get(Integer.parseInt(s1.getValue1()));
                     
                     if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies && eggSp!=taxonSpecies)
                         continue;
                     
                         
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<Pair<String,String>>());
                     }
                     
                     HashSet<Pair<String,String>> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                            og.add(new Pair(s3.getValue0(),taxonId));
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                     
                  HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                     
                     geneOGs.put(s, og1);
                 }
                 else{

                 }
                System.out.println("\n");    
              }  

             }

             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);

             int numProcG=0,np=0;
              for(String s:genes){
                  
                  HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s) && !ProteinNOG.containsKey(s))
                     continue;
                 else if(!map.mappings.containsKey(s) && ProteinNOG.containsKey(s)){
                         if(ProteinNOG.containsKey(s)){
                             
                     if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }        
                             
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s);
                      HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     int eggSp1=-1;

                      if(taxTranslation.containsKey(taxonId))
                              eggSp1 = taxTranslation.get(taxonId);
                      
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp1==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                 }
                         continue;
                 }
           
                  
                  System.out.println("Gene: "+s);
                 for(Pair<String,String> s1:matches){ 
                     System.out.println(" "+s1.getValue0());
                  if(ProteinNOG.containsKey(s1.getValue0())){
                            
                        int eggSp=-1;
                     if(taxTranslation.containsKey(Integer.parseInt(s1.getValue1())))
                              eggSp = taxTranslation.get(Integer.parseInt(s1.getValue1()));
                      
                        if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies && eggSp!=taxonSpecies)
                             continue;
                        
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<Pair<String,String>>());
                     }
                     
                     HashSet<Pair<String,String>> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                            og.add(new Pair(s3.getValue0(),taxonId));
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                     
                  HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                     
                     geneOGs.put(s, og1);
                     
                 }
                  else{
                  }
               } 
                 System.out.println("\n");
             }
             taxCountMap.put(taxonId, count1);
         }
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 Iterator<String> it=geneOGs.keySet().iterator();
                 
                 while(it.hasNext()){
                     String gene=it.next();
                     
                     HashSet<Pair<String,String>> f=geneOGs.get(gene);
                     int count=f.size();
              
                     if(count>0){
                     fw.write(gene+" ");
                                        
                     for(Pair<String,String> s1:f){
                         if(count>1){
                              fw.write(String.valueOf(s1.getValue0())+":"+String.valueOf(s1.getValue1())+" ");
                              count--;
                         }
                         else fw.write(String.valueOf(s1.getValue0())+":"+String.valueOf(s1.getValue1())+"\n");
                     }
                 }
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
          
          
           public void createGeneOGMappingNewMetazoa(File inputFolder, File mappingFile, /*File ensembleTIDsFile,*/ HashMap<Integer,Integer> taxTranslation, int useNC ,File outputFile){
              
              HashMap<String,HashSet<Pair<String,String>>> geneOGs=new HashMap<>();

            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>();
            
            Mappings map=new Mappings();
            map.loadMappingsMetazoa(mappingFile); //create metazoa version
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                System.out.println("geneOGMap size: ");
                System.out.println(geneOGs.size());
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 
               if(input.getName().contains(".nonchromosomal.") && useNC==0)
                     continue;

                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0, proteinSection=0;
                       int wrong=0;
                       String gen="";
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              if(taxTranslation.containsKey(taxonId))
                              taxonSpecies = taxTranslation.get(taxonId);
                          }
                           else if(line.contains("gene") && !line.contains("=") && line.contains("..") && (taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183) /*&& taxonId!=6238 && taxonId!=9595 && taxonId!=9606 && taxonId!=69293*/){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          else if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
                         //28377, 9913, 9595, 9606, 69293, 6239, 6238, 9483, 9615, 60711, 7719, 7955, 7227, 46245, 7245, 9796, 9685, 9031, 69293, 9544, 9103, 13616, 10090, 9258, 9986, 8090, 9598, 9601, 10116, 9823, 59729, 99883   
                         // gene - 7165, 7460, 7425, 6183, 
                          if(line.contains("/protein_id") && proteinSection==1){
                               gen = line.split("=")[1].trim();
                              gen = gen.replaceAll("\"", "");
                              gen = gen.replace("-PA", "");
                              gen = gen.replace("-PB", "");
                              gen = gen.replace("-tr", "");
                              if(taxonId!=6183 && taxonId!=6239)
                                  gen = (gen.split("\\.")[0]).trim();
                              genes.add(gen);
                              proteinSection=0;
                          }
                          
                          if(line.contains("/gene=") && geneSection==1){
                               gen=line.split("=")[1].trim();
                               gen = gen.replaceAll("\"", "");
                              gen = gen.replace("-PA", "");
                              gen = gen.replace("-PB", "");
                              gen = gen.replace("-tr", "");
                              if(taxonId!=6183)
                                  gen = (gen.split("\\.")[0]).trim();
                              //System.out.println("gene: "+gen);
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
         if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;

             int numProcG=0, np=0;
             for(String s:genes){
                 HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 int ind=0;
                 if(map.mappings.containsKey(s+"-PA"))
                     ind=1;
                 if(map.mappings.containsKey(s+"-PB"))
                     ind=2;
                 if(map.mappings.containsKey(s+"-tr"))
                     ind=3;
                
                if(ind==0){ 
                 if(!map.mappings.containsKey(s) && !ProteinNOG.containsKey(s))
                     continue;
                 else if(!map.mappings.containsKey(s) && ProteinNOG.containsKey(s)){
                         if(ProteinNOG.containsKey(s)){
                      
                  if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                           
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s);
                      HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     int eggSp1=-1;
                     
                      if(taxTranslation.containsKey(taxonId))
                              eggSp1 = taxTranslation.get(taxonId);
                      
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp1==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                 }
                         continue;
                 }
                }
                else{
                    if(ind==1)
                          matches = map.mappings.get(s+"-PA");
                 
                    else if(ind==2)
                         matches = map.mappings.get(s+"-PB");
                    else if(ind==3)
                         matches = map.mappings.get(s+"-tr");
                }
                 System.out.println("Gene: "+s);//add changed check
                 for(Pair<String,String> s1:matches){
                 System.out.println(" "+s1.getValue0());
                 
                 String tmp = s1.getValue0();
                 tmp=tmp.replace("-PA", "");
                 tmp=tmp.replace("-PB", "");
                 tmp=tmp.replace("-tr", "");
                 s1=s1.setAt0(tmp);
                 if(ProteinNOG.containsKey(s1.getValue0())){
                     
                     
                      int eggSp=-1;
                     if(taxTranslation.containsKey(Integer.parseInt(s1.getValue1())))
                              eggSp = taxTranslation.get(Integer.parseInt(s1.getValue1()));
                     
                     if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies && eggSp!=taxonSpecies)
                         continue;
                     
                         
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<Pair<String,String>>());
                     }
                     
                     HashSet<Pair<String,String>> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                            og.add(new Pair(s3.getValue0(),taxonId));
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                     
                  HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                     
                     geneOGs.put(s, og1);
                 }
                 else{

                 }
                System.out.println("\n");    
              }  

             }

             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);

             int numProcG=0,np=0;
              for(String s:genes){
                   int ind=0;
                   
                 if(map.mappings.containsKey(s+"-PA"))
                     ind=1;
                 if(map.mappings.containsKey(s+"-PB"))
                     ind=2;
                 if(map.mappings.containsKey(s+"-tr"))
                     ind=3;
                 
                  HashSet<Pair<String,String>> matches = map.mappings.get(s);
               
                  if(ind==0){
                 if(!map.mappings.containsKey(s) && !ProteinNOG.containsKey(s))
                     continue;
                 else if(!map.mappings.containsKey(s) && ProteinNOG.containsKey(s)){
                         if(ProteinNOG.containsKey(s)){
                             
                     if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }        
                             
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s);
                      HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     int eggSp1=-1;

                      if(taxTranslation.containsKey(taxonId))
                              eggSp1 = taxTranslation.get(taxonId);
                      
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp1==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                 }
                         continue;
                 }
                  }
                  else{
                    if(ind==1)
                          matches = map.mappings.get(s+"-PA");
                 
                    else if(ind==2)
                         matches = map.mappings.get(s+"-PB");
                    else if(ind==3)
                         matches = map.mappings.get(s+"-tr");
                }
           
                  
                  System.out.println("Gene: "+s);
                 for(Pair<String,String> s1:matches){ 
                     System.out.println(" "+s1.getValue0());
                      String tmp = s1.getValue0();
                 tmp=tmp.replace("-PA", "");
                 tmp=tmp.replace("-PB", "");
                 tmp=tmp.replace("-tr", "");
                 s1=s1.setAt0(tmp);
                  if(ProteinNOG.containsKey(s1.getValue0())){
                            
                        int eggSp=-1;
                     if(taxTranslation.containsKey(Integer.parseInt(s1.getValue1())))
                              eggSp = taxTranslation.get(Integer.parseInt(s1.getValue1()));
                      
                        if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies && eggSp!=taxonSpecies)
                             continue;
                        
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<Pair<String,String>>());
                     }
                     
                     HashSet<Pair<String,String>> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                            og.add(new Pair(s3.getValue0(),taxonId));
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<Pair<String,String>>());
                     }
                     
                  HashSet<Pair<String,String>> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp==taxonSpecies)
                         og1.add(new Pair(s3.getValue0(),taxonId));
                     
                     geneOGs.put(s, og1);
                     
                 }
                  else{
                  }
               } 
                 System.out.println("\n");
             }
             taxCountMap.put(taxonId, count1);
         }
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 Iterator<String> it=geneOGs.keySet().iterator();
                 
                 while(it.hasNext()){
                     String gene=it.next();
                     
                     HashSet<Pair<String,String>> f=geneOGs.get(gene);
                     int count=f.size();
              
                     if(count>0){
                     fw.write(gene+" ");
                                        
                     for(Pair<String,String> s1:f){
                         if(count>1){
                              fw.write(String.valueOf(s1.getValue0())+":"+String.valueOf(s1.getValue1())+" ");
                              count--;
                         }
                         else fw.write(String.valueOf(s1.getValue0())+":"+String.valueOf(s1.getValue1())+"\n");
                     }
                 }
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
          
        
         public void createGeneOGMapping(File inputFolder, File mappingFile, HashMap<Integer,Integer> taxTranslation ,File outputFile){
             BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneOGs=new HashMap<>();

            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            HashMap<Integer,HashSet<String>> taxGeneMap = new HashMap<>();
            HashMap<String,HashSet<Integer>> nogCount= new HashMap<>();
            
            Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                System.out.println("geneOGMap size: ");
                System.out.println(geneOGs.size());
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              if(taxTranslation.containsKey(taxonId))
                              taxonSpecies = taxTranslation.get(taxonId);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
         if(!taxCountMap.containsKey(taxonId)){
             int count=0, count1=0;
             //System.out.println("genes size: "+genes.size());
             int numProcG=0, np=0;
             for(String s:genes){
                 
                 HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s))
                     continue;
           
                 System.out.println("Gene: "+s);//add changed check
                 for(Pair<String,String> s1:matches){
                 System.out.println(" "+s1.getValue0());
                 if(ProteinNOG.containsKey(s1.getValue0())){
                     
                     int eggSp=-1;
                     if(taxTranslation.containsKey(Integer.parseInt(s1.getValue1())))
                              eggSp = taxTranslation.get(Integer.parseInt(s1.getValue1()));
                     
                     if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies && eggSp != taxonSpecies)
                         continue;
                         
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<String>());
                     }
                     
                     HashSet<String> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp == taxonSpecies)
                            og.add(s3.getValue0());
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<String>());
                     }
                     
                  HashSet<String> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies || eggSp == taxonSpecies)
                         og1.add(s3.getValue0());
                     
                     geneOGs.put(s, og1);
                 }
                 else{
                    /* Iterator<String> it1=ProteinNOG.keySet().iterator();
                     
                     while(it1.hasNext()){
                         String g=it1.next();
                         if(!map.mappings.containsKey(g))
                             continue;
                         HashSet<String> matches1=map.mappings.get(g);

                             if(matches1.contains(s1)){
                                 HashSet<String> ogs=ProteinNOG.get(g);
                     
                                 if(!geneOGs.containsKey(g)){
                                        geneOGs.put(g, new HashSet<String>());
                                     }
                     
                                        HashSet<String> og=geneOGs.get(g);
                     
                                            for(String s3:ogs)
                                                      og.add(s3);
                                           }
                     }*/
                     //continue;
                 }
                System.out.println("\n");    
              }  
                // System.out.println("genes remaining: "+(genes.size()-np));
             }

             taxCountMap.put(taxonId, count1);   
         }
         else{
             int count=taxCountMap.get(taxonId), count1=taxCountMap.get(taxonId);
         //    System.out.println("genes size: "+genes.size());
             int numProcG=0,np=0;
              for(String s:genes){
                  
                  HashSet<Pair<String,String>> matches = map.mappings.get(s);
                 
                 if(!map.mappings.containsKey(s))
                     continue;
                  
                  System.out.println("Gene: "+s);
                 for(Pair<String,String> s1:matches){ 
                     System.out.println(" "+s1.getValue0());
                  if(ProteinNOG.containsKey(s1.getValue0())){
                                   
                        if(Integer.parseInt(s1.getValue1())!=taxonId && Integer.parseInt(s1.getValue1())!=taxonSpecies)
                         continue;
                        
                     HashSet<Pair<String,String>> ogs=ProteinNOG.get(s1.getValue0());
                     
                     if(!geneOGs.containsKey(s1.getValue0())){
                         geneOGs.put(s1.getValue0(), new HashSet<String>());
                     }
                     
                    /* HashSet<String> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         og.add(s3.getValue0());
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<String>());
                     }
                     
                    og=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         og.add(s3.getValue0());*/
                     
                     
                     
                     HashSet<String> og=geneOGs.get(s1.getValue0());
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies)
                            og.add(s3.getValue0());
                     
                      if(!geneOGs.containsKey(s)){
                         geneOGs.put(s, new HashSet<String>());
                     }
                     
                  HashSet<String> og1=geneOGs.get(s);
                     
                     for(Pair<String,String> s3:ogs)
                         if(Integer.parseInt(s3.getValue1())==taxonId || Integer.parseInt(s3.getValue1())==taxonSpecies)
                         og1.add(s3.getValue0());
                     
                     geneOGs.put(s, og1);
                     
                 }
                  else{
                       /*Iterator<String> it1=ProteinNOG.keySet().iterator();
                     
                     while(it1.hasNext()){
                         String g=it1.next();
                         if(!map.mappings.containsKey(g))
                             continue;
                         HashSet<String> matches1=map.mappings.get(g);

                             if(matches1.contains(s1)){
                                 HashSet<String> ogs=ProteinNOG.get(g);
                     
                                 if(!geneOGs.containsKey(g)){
                                        geneOGs.put(g, new HashSet<String>());
                                     }
                     
                                        HashSet<String> og=geneOGs.get(g);
                     
                                            for(String s3:ogs)
                                                      og.add(s3);
                                           }*/
                    // }
                     // continue;
                  }
               } 
                 System.out.println("\n");
             }
             taxCountMap.put(taxonId, count1);
         }
             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
            
            FileWriter fw = null;
            
              try{
                fw = new FileWriter(outputFile); 
                
                
                 Iterator<String> it=geneOGs.keySet().iterator();
                 
                 while(it.hasNext()){
                     String gene=it.next();
                     
                     HashSet<String> f=geneOGs.get(gene);
                     
                     fw.write(gene+" ");
                     int count=f.size();
                     
                     for(String s1:f){
                         if(count>1){
                              fw.write(s1+" ");
                              count--;
                         }
                         else fw.write(s1+"\n");
                     }
                   
                    // fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        
        public void writeGeneFunctionMapping(File inputFolder, File mappingFile, File uniprotMapping, File outputFile){
            Mappings map=new Mappings();
            map.loadMappings(mappingFile);
            map.loadUniprotMappings(uniprotMapping);
            
            HashMap<String,HashSet<Pair<String,String>>> geneFunctions = new HashMap<>();
            
             FileWriter fw = null;
           
              String extensions[] = {"dat"};
              boolean recursive = true;
             
             try{
                fw = new FileWriter(outputFile);       
                
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();

                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);

                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                       int geneSection=0, cdsSection=0;
                       int wrong=0;
                       int countTotal=0, countAssigned=0, countInEggnog=0;
                        String gen=""; HashSet<String> functions= new HashSet<>();
                      while ((line = reader1.readLine()) != null){//parse a .dat file
                         
                          //System.out.println(line);
                          //System.out.println("geneSection: "+geneSection);
                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              /*if(taxonId!=559292){
                                 wrong=1;
                                  break;
                              }*/
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          if(line.contains("/gene=") && geneSection==1){
                                 gen=line.split("=")[1].trim();
                                 
                                 if(!geneFunctions.containsKey(gen)){
                                        geneFunctions.put(gen, new HashSet<Pair<String,String>>());
                                 }
                                 
                          }
                          
                          if(line.contains("CDS") && !line.contains("=") && line.contains("..") )
                              cdsSection=1;
                          
                          if(line.contains(" /db_xref=") && cdsSection==1){
                                    if(line.contains("GO:")){
                                        String tmp[]=line.split("=");
                                         HashSet<Pair<String,String>> functs = geneFunctions.get(gen);
                                         functs.add(new Pair(tmp[1].trim(),taxonId));
                                        functions.add(tmp[1]);    
                                    }            
                          }
                          else if(line.contains("/translation") && cdsSection==1){
                                if(functions.size()>0){ 
                                        fw.write(gen.replaceAll("\"","")+": ");
                              
                              for(String fun:functions)
                                  fw.write(fun.replaceAll("\"", "")+" ");
                              
                              fw.write("\n");
                              
                              functions.clear();
                              cdsSection=0;
    
                                }
                          }
                      }
                      
                      System.out.println("Count total: "+countTotal);
                      System.out.println("Count in Eggnog total: "+countInEggnog);
                      System.out.println("CountAssigned: "+countAssigned);
                      System.out.println("Fraction with function: "+((double)countAssigned/countTotal));
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close();  
            }
               
                fw.close();
            }
            catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        public void writeGeneFunctionMappingNew(File inputFolder, File mappingFile, File uniprotMapping, int useNC ,File outputFile){
            Mappings map=new Mappings();
            map.loadMappings(mappingFile);
            map.loadUniprotMappings(uniprotMapping);
            
            HashMap<String,HashSet<Pair<String,String>>> geneFunctions = new HashMap<>();
            
             FileWriter fw = null;
           
              String extensions[] = {"dat"};
              boolean recursive = true;
             
             try{
                fw = new FileWriter(outputFile);       
                
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();

                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);

                  if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;
                 
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                       int geneSection=0, cdsSection=0;
                       int wrong=0;
                       int countTotal=0, countAssigned=0, countInEggnog=0;
                        String gen=""; HashSet<String> functions= new HashSet<>();
                      while ((line = reader1.readLine()) != null){//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                                 gen=line.split("=")[1].trim();
                                 
                                 if(gen.contains(".")){
                                  String tmp[] = gen.split("\\.");
                                  if(tmp.length==2)
                                      gen=tmp[0].trim();
                                  else if(tmp.length==3 && gen.trim().matches("SP[^.]+\\.[^.]+\\.\\d+")){
                                      gen=tmp[0].trim()+"."+tmp[1].trim();
                                  }
                                  else if(tmp.length==3)
                                      gen=tmp[1].trim();
                              }
                                 
                                 if(!geneFunctions.containsKey(gen)){
                                        geneFunctions.put(gen, new HashSet<Pair<String,String>>());
                                 }
                          }
                          
                          if(line.contains("CDS") && !line.contains("=") && line.contains("..") )
                              cdsSection=1;
                          
                          if(line.contains(" /db_xref=") && cdsSection==1){
                                    if(line.contains("GO:")){
                                        String tmp[]=line.split("=");
                                         HashSet<Pair<String,String>> functs = geneFunctions.get(gen);
                                         functs.add(new Pair(tmp[1].trim(),taxonId));
                                        functions.add(tmp[1]);    
                                    }            
                          }
                          else if(line.contains("/translation") && cdsSection==1){

                              cdsSection=0;       
                          }
                      }
                      

                      System.out.println("Count total: "+countTotal);
                      System.out.println("Count in Eggnog total: "+countInEggnog);
                      System.out.println("CountAssigned: "+countAssigned);
                      System.out.println("Fraction with function: "+((double)countAssigned/countTotal));
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close();  
            }
            
            Iterator<String> it = geneFunctions.keySet().iterator();
                              
                              while(it.hasNext()){
                              String gene = it.next();
                              HashSet<Pair<String,String>> functio = geneFunctions.get(gene);
                              int count=functio.size();
                              
                              if(count==0)
                                  continue;
                              
                               fw.write(gene+" ");
                               
                                for(Pair<String,String> p:functio){
                                    if(count>1){
                                        fw.write(String.valueOf(p.getValue0())+"::"+String.valueOf(p.getValue1())+" ");
                                        count--;
                                    }
                                    else {
                                       fw.write(String.valueOf(p.getValue0())+"::"+String.valueOf(p.getValue1())+"\n"); 
                                    }
                                }        
                  }

                fw.close();
            }
            catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        
         public void writeGeneFunctionMappingNewMetazoa(File inputFolder, File mappingFile, File uniprotMapping, int useNC ,File outputFile){
            Mappings map=new Mappings();
            map.loadMappingsMetazoa(mappingFile);
            //map.loadUniprotMappings(uniprotMapping);
            
            HashMap<String,HashSet<Pair<String,String>>> geneFunctions = new HashMap<>();
            
             FileWriter fw = null;
           
              String extensions[] = {"dat"};
              boolean recursive = true;
             
             try{
                fw = new FileWriter(outputFile);       
                
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();

                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);

                  if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;
                 
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                       int geneSection=0, cdsSection=0, proteinSection=0;
                       int wrong=0;
                       int countTotal=0, countAssigned=0, countInEggnog=0;
                        String gen=""; HashSet<String> functions= new HashSet<>();
                      while ((line = reader1.readLine()) != null){//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..") && (taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183) /*&& taxonId!=6238 && taxonId!=9595 && taxonId!=9606 && taxonId!=69293*/){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          else if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
                         //28377, 9913, 9595, 9606, 69293, 6239, 6238, 9483, 9615, 60711, 7719, 7955, 7227, 46245, 7245, 9796, 9685, 9031, 69293, 9544, 9103, 13616, 10090, 9258, 9986, 8090, 9598, 9601, 10116, 9823, 59729, 99883   
                         // gene - 7165, 7460, 7425, 6183, 
                          if(line.contains("/protein_id") && proteinSection==1){
                               gen = line.split("=")[1].trim();
                              gen = gen.replaceAll("\"", "");
                              gen = gen.replace("-PA", "");
                              gen = gen.replace("-tr", "");
                              if(taxonId!=6183 && taxonId!=6239)
                                  gen = gen.split("\\.")[0];
                             // genes.add(gen);
                              proteinSection=0;
                               if(!geneFunctions.containsKey(gen)){
                                        geneFunctions.put(gen, new HashSet<Pair<String,String>>());
                                 }
                          }
                          
                          if(line.contains("/gene=") && geneSection==1){
                               gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              //System.out.println("gene: "+gen);
                             // genes.add(gen);
                              geneSection=0;
                               if(!geneFunctions.containsKey(gen)){
                                        geneFunctions.put(gen, new HashSet<Pair<String,String>>());
                                 }
                          }
                          
                          if(line.contains("CDS") && !line.contains("=") && line.contains("..") )
                              cdsSection=1;
                          
                          if(line.contains(" /db_xref=") && cdsSection==1){
                                    if(line.contains("GO:")){
                                        String tmp[]=line.split("=");
                                         HashSet<Pair<String,String>> functs = geneFunctions.get(gen);
                                         functs.add(new Pair(tmp[1].trim(),taxonId));
                                        functions.add(tmp[1]);    
                                    }            
                          }
                          else if(line.contains("/translation") && cdsSection==1){

                              cdsSection=0;       
                          }
                      }
                      

                      System.out.println("Count total: "+countTotal);
                      System.out.println("Count in Eggnog total: "+countInEggnog);
                      System.out.println("CountAssigned: "+countAssigned);
                      System.out.println("Fraction with function: "+((double)countAssigned/countTotal));
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close();  
            }
            
            Iterator<String> it = geneFunctions.keySet().iterator();
                              
                              while(it.hasNext()){
                              String gene = it.next();
                              HashSet<Pair<String,String>> functio = geneFunctions.get(gene);
                              int count=functio.size();
                              
                              if(count==0)
                                  continue;
                              
                               fw.write(gene+" ");
                               
                                for(Pair<String,String> p:functio){
                                    if(count>1){
                                        fw.write(String.valueOf(p.getValue0())+"::"+String.valueOf(p.getValue1())+" ");
                                        count--;
                                    }
                                    else {
                                       fw.write(String.valueOf(p.getValue0())+"::"+String.valueOf(p.getValue1())+"\n"); 
                                    }
                                }        
                  }

                fw.close();
            }
            catch(IOException e){
              e.printStackTrace();
            }  
        }
        
        public void writeGenesEnsemble(File inputFolder, File output){
             String extensions[] = {"dat"};
         boolean recursive = true;
           HashSet<String> genes=new HashSet<>();
           
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
      
         BufferedReader reader1;
          
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                // System.out.println("File = " + input.getAbsolutePath());
               //  System.out.println("File = " + input.getName());
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                     
                       int geneSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                         
                          //System.out.println(line);
                          //System.out.println("geneSection: "+geneSection);
                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                             /* if(taxonId!=559292){
                                 wrong=1;
                                  break;
                              }*/
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              //System.out.println("gene: "+gen);
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                    /*  if(wrong==1){
                          wrong=0;
                          continue;
                      }*/
                     // System.out.println("num genes: "+genes.size());
                        reader1.close();       
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
         FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
                
                for(String s:genes)
                    fw.write(s+"\n");
               
                fw.close();
            }
            catch(IOException e){
              e.printStackTrace();
            }
        }
        
        
        public void writeGenesEnsembleMetazoa(File inputFolder, File output){
             String extensions[] = {"dat"};
         boolean recursive = true;
           HashSet<String> genes=new HashSet<>();
           
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
      
         BufferedReader reader1;
          
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                     
                       int geneSection=0, proteinSection=0;
                       int wrong=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                         
                          //System.out.println(line);
                          //System.out.println("geneSection: "+geneSection);
                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                             /* if(taxonId!=559292){
                                 wrong=1;
                                  break;
                              }*/
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..") && (taxonId==7165 || taxonId == 7460 || taxonId == 7425 || taxonId == 6183) /*&& taxonId!=6238 && taxonId!=9595 && taxonId!=9606 && taxonId!=69293*/){
                              geneSection=1;
                              //System.out.println("Gene section entered!");
                        }
                          else if(line.contains("CDS") && taxonId != 7165 && taxonId!=7460 && taxonId!=7425 && taxonId!=6183 /*&& (taxonId == 6239 || taxonId == 6238 || taxonId==9595 || taxonId == 9606 || taxonId == 69293)*/){
                              proteinSection=1;
                          }
                         //28377, 9913, 9595, 9606, 69293, 6239, 6238, 9483, 9615, 60711, 7719, 7955, 7227, 46245, 7245, 9796, 9685, 9031, 69293, 9544, 9103, 13616, 10090, 9258, 9986, 8090, 9598, 9601, 10116, 9823, 59729, 99883   
                         // gene - 7165, 7460, 7425, 6183, 
                          if(line.contains("/protein_id") && proteinSection==1){
                              String gen = line.split("=")[1].trim();
                              gen = gen.replaceAll("\"", "");
                              gen = gen.replace("-PA", "");
                              gen = gen.replace("-tr", "");
                              if(taxonId!=6183 && taxonId!=6239)
                                  gen = gen.split("\\.")[0];
                              genes.add(gen);
                              proteinSection=0;
                          }
                          
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              if(taxonId!=6183)
                                  gen = gen.split("\\.")[0];
                              //System.out.println("gene: "+gen);
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                    /*  if(wrong==1){
                          wrong=0;
                          continue;
                      }*/
                     // System.out.println("num genes: "+genes.size());
                        reader1.close();       
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
         FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
                
                for(String s:genes)
                    fw.write(s+"\n");
               
                fw.close();
            }
            catch(IOException e){
              e.printStackTrace();
            }
        }
        
        
        public void createReduced_mappingsEns(File mapping, File geneFile, File output){
         
            FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
             
            EnsembleGenes eGene=new EnsembleGenes();
            eGene.loadGenes(geneFile);
            
            BufferedReader reader;
            
                try{
                     Path path =Paths.get(mapping.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");

                            String id1=tmp[1].trim();
                            if(id1.contains(".")){
                                String tmpt[]=id1.split("\\.");
                                if(tmpt.length<2)
                                    continue;
                                id1=tmpt[0].trim();
                            }
                                
                            String id2=tmp[2].trim();
                             if(id2.contains(".")){
                                String tmpt[]=id2.split("\\.");
                                if(tmpt.length<2){
                                    continue;
                                   /* System.out.println("id1: "+id1);
                                    System.out.println("id2: "+id2);
                                    System.out.println(line);*/
                                }
                                
                                id2=tmpt[0].trim();
                            }

                             
                             /*int tid = Integer.parseInt(pid1[0]);
                                if( tid == 7165 || tid == 7425 || tid == 7460){
                                    pid1[1] = pid1[1].replace("-PA", "");
                                }
                                else if(tid == 6183){
                                    pid1[1] = pid1[1].replace("__mRNA", "");
                                }
                                else if(tid == 6239){
                                    pid1[1] = pid1[1]+"."+pid1[2];
                                }*/
                             
                            if(proteinsEggnog.contains(id1) || proteinsEggnog.contains(id1.toLowerCase()) || proteinsEggnog.contains(id2.toLowerCase()) || proteinsEggnog.contains(id2) || eGene.genes.contains(id1) || eGene.genes.contains(id2)){
                                fw.write(line+"\n");

                               // HashSet<String> containedPIDsInMap=new HashSet<>();
                                if(proteinsEggnog.contains(id1))
                                    containedPIDsInMap.add(id1);
                                else if(proteinsEggnog.contains(id2)) containedPIDsInMap.add(id2);
                                
                                // System.out.println(id1+" "+id2);

                            String mapType=tmp[3];
                            String types[] = mapType.split("_");
                           if(types.length<2)
                               continue;
                            if(/*types[0].toLowerCase().contains("blast") &&*/ types[1].toLowerCase().contains("uniprot")){
                              //  System.out.println("conditions met1");
                                String uni=tmp[2];
                                if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else if(proteinsEggnog.contains(id2))
                                    assignedPIDsToUID.add(id2);
                            }
                            else if(types.length==2){
                                    if(types[0].toLowerCase().contains("uniprot") /*&& (types[1].toLowerCase().contains("blast"))*/){
                                   
                              // System.out.println("conditions met2");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else if(proteinsEggnog.contains(id2))
                                    assignedPIDsToUID.add(id2);
                            } 
                        }
                            else if(types.length==3){
                                 if(types[0].toLowerCase().contains("uniprot") /*&& (types[1].toLowerCase().contains("blast") || types[2].toLowerCase().contains("blast")/*)*/){
                               //      System.out.println("conditions met3");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else if(proteinsEggnog.contains(id2))
                                    assignedPIDsToUID.add(id2);
                                 }
                            }
                      }
                    }
                fw.close();
                reader.close();
                }
                catch(Exception e){e.printStackTrace();}    
                
                System.out.println("Fraction of eggnog pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/proteinsEggnog.size()));
                 System.out.println("Fraction of eggnog pids contained in mapping: "+((double)containedPIDsInMap.size()/proteinsEggnog.size()));
        }
        
        
        
        public void createReduced_mappings(File mapping, File output){
         
            FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
            
            BufferedReader reader;
            
                try{
                     Path path =Paths.get(mapping.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                      HashSet<String> extendedMapping = new HashSet<>();
                      int useExtendedMapping=1;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");

                            String id1=tmp[1].trim();
                            String id2=tmp[2].trim();

                            if(proteinsEggnog.contains(id1) || proteinsEggnog.contains(id2) || extendedMapping.contains(id1) || extendedMapping.contains(id2)){
                                fw.write(line+"\n");
                                
                                if(useExtendedMapping==1){
                                extendedMapping.add(id1); extendedMapping.add(id2);
                                }
                               // HashSet<String> containedPIDsInMap=new HashSet<>();
                                if(proteinsEggnog.contains(id1))
                                    containedPIDsInMap.add(id1);
                                else containedPIDsInMap.add(id2);
                                
                                // System.out.println(id1+" "+id2);

                            String mapType=tmp[3];
                            String types[] = mapType.split("_");
                           
                            if(types[0].toLowerCase().contains("blast") && types[1].toLowerCase().contains("uniprot")){
                              //  System.out.println("conditions met1");
                                String uni=tmp[2];
                                if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                            }
                            else if(types.length==2){
                                    if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast"))){
                                   
                              // System.out.println("conditions met2");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                            } 
                        }
                            else if(types.length==3){
                                 if(types[0].toLowerCase().contains("uniprot") && (types[1].toLowerCase().contains("blast") || types[2].toLowerCase().contains("blast"))){
                               //      System.out.println("conditions met3");
                                String uni=tmp[1];
                                 if(uni.contains("_"))
                                    uni=uni.split("_")[0];
                                uniprot_acs.add(uni.trim());
                                if(proteinsEggnog.contains(id1))
                                        assignedPIDsToUID.add(id1);
                                else
                                    assignedPIDsToUID.add(id2);
                                 }
                            }
                      }
                    }
                fw.close();
                reader.close();
                }
                catch(Exception e){e.printStackTrace();}    
                
                System.out.println("Fraction of eggnog pids assigned to uniprotID: "+((double)assignedPIDsToUID.size()/proteinsEggnog.size()));
                 System.out.println("Fraction of eggnog pids contained in mapping: "+((double)containedPIDsInMap.size()/proteinsEggnog.size()));
        }
        
        public void createReduced_uniprot(File uniprot, File output){
             FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
            
            BufferedReader reader;
            
                try{
                     Path path =Paths.get(uniprot.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            if(tmp.length<2)
                                continue;
                            String id1=tmp[1].trim();
                            
                            if(uniprot_acs.contains(id1))
                                fw.write(line+"\n");      
                    }
                fw.close();
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
        }
        
        public void filterNOGs(File input, File output){
             FileWriter fw = null;
            
             try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
            
            BufferedReader reader;
            
                try{
                     Path path =Paths.get(input.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length<2)
                                continue;

                            fw.write(line+"\n");      
                    }
                fw.close();
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
        }
        
        public void propagateFunctions(File input, String ontologyPath ,File output){
            
            HashMap<String,HashSet<String>> nogGOs=new HashMap<>();
            
              FileWriter fw = null;
   
            BufferedReader reader;
            
                try{
                     Path path =Paths.get(input.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length<2)
                                continue;
                            
                            String nog= tmp[0].replace(":", "").trim();
                            
                            nogGOs.put(nog, new HashSet<String>());
                            
                            for(int i=1;i<tmp.length;i++){
                                if(tmp[i].trim().equals("GO:0008150") || tmp[i].trim().equals("GO:0003674") || tmp[i].trim().equals("GO:0005575"))
                                    continue;
                                nogGOs.get(nog).add(tmp[i].trim());  
                            }
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
            
            
             OntologyTools.GeneOntology myGO=null;

             try{
                   myGO= new OntologyTools.GeneOntology(ontologyPath);
                 }
             catch(IOException e){
                      e.printStackTrace();
                 }

             
                   try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
             
           Iterator<String> it=nogGOs.keySet().iterator();
           
           while(it.hasNext()){
             String nog=it.next();
            HashSet<String> gos=nogGOs.get(nog);
            HashSet<Integer> extendedGOs=new HashSet<>();

            for(String goS:gos){
                String tmp[]=goS.split(":");
                int go=Integer.parseInt(tmp[1]);
                 Collection<OntologyTools.GOTerm> parents = myGO.get(go).getAllParents();
             for ( OntologyTools.GOTerm curPar : parents ) {
                     // if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
                             extendedGOs.add(curPar.getId());
                }
                extendedGOs.add(go);
             }
            
            
            for(int g:extendedGOs){
                String GOName="GO:";
                 int nNumGo=g;
                 if(nNumGo == 8150 || nNumGo == 3674 || nNumGo == 5575)
                     continue;
             int numTmpGo=nNumGo, nDigGo=0;
              
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=nNumGo;
              nogGOs.get(nog).add(GOName);
            }    
           }
           
           it=nogGOs.keySet().iterator();
           try{
                 while(it.hasNext()){
                      String nog=it.next();
                     int size=nogGOs.get(nog).size();
                     
                     if(size==0)
                         continue;
                     
                    
                    fw.write(nog+": ");
                 
                  
                    
                 for(String s:nogGOs.get(nog)){
                     if(size>1)
                        fw.write(s+" ");
                     else if(size==1)
                         fw.write(s+"\n");
                     size--;
                 }
              }
              fw.close();
           }
           catch(IOException e){
               e.printStackTrace();
           } 
        }
        
    public void finalNOGSelection(File NOGsFilt, File NOGsRed, File output){
            HashSet<String> NOgred=new HashSet<>();
            
            BufferedReader reader;
             FileWriter fw = null;
             
                try{
                     Path path =Paths.get(NOGsRed.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                           NOgred.add(line.trim());
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
        
                 try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }

                try{
                     Path path =Paths.get(NOGsFilt.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                           String tmp[]=line.split(" ");
                           tmp[0]=tmp[0].replace(":", "").trim();
                           if(!NOgred.contains(tmp[0]))
                               continue;
                           else fw.write(line+"\n");
                    }
                reader.close();
                fw.close();
                }
                catch(Exception e){e.printStackTrace();}          
        }
    
    public void countGOFunctions(File input, File output){
        
        HashMap<String,Integer> GOCount=new HashMap<>();
            
            BufferedReader reader;
             FileWriter fw = null;
             
                try{
                     Path path =Paths.get(input.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                           String tmp[]=line.split(" ");
                           tmp[0]=tmp[0].replace(":", "").trim();
                           for(int i=1;i<tmp.length;i++){
                               String go=tmp[i].trim();
                               if(!GOCount.containsKey(go))
                                   GOCount.put(go, 1);
                               else{
                                   int c=GOCount.get(go);
                                   c++;
                                   GOCount.put(go, c);
                               }
                           }
                               
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
        
                 try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }

                try{
                    Iterator<String> it=GOCount.keySet().iterator();
                    while(it.hasNext()){
                        String go=it.next();
                        int count=GOCount.get(go);
                        fw.write(go+" "+count+"\n");
                    }
                fw.close();
                }
                catch(Exception e){e.printStackTrace();}    
        
    }
    
    public void filterAndReduceGOs(File nogGO, File goRed, File output){
           HashSet<Integer> GOred=new HashSet<>();
            
            BufferedReader reader;
             FileWriter fw = null;
             
                try{
                     Path path =Paths.get(goRed.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                          GOred.add(Integer.parseInt(line.trim().replace("GO:", "")));
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
        
                 try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }

                try{
                     Path path =Paths.get(nogGO.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                     String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                         
                           String tmp[]=line.split(" ");
                           tmp[0]=tmp[0].replace(":", "").trim();
                           
                           HashSet<String> gos=new HashSet<>();
                           
                           for(int i=1;i<tmp.length;i++){
                               if(GOred.contains(Integer.parseInt(tmp[i].trim().replace("GO:", ""))))
                                   gos.add(tmp[i].trim());
                           }
                           
                           if(gos.size()<1)
                               continue;
                           
                           fw.write(tmp[0]+": ");
                           
                           int size=gos.size();
                           
                           for(String s:gos){
                               if(size>1)
                               fw.write(s+" ");
                               else
                                   fw.write(s+"\n");
                               size--;
                           }
                           
                           /*if(!NOgred.contains(tmp[0]))
                               continue;
                           else fw.write(line+"\n");*/
                    }
                reader.close();
                fw.close();
                }
                catch(Exception e){e.printStackTrace();}  
    }
    
     public void countCoveredLocations(File inputFolder, File geneOGFile, File mappingFile){
             BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneOgs=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneOGFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOgs.containsKey(gene)){
                                geneOgs.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOgs.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].trim();
                                
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOgs.put(gene, ogs);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
           
            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
           Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 if(input.getName().contains(".nonchromosomal."))
                     continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0; int locOK=0, locEx=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                               locEx++;
                              if(geneOgs.containsKey(gen))
                                  locOK++;
                              else{
                                  if(!map.mappings.containsKey(gen))
                                      continue;
                                  for(Pair<String,String> s:candidates)
                                      if(geneOgs.containsKey(s.getValue0())){
                                          locOK++;
                                          break;
                                      }
                              }
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
                        System.out.println("Perc locations OK: "+((double)locOK/(double)locEx));

             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  

        }
      
      public void countCoveredLocationsNew(File inputFolder, File geneOGFile, HashMap<Integer,Integer> taxTranslation ,File mappingFile){
             BufferedReader reader;

              HashMap<String,HashSet<String>> geneOgs=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneOGFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            if(tmp.length>2)
                                continue;
                            String gene = tmp[0].trim();
                            
                            if(!geneOgs.containsKey(gene)){
                                geneOgs.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOgs.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].trim();
                                
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOgs.put(gene, ogs);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
           
            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
           Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 if(input.getName().contains(".nonchromosomal."))
                     continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0; int locOK=0, locEx=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              if(taxTranslation.containsKey(taxonId))
                              taxonSpecies = taxTranslation.get(taxonId);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                               locEx++;
                              if(/*ProteinNOG.containsKey(gen)*/geneOgs.containsKey(gen)){
                                  HashSet<Pair<String,String>> d=ProteinNOG.get(gen);
                                  
                                  /*for(Pair<String,String> p:d){
                                      if(Integer.parseInt(p.getValue1())== taxonId){
                                          locOK++;
                                          break;
                                      }
                                  }*/
                                  
                                  locOK++;
                              }
                              /*else{
                                  if(!map.mappings.containsKey(gen))
                                      continue;
                                  for(Pair<String,String> s:candidates)
                                      if(geneOgs.containsKey(s.getValue0())){
                                          locOK++;
                                          break;
                                      }
                              }*/
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
                        System.out.println("Perc locations OK: "+((double)locOK/(double)locEx));

             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  

        }
      
      public HashSet<Integer> countEggnogTax(File eggnogFolder, HashMap<Integer,Integer> translationMap){
          HashSet<String> taxId=new HashSet();
          
           BufferedReader reader;
            
             try {
                     Path path =Paths.get(eggnogFolder.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            String protIds[]=tmp[5].split(",");
                           
                            for(int i=0;i<protIds.length;i++){
                                String pid=protIds[i].trim();
                                String pid1[]=pid.split("\\.");
                                taxId.add(pid1[0].trim());
                            }      
                            //System.out.println("Mapping: "+" "+NOG+" "+NOGProtein.get(NOG).size());
                    }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
          
          System.out.println("Num taxIDs: "+taxId.size());
          
          HashSet<Integer> taxIdNumeric = new HashSet();
          
          for(String s:taxId)
              taxIdNumeric.add(Integer.parseInt(s));
          
          
          HashSet<Integer> taxSpecies=new HashSet();
          
          for(int i:taxIdNumeric)
              if(translationMap.containsKey(i)){
              taxSpecies.add(translationMap.get(i));
              System.out.println(i);
              }
          
          System.out.println("Num species taxIDs: "+taxSpecies.size());
          
          return taxIdNumeric;
      }
     
     
      public void countTaxSpecies(File inputFolder, HashMap<Integer,Integer> taxTranslation,HashSet<Integer> eggnogTaxes){
             BufferedReader reader;

            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
            
            HashSet<Integer> taxIdEnsemble= new HashSet();
            HashSet<Integer> taxIdSpecies = new HashSet();
            
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 if(input.getName().contains(".nonchromosomal."))
                     continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1, taxonSpecies=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0; int locOK=0, locEx=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "").trim();
                              taxonId=Integer.parseInt(tn[1]);
                              taxIdEnsemble.add(taxonId);
                              /*if(taxTranslation.containsKey(taxonId)){
                              taxonSpecies = taxTranslation.get(taxonId);
                              taxIdSpecies.add(taxonSpecies);
                              break;
                              }*/
                               taxIdSpecies.add(taxTranslation.get(taxonId));
                               break;
                          }
                    }
                        reader1.close(); 

             fCount++;      
            }
            
             int overlap1=0, overlap2=0,overlap3=0;        
              
            System.out.println("Tax ids not found in ensemble files");
                for(int sp:eggnogTaxes){
                    if(taxIdEnsemble.contains(sp))
                        overlap1++;
                    if(taxIdSpecies.contains(sp)){
                        overlap2++;
                        System.out.println("Species: "+sp);
                    }
                    if(taxIdEnsemble.contains(sp) || taxIdSpecies.contains(sp))
                        overlap3++;
                    if(!taxIdEnsemble.contains(sp) && !taxIdSpecies.contains(sp))
                        System.out.println(sp);
                }
                
                System.out.println("Overlap on subspecies level: "+overlap1);
                System.out.println("Overlap on species level: "+overlap2);
                System.out.println("Overlap on either level: "+overlap3);
                
                HashSet<Integer> eggnogSpecies = new HashSet<>();
                
                for(int sp:eggnogTaxes){
                
                 if(taxTranslation.containsKey(sp)){
                              int taxS = taxTranslation.get(sp);
                              eggnogSpecies.add(taxS);
                 }
               }
                
                System.out.println("Number on species level in eggnog: "+eggnogSpecies.size());
                       
        } catch (Exception e) {
            e.printStackTrace();
        }  
}
          
      public void createFunctionalNeighbourhoodDataset(File inputFolder, File geneOGFile, File mappingFile, File ogFuncFile){
             BufferedReader reader;
              
              HashMap<String,HashSet<String>> geneOgs=new HashMap<>();
            
             try {
                     Path path =Paths.get(geneOGFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOgs.containsKey(gene)){
                                geneOgs.put(gene, new HashSet<String>());
                            }
                            
                            HashSet<String> ogs=geneOgs.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og=tmp[i].trim();
                                
                               ogs.add(og);
                            }
                           // System.out.println("GF: "+gene+" "+func.size());
                            geneOgs.put(gene, ogs);
                    }
         //System.out.println("Num pids: "+proteinsEggnog.size());
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
            
             //System.out.println("Gene-function file loaded");
           
            HashMap<Integer,Integer> taxCountMap = new HashMap<>();
            
         String extensions[] = {"dat"};
         boolean recursive = true;
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;
           Mappings map=new Mappings();
            map.loadMappings(mappingFile); 
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                 if(input.getName().contains(".nonchromosomal."))
                     continue;
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                      HashSet<String> genes=new HashSet<>();
                      HashSet<String> geneswithFunction=new HashSet<>();
                       int geneSection=0;
                       int wrong=0; int locOK=0, locEx=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                          }
                          else if(line.contains("gene") && !line.contains("=") && line.contains("..")){
                              geneSection=1;
                        }
                          if(line.contains("/gene=") && geneSection==1){
                              String gen=line.split("=")[1].trim();
                              HashSet<Pair<String,String>> candidates=map.mappings.get(gen);
                               locEx++;
                              if(geneOgs.containsKey(gen))
                                  locOK++;
                              else{
                                  if(!map.mappings.containsKey(gen))
                                      continue;
                                  for(Pair<String,String> s:candidates)
                                      if(geneOgs.containsKey(s.getValue0())){
                                          locOK++;
                                          break;
                                      }
                              }
                              genes.add(gen);
                              geneSection=0;
                          }
                    }
                      if(wrong==1){
                          wrong=0;
                          continue;
                      }

                        reader1.close(); 
                        System.out.println("Genes computed...");
                        System.out.println("Perc locations OK: "+((double)locOK/(double)locEx));

             fCount++;      
            }
        } catch (Exception e) {
            e.printStackTrace();
        }  
      }  
}
