/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create Caffa validation set
 */
public class CreateCaffaValidationSetTotal {
    static public void main(String[] args){
        
        File inputProkaryot = new File("C:\\Users\\matej\\Downloads\\Phenotypes and data Barcelona\\CAFA-2013-targets\\TargetListProkarya.txt");
        File inputMetazoa = new File("C:\\Users\\matej\\Downloads\\Phenotypes and data Barcelona\\CAFA-2013-targets\\TargetListMetazoa.txt");
        File inputFungy = new File("C:\\Users\\matej\\Downloads\\Phenotypes and data Barcelona\\CAFA-2013-targets\\TargetListFungi.txt");
        File inputGeneInfo = new File("C:\\Users\\matej\\Downloads\\Phenotypes and data Barcelona\\CAFA-2013-targets\\caffaSeq.fasta.emapper.annotations");
        
        HashSet<String> prokaryotGenes = new HashSet<>();
        HashSet<String> metazoaGenes = new HashSet<>();
        HashSet<String> fungyGenes = new HashSet<>();
        
        HashMap<String, Pair<String,ArrayList<String>>> prokaryot = new HashMap<>(); 
        HashMap<String, Pair<String,ArrayList<String>>> fungi = new HashMap<>();
        HashMap<String, Pair<String,ArrayList<String>>> metazoa = new HashMap<>();
        
        Path p=Paths.get(inputProkaryot.getAbsolutePath());
        
        try{
        BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
      
        String line = "";
        
        while((line = read.readLine())!=null){
            prokaryotGenes.add(line.trim());
         }
            read.close();
            }
        catch(IOException e){
            e.printStackTrace();
        }
        
         p=Paths.get(inputFungy.getAbsolutePath());
        
        try{
        BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
      
        String line = "";
        
        while((line = read.readLine())!=null){
            fungyGenes.add(line.trim());
         }
            read.close();
            }
        catch(IOException e){
            e.printStackTrace();
        }
    
     p=Paths.get(inputMetazoa.getAbsolutePath());
        
    try{
        BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
      
        String line = "";
        
        while((line = read.readLine())!=null){
            metazoaGenes.add(line.trim());
         }
            read.close();
            }
        catch(IOException e){
            e.printStackTrace();
        }
    
        p = Paths.get(inputGeneInfo.getAbsolutePath());
        
        try{
            BufferedReader read = Files.newBufferedReader(p, StandardCharsets.UTF_8);
            String line = "";
            
            while((line = read.readLine())!=null){
                String tmp[] = line.split("\t");
                
                String gene = tmp[0].trim();
                String ename = tmp[1].trim();
                String ogs = "";
         
                for(int i=0;i<tmp.length;i++){
                    if(tmp[i].contains("@")){
                        ogs = tmp[i];
                        break;
                    }
                }
                     
                if(ogs.equals(""))
                    continue;
                
                if(prokaryotGenes.contains(gene)){
                    String t1[] = ogs.split(",");
                    ArrayList<String> cogs = new ArrayList<>();
                    for(String s:t1){
                        if(s.contains("COG")){
                                String t[] = s.split("@");
                                cogs.add(t[0].trim());
                            }
                    }
                    Pair<String,ArrayList<String>> pp = new Pair(ename,cogs);
                            if(!prokaryot.containsKey(gene))
                                prokaryot.put(gene, pp);
                }
                else if(fungyGenes.contains(gene)){
                    String t1[] = ogs.split(",");
                    ArrayList<String> cogs = new ArrayList<>();
                    for(String s:t1){
                        if(s.contains("fuNOG")){
                                String t[] = s.split("@");
                                cogs.add("ENOG41"+t[0].trim());
                            }
                    }
                    Pair<String,ArrayList<String>> pp = new Pair(ename,cogs);
                            if(!fungi.containsKey(gene))
                                fungi.put(gene, pp);
                }
                else if(metazoaGenes.contains(gene)){
                    String t1[] = ogs.split(",");
                    ArrayList<String> cogs = new ArrayList<>();
                    for(String s:t1){
                        if(s.contains("meNOG")){
                                String t[] = s.split("@");
                                cogs.add("ENOG41"+t[0].trim());
                            }
                    }
                    Pair<String,ArrayList<String>> pp = new Pair(ename,cogs);
                            if(!metazoa.containsKey(gene))
                                metazoa.put(gene, pp);
                }
            }
            
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        System.out.println("Prokaryot "+prokaryot.keySet().size());
        System.out.println("Fungi "+fungi.keySet().size());
        System.out.println("Metazoa "+metazoa.keySet().size());
        
        File outputProkaryot = new File("prokaryotInfo.txt");
        File outputFungi = new File("fungiInfo.txt");
        File outputMetazoa = new File("metazoaInfo.txt");
        
        try{
              FileWriter fw = new FileWriter(outputProkaryot);
              
              Iterator<String> it = prokaryot.keySet().iterator();
              
              while(it.hasNext()){
                  String cafName = it.next();
                  Pair<String,ArrayList<String>> pp = prokaryot.get(cafName);
                  
                  String out = cafName+"\t"+pp.getValue0()+"\t";
                  
                  ArrayList<String> ogs = pp.getValue1();
                  
                  if(ogs.isEmpty())
                      continue;
                  
                  for(int i=0;i<ogs.size();i++){
                      if((i+1)<ogs.size())
                         out+=ogs.get(i)+"\t";
                      else out+=ogs.get(i)+"\n";
                  }
                  fw.write(out);
              }
              fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        try{
              FileWriter fw = new FileWriter(outputFungi);
              
              Iterator<String> it = fungi.keySet().iterator();
              
              while(it.hasNext()){
                  String cafName = it.next();
                  Pair<String,ArrayList<String>> pp = fungi.get(cafName);
                  
                  String out = cafName+"\t"+pp.getValue0()+"\t";
                  
                  ArrayList<String> ogs = pp.getValue1();
                  
                  if(ogs.isEmpty())
                      continue;
                  
                  for(int i=0;i<ogs.size();i++){
                      if((i+1)<ogs.size())
                         out+=ogs.get(i)+"\t";
                      else out+=ogs.get(i)+"\n";
                  }  
                  fw.write(out);
              }
              fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
         try{
              FileWriter fw = new FileWriter(outputMetazoa);
              
              Iterator<String> it = metazoa.keySet().iterator();
              
              while(it.hasNext()){
                  String cafName = it.next();
                  Pair<String,ArrayList<String>> pp = metazoa.get(cafName);
                  
                  String out = cafName+"\t"+pp.getValue0()+"\t";
                  
                  ArrayList<String> ogs = pp.getValue1();
                  
                  if(ogs.isEmpty())
                      continue;
                  
                  for(int i=0;i<ogs.size();i++){
                      if((i+1)<ogs.size())
                         out+=ogs.get(i)+"\t";
                      else out+=ogs.get(i)+"\n";
                  }  
                  fw.write(out);
              }
              fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
}
