/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
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
  *@mail matmih1@gmail.com
  *@description class to store information about synthenic blocks in eukaryotic organisms
 */
public class SynthenicBlocks {
    public HashMap<Integer,ArrayList<Pair<Integer,Integer>>> blocks;
    
    
    public SynthenicBlocks(){
        blocks = new HashMap<>();
    }
    
    
    public void readBlocks(File inputFolder, String [] extensions, HashMap<String,Integer> taxIDMap){
          
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions , true);

        
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                      File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                   String fileName=input.getName().substring(0,input.getName().length()-4);
                   
                   Path p = Paths.get(input.getAbsolutePath());
                   BufferedReader read = Files.newBufferedReader(p, StandardCharsets.UTF_8);
                   
                   String line = "";
                   
                   int fileMode = 0;
                   
                   if(input.getName().contains("Vertebrata_"))
                       fileMode = 1;
                   
                   if(input.getName().contains("mm10_UCNE_orthologs_mouse")){
                       fileMode = 2;
                   }
                   
                   if(input.getName().contains("GalGalhg19")){
                       fileMode = 3;
                   }
                   
                   int count =0;
                   while((line = read.readLine())!=null){
                       
                        if(fileMode == 0){
                            if(count ==0 ){
                                count = 1;
                                continue;
                            }
                            
                            String tmp[] = line.split(",");
                            
                            if(!line.contains(","))
                                tmp=line.split("\t");
                            
                            Pair<Integer,Integer> pPair = new Pair(Integer.parseInt(tmp[5].trim()),Integer.parseInt(tmp[6].trim()));
                            int taxID = taxIDMap.get(fileName);
                            
                            if(!blocks.containsKey(taxID))
                                   blocks.put(taxID, new ArrayList<Pair<Integer,Integer>>());
                            blocks.get(taxID).add(pPair);
                        }
                        else if(fileMode == 1){
                             String tmp[] = line.split("\t");
                              Pair<Integer,Integer> pPair = new Pair(Integer.parseInt(tmp[1].trim()),Integer.parseInt(tmp[2].trim()));
                            int taxID = taxIDMap.get(fileName);
                            
                            if(!blocks.containsKey(taxID))
                                   blocks.put(taxID, new ArrayList<Pair<Integer,Integer>>());
                            blocks.get(taxID).add(pPair);
                        }
                        else if(fileMode == 2){
                            
                            if(count ==0 ){
                                count = 1;
                                continue;
                            }
                            
                             String tmp[] = line.split(",");
                            
                            if(!line.contains(","))
                                tmp=line.split("\t");
                            
                            Pair<Integer,Integer> pPair = new Pair(Integer.parseInt(tmp[4].trim()),Integer.parseInt(tmp[5].trim()));
                            int taxID = taxIDMap.get(fileName);
                            
                            if(!blocks.containsKey(taxID))
                                   blocks.put(taxID, new ArrayList<Pair<Integer,Integer>>());
                            blocks.get(taxID).add(pPair);
                            
                        }
                        else if(fileMode == 3){
                            
                            if(count ==0 ){
                                count = 1;
                                continue;
                            }
                            
                             String tmp[] = line.split("\t");
                            
                            Pair<Integer,Integer> pPair = new Pair(Integer.parseInt(tmp[1].trim()),Integer.parseInt(tmp[2].trim()));
                            int taxID = taxIDMap.get(fileName);
                            
                            if(!blocks.containsKey(taxID))
                                   blocks.put(taxID, new ArrayList<Pair<Integer,Integer>>());
                            blocks.get(taxID).add(pPair);
                            
                        }
                   }
                   
                   read.close();
              }
            }
          catch(IOException e){
              e.printStackTrace();
          }
    }
    
}
