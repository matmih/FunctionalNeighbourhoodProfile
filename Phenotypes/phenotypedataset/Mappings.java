/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import static phenotypedataset.LoadGeneData.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class storing different gene-OG, cog, GO index mappings.
 */
public class Mappings {
    
    HashMap<String,HashSet<String>> geneOgMapping = new HashMap<>();
    HashMap<String, HashSet<String>> ogFunctionMapping = new HashMap<>();
    HashMap<String,Integer> GOtoIndex = new HashMap<>();
    HashMap<Integer,String> indexToGO = new HashMap<>();
    HashMap<Integer,String> cogIndex = new HashMap<>();
    HashSet<String> allGOs = new HashSet<>();
    
    public void loadGeneFunctionMapping(File input){
        int count=0;
           BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int countOG=0;
      
      while ((line = reader.readLine()) != null) {
        String tmp[] = line.split(":");
        String og = tmp[0].trim();
        
        if(!ogFunctionMapping.containsKey(tmp[0].trim())){
            ogFunctionMapping.put(tmp[0].trim(), new HashSet<>());
            cogIndex.put(countOG++, tmp[0].trim());
        }
          
        String gos[] = tmp[1].split(",");
        
        for(int i=0;i<gos.length;i++){
            String go = gos[i].trim();            
            ogFunctionMapping.get(og).add(gos[i].trim());
            allGOs.add(go);
            if(!GOtoIndex.containsKey(go)){
                GOtoIndex.put(go, count++);
                indexToGO.put(count-1, go);
            }
        }
        
      }
         }
         catch(Exception e){
             e.printStackTrace();
         }
         
         System.out.println("Num GOs:");
         System.out.println(allGOs.size());      
    }
    
    public void loadGeneOGMapping(File input){
         BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      while ((line = reader.readLine()) != null) {
        String tmp[] = line.split(":");
        String gene = tmp[0].trim();
        
        if(!geneOgMapping.containsKey(tmp[0].trim()))
            geneOgMapping.put(tmp[0].trim(), new HashSet<>());
          
        String ogs[] = tmp[1].split(",");
        
        for(int i=0;i<ogs.length;i++){
            geneOgMapping.get(gene).add(ogs[i].trim());
        }
        
      }
         }
         catch(Exception e){
             e.printStackTrace();
         }
         
         System.out.println("Num Genes:");
         System.out.println(geneOgMapping.keySet().size());
    }
    
}
