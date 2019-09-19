/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;
import org.javatuples.Triplet;
import static recursivefilesearch.COG.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about gene-OG mappings
 */
 
public class GenesAndOgsOrganism {
    public HashSet<String> organismOgs;
    public HashMap<String,HashSet<String>> geneOGMappingOrganism;
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
        geneOGMappingOrganism.putAll(c.geneOGMappingOrganism);
        OGGeneMappingOrganism.putAll(c.OGGeneMappingOrganism);
    }
    
    public void traverseGenomes(File folderPath, int taxIdT, String[] extensions, Boolean recursive, String taxIDFilePath, GeneOGMapping geneOGMap){
         try{
                Collection files = FileUtils.listFiles(folderPath, extensions, recursive);            
               
                BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
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
               
                for (Iterator iterator = files.iterator(); iterator.hasNext();) { 
                File input = (File) iterator.next();

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                if(taxid!=taxIdT)
                    continue;
                 
                getGenesAndOgs(input.getAbsolutePath(),geneOGMap);
                
                }
                
           }
           catch(Exception e){
               e.printStackTrace();
           }
    }
    
    public void getGenesAndOgs(String input, GeneOGMapping geneOGmap){

        int lineNum=0, minCoordinate, maxCoordinate;
        HashSet<String> genes = new HashSet<>();

        File inF = new File(input);
        
        BufferedReader reader;
         try {
      Path path =Paths.get(inF.getAbsolutePath());
      System.out.println("Path: "+inF.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          
          if(lineNum==1){
              String[] gensize=line.split(" ");
              String c=gensize[gensize.length-1];
              System.out.println("String c COG: "+c);
              String [] coord=c.split("\\.\\.");
              System.out.println("coord size: "+coord.length);
              minCoordinate=Integer.parseInt(coord[0]);
              maxCoordinate=Integer.parseInt(coord[1]);
          }
          
          if(lineNum>3){
          String[] st=line.split("\t");
          String[] coord=st[0].split("\\..");
          int x=Integer.parseInt(coord[0]), y=Integer.parseInt(coord[1]);

          String geneName="";

             geneName=st[3];
         genes.add(geneName);
         if(geneOGmap.geneOGsMap.containsKey(geneName)){
                organismOgs.addAll(geneOGmap.geneOGsMap.get(geneName));
         HashSet<String> ogs = new HashSet();
         ogs.addAll(geneOGmap.geneOGsMap.get(geneName));
         if(!geneOGMappingOrganism.containsKey(geneName))
             geneOGMappingOrganism.put(geneName, ogs);
         else{
             geneOGMappingOrganism.get(geneName).addAll(ogs);
         }
         
         for(String og:geneOGmap.geneOGsMap.get(geneName))
                if(!OGGeneMappingOrganism.containsKey(og)){
                    OGGeneMappingOrganism.put(og,new HashSet<String>());
                    OGGeneMappingOrganism.get(og).add(geneName);
                }
                else{
                    OGGeneMappingOrganism.get(og).add(geneName);
                }
            }
          }
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
        
    }
