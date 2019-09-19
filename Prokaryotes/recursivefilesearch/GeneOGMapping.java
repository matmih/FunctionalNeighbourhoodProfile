/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about geneOGmapping
 */
public class GeneOGMapping {
    public HashMap<String,HashSet<String>> geneOGsMap=new HashMap<>();
    public HashMap<String,HashSet<String>> OGGenesMap=new HashMap<>();
 
    public GeneOGMapping(){
        
    }
    
    public GeneOGMapping(GeneOGMapping cp){
        Iterator<String> it;
        it=cp.geneOGsMap.keySet().iterator();
        
        while(it.hasNext()){
            String k=it.next();
            HashSet<String> tmp=cp.geneOGsMap.get(k);
            HashSet<String> tmp1=new HashSet<>();
            tmp1.addAll(tmp);
            geneOGsMap.put(k, tmp1);
                   } 
    }
    
    public HashSet<String> getGenes(String Og){
        
        Iterator<String> geneIt= geneOGsMap.keySet().iterator();
        HashSet<String> associatedGenes=new HashSet<>();
        
        while(geneIt.hasNext()){
            String gene=geneIt.next();
            HashSet<String> ogs=geneOGsMap.get(gene);
            if(ogs.contains(Og))
                      associatedGenes.add(gene);
        }
        return associatedGenes;
    }
    
    public void reduceMap(ReducedOgsGOs red){
        
        Iterator<String> it=geneOGsMap.keySet().iterator();
        
        while(it.hasNext()){
            String g=it.next();
            HashSet<String> ogs=geneOGsMap.get(g);
            
            Iterator<String> it1=ogs.iterator();
            
           while(it1.hasNext()){
               String og=it1.next();
               if(!red.redCogs.contains(og))
                   it1.remove();
           }
           
           if(ogs.size()==0)
               it.remove();
           else{
               geneOGsMap.put(g, ogs);
           }
        }
    }
    
    
    public void loadGOGMapping(File input){

        BufferedReader reader;
         try {
                Path path =Paths.get(input.getAbsolutePath());

                reader = Files.newBufferedReader(path,ENCODING);
                String line = null;
                int lineC=0;
      while ((line = reader.readLine()) != null) {
            if(lineC==0){
                lineC=1;
                continue;
            }
            
          String[] st=line.split("\t");

          String geneName=st[0];
          String ogNames[]=st[2].split(",");
        
         for(int i=0;i<ogNames.length;i++){
          String ogName=ogNames[i];
          String GOName="";
          int ogNumber=Integer.parseInt(ogName);
          
          int numTmpGo=ogNumber, nDigGo=0;
          
          if(numTmpGo<0){
              numTmpGo=-numTmpGo;
              GOName+="NOG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<6){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=(-ogNumber);
              
              if(!geneOGsMap.containsKey(geneName))
                    geneOGsMap.put(geneName, new HashSet<String>());
              
              geneOGsMap.get(geneName).add(GOName);
              
          }
          else{
              GOName+="COG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<4){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=ogNumber;
              
              if(!geneOGsMap.containsKey(geneName))
                    geneOGsMap.put(geneName, new HashSet<String>());
              
              geneOGsMap.get(geneName).add(GOName);
          }
         }
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         }
    
    
    public void loadGOGMapping(File input, int flag){//flag allows computing COGgene mapping

        if(flag==1)
            OGGenesMap = new HashMap<>();
        
        BufferedReader reader;
         try {
                Path path =Paths.get(input.getAbsolutePath());

                reader = Files.newBufferedReader(path,ENCODING);
                String line = null;
                int lineC=0;
      while ((line = reader.readLine()) != null) {
            if(lineC==0){
                lineC=1;
                continue;
            }
            
          String[] st=line.split("\t");

          String geneName=st[0];
          String ogNames[]=st[2].split(",");
        
         for(int i=0;i<ogNames.length;i++){
          String ogName=ogNames[i];
          String OGName="";
          int ogNumber=Integer.parseInt(ogName);
          
          int numTmpGo=ogNumber, nDigGo=0;
          
          if(numTmpGo<0){
              numTmpGo=-numTmpGo;
              OGName+="NOG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<6){
                 OGName+="0";
                  nDigGo++;
              }
              
              OGName+=(-ogNumber);
              
              if(!geneOGsMap.containsKey(geneName))
                    geneOGsMap.put(geneName, new HashSet<String>());
              
              geneOGsMap.get(geneName).add(OGName);
              
              if(flag==1){
                  if(!OGGenesMap.containsKey(OGName))
                    OGGenesMap.put(OGName, new HashSet<String>());
                  OGGenesMap.get(OGName).add(geneName);
                          }             
          }
          else{
              OGName+="COG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<4){
                  OGName+="0";
                  nDigGo++;
              }
              
              OGName+=ogNumber;
              
              if(!geneOGsMap.containsKey(geneName))
                    geneOGsMap.put(geneName, new HashSet<String>());
              
              geneOGsMap.get(geneName).add(OGName);
              
              if(flag==1){
                  if(!OGGenesMap.containsKey(OGName))
                    OGGenesMap.put(OGName, new HashSet<String>());
                  OGGenesMap.get(OGName).add(geneName);
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
    
    
    public void printGeneOGMapping(){
        
        for(String s:geneOGsMap.keySet()){
            HashSet<String> tmp= geneOGsMap.get(s);
            
            System.out.print(s+": ");
                    
            for(String og:tmp){
                System.out.print(og+" ");
            }
            System.out.println("\n");
        }    
    }
    
    public void saveGeneOGMapping(String outFile){
        
        Iterator<String> keySetIterator = geneOGsMap.keySet().iterator();
        Path path= Paths.get(outFile);
        try{
        BufferedWriter writer = Files.newBufferedWriter(path,ENCODING);
        while(keySetIterator.hasNext()){
            String go=keySetIterator.next();
            HashSet<String> sTmp= geneOGsMap.get(go);
            
            writer.write(go+" ");
            
            for(String s:sTmp)
                writer.write(s+" ");;
                writer.write("\n");
        }
        writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
         
    }
