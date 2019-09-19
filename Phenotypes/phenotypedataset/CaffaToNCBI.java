/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import static phenotypedataset.LoadGeneData.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to map caffa PID to NCBI id
 */
public class CaffaToNCBI {
    
    static public void main(String [] args){
        
        File inputCaffa = new File("CaffaGene.txt");
        File inputSwissToUniprot = new File("SwissProtToUniProt.tab");
        File refUniprot = new File("gene_refseq_uniprotkb_collab.txt");
        
        
        HashMap<String, String> cafaIdToSwissProtMap = new HashMap<>();
        HashMap<String, String> inputSwissToUniprotMap = new HashMap<>();
        HashMap<String,String> UniProtToNCBIMap = new HashMap<>();
        
        File output = new File("CaffaToNCBI.txt");
       
        String extension[] = {"tfa"};
        BufferedReader reader;
       
         Path path =Paths.get(inputCaffa.getAbsolutePath());
        
         try{
        reader = Files.newBufferedReader(path,ENCODING);
                 String line="";
               while ((line = reader.readLine()) != null) {
                    
                   String tmp[] = line.split("\t");
                   cafaIdToSwissProtMap.put(tmp[0].trim(), tmp[1].trim());
          
                 }
               reader.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }
        
         path =Paths.get(inputSwissToUniprot.getAbsolutePath());
            try{
        reader = Files.newBufferedReader(path,ENCODING);
                 String line="";
                 int count=0;
               while ((line = reader.readLine()) != null) {
                    if(count==0){
                        count=1;
                        continue;
                    }
                   String tmp[] = line.split("\t");
                   inputSwissToUniprotMap.put(tmp[0].trim(), tmp[1].trim());
          
                 }
               reader.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }
            
            path =Paths.get(refUniprot.getAbsolutePath());
            try{
        reader = Files.newBufferedReader(path,ENCODING);
                 String line="";
                 int count=0;
               while ((line = reader.readLine()) != null) {
                    if(count==0){
                        count=1;
                        continue;
                    }
                   String tmp[] = line.split("\t");
                   UniProtToNCBIMap.put(tmp[1].trim(), tmp[0].trim());
          
                 }
               reader.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }  
            
            HashMap<String,String> caffaToNCBIMap = new HashMap<>();
            
            Iterator<String> it = cafaIdToSwissProtMap.keySet().iterator();
            
            while(it.hasNext()){
                String caffa = it.next();
                String swiss = cafaIdToSwissProtMap.get(caffa);
                String uniprot = inputSwissToUniprotMap.get(swiss);
                String ncbi =  UniProtToNCBIMap.get(uniprot);
                
                if(!inputSwissToUniprotMap.containsKey(swiss))
                    continue;
                if(!UniProtToNCBIMap.containsKey(uniprot))
                    continue;
               caffaToNCBIMap.put(caffa, ncbi);
            }
            
        
      try{
          FileWriter  fw = new FileWriter(output);
         
          String line = "";
          
           it = caffaToNCBIMap.keySet().iterator();
          
           while(it.hasNext()){
               String caffa = it.next();
          
               fw.write(caffa+"\t"+caffaToNCBIMap.get(caffa)+"\n");
           }
           
            fw.close();
                  } catch (Exception e) {
                         e.printStackTrace();
                   }         
        
    }
    
}
