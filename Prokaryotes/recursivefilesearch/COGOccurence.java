/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
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
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class used to compute different statistics about COG occurrence
 */
 
public class COGOccurence {
    HashMap<String,Integer> occurence=new HashMap<>();
    
    void initialize(COGGOMap cgmap){
        for(String CGName:cgmap.CogGOmap.keySet()){
            occurence.put(CGName, 0);
        }
    }
    
    String findMinimumOccuring(HashSet<String> ogs){
        String tmp="";
        int minCount=Integer.MAX_VALUE;
        
        for(String s:ogs){
            int c = occurence.get(s);
            if(c<minCount){
                minCount = c;
                tmp = s;
            }
        }
     
        return tmp;
    }
    
    void createOccurence(COG cog, COGGOMap cgmap){
    
       for(String CGName:cgmap.CogGOmap.keySet()){
        for(Triplet cogTmp:cog.ancogs){
            String cname=cogTmp.getValue2().toString();
            if(CGName.contentEquals(cname)){
                int val=occurence.get(CGName);
                val++;
                occurence.put(CGName, val);
                break;
            }
                
        }
       }
}
    
    
    void createOccurence1(COG cog, GeneOGMapping geneOGMap){
        
        String gene="";
        
        for(int i=0;i<cog.ancogs.size();i++){
            gene = cog.ancogs.get(i).getValue2();
            HashSet<String> ogs = geneOGMap.geneOGsMap.get(gene);
            
            for(String s:ogs){
                if(!occurence.containsKey(s)){
                    occurence.put(s, 1);
                }
                else{
                    int c = occurence.get(s);
                    c++;
                    occurence.put(s, c);
                }
            }
            
        }     
    }
    
    
    void computeCoocurence(File folderPath, int taxIdT, String[] extensions, Boolean recursive, String taxIDFilePath, GeneOGMapping geneOGMap){
           try{
                Collection files = FileUtils.listFiles(folderPath, extensions, recursive);            
               
                BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
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
                 
                  COG cg=new COG();
                  cg.findCogs(input,1,geneOGMap); //change to 0 to get COGs from file
                  createOccurence1(cg,geneOGMap);
                }
                
           }
           catch(Exception e){
               e.printStackTrace();
           }
    }
    
    void clearMapping(){
        occurence.clear();
    }
    
    void writeToFile(String output){
         FileWriter fw;
                try
                {
           fw = new FileWriter(output); 

           Iterator<String> keySetIterator = occurence.keySet().iterator();
         
          while(keySetIterator.hasNext()){
              String COG=keySetIterator.next();
              int value=occurence.get(COG);
              
              fw.write(COG+" "+value+"\n");
              
          }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }

}