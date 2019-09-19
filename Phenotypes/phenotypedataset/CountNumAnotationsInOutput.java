/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import static phenotypedataset.LoadGeneData.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description helper class to count number of genes assigned to COGs/KEggs and GOs
 */
public class CountNumAnotationsInOutput {
    public static void main(String[] args) {
        
        File input = new File("C:\\Users\\mmihelcic\\Desktop\\Sequence alignement and annotatoin\\SequencesClean\\SequencesAll.txt");
        
        int genesWithCogs=0;
        int genesWithGO = 0;
        int genesWithKegg = 0;
        int count=0;
        
         BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          
          if(count<4){
              count++;
              continue;
          }
          
          String tmp[] = line.split("\t");
          
          if(tmp.length>1){
              
              String validPatternC =   "COG\\d{4}";
    Pattern patternC = Pattern.compile(validPatternC);
    Matcher matcherC = null; 
              
              for(int i=0;i<tmp.length;i++){
                  matcherC=patternC.matcher(tmp[i]);  

                  if(matcherC.find() == true){
                      genesWithCogs++;
                      break;
                  }
              }
              
               for(int i=0;i<tmp.length;i++){
                  if(tmp[i].contains("GO:")){
                      genesWithGO++;
                      break;
                  }
              }
               
   String validPattern =   "K\\d+";
    Pattern pattern = Pattern.compile(validPattern);
    Matcher matcher = null; 
    
               for(int i=0;i<tmp.length;i++){
                   matcher=pattern.matcher(tmp[i]);  
                  if(matcher.find() == true){
                      genesWithKegg++;
                      break;
                  }
              }    
          }   
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
         System.out.println("genesWithCog: "+genesWithCogs);
         System.out.println("genesWithGO: "+genesWithGO);
         System.out.println("genesWithKegg: "+genesWithKegg);
         
    }
}
