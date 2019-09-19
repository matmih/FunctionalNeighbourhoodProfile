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
import java.util.ArrayList;
import java.util.HashSet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to load phenotypes and strains from the corresponding files
 */
 
public class PhenotypesAndStrains {
    
    ArrayList<String> phenotypes = new ArrayList<>();
    ArrayList<String> strains = new ArrayList<>();
    
    
    public void loadPhenotypesAndStrains(String pathStrains, String pathPhenotypes){
        
         BufferedReader reader;
         try {
             
       File inputStrains = new File(pathStrains);    
      Path path =Paths.get(inputStrains.getAbsolutePath());
      System.out.println("Path: "+inputStrains.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null){
          
       strains.add(line.trim());
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
        
         
         try {
             
       File inputStrains = new File(pathPhenotypes);    
      Path path =Paths.get(inputStrains.getAbsolutePath());
      System.out.println("Path: "+inputStrains.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null){
          
       phenotypes.add(line.trim());
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            } 
         
    }
    
}
