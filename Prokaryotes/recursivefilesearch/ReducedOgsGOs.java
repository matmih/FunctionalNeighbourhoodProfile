/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import static recursivefilesearch.GOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to reduce the mapping to a selected subset of functions and attributes
 */
 
public class ReducedOgsGOs {
    public HashSet<String> redCogs=new HashSet<>();
    public HashSet<String> redGOs=new HashSet<>();
    
    public void LoadReducedOGGos(File inputCogs, File inputGOs){
               BufferedReader reader;
          try{
            Path path =Paths.get(inputCogs.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      while ((line = reader.readLine()) != null) {
          redCogs.add(line.trim());
           
        }
      reader.close();
      
      
      path =Paths.get(inputGOs.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
       line = null;
      
      while ((line = reader.readLine()) != null) {
          redGOs.add(line.trim());
        }
      reader.close();
      
      
        }
        catch(Exception e){
            e.printStackTrace();
        }
      }
    
}
