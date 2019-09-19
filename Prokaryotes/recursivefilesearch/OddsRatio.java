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
import java.util.HashMap;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store ORs and read them from a file
 */
 
public class OddsRatio {
    HashMap<Integer,Double> OR=new HashMap<>();
    
    void readOR(File input){
         BufferedReader reader;
         
        try {
            Path path =Paths.get(input.getAbsolutePath());

            reader = Files.newBufferedReader(path,ENCODING);
            String line = null;
             while ((line = reader.readLine()) != null) {
                 String tmp[]=line.split(" ");
                     OR.put(Integer.parseInt(tmp[0]), Double.parseDouble(tmp[1]));
      }
             reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
}
