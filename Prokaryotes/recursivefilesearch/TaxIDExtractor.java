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
import static recursivefilesearch.COG.ENCODING;

/**
 *
 * @author matej
 */
public class TaxIDExtractor {
    String TaxID="";
    
    void ExtractID(File input){
        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          if(line.contains("/db_xref=\"taxon:")){
              System.out.println(line);
              String tmp[]=line.split(":");
              TaxID=new String(tmp[1].substring(0,tmp[1].length()-1));
          }
      }
       reader.close();
         }
            catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
      }
    }