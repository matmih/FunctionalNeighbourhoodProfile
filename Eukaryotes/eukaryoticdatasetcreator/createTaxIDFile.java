/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create taxID file maping for eukaryotic organisms
 */
public class createTaxIDFile {
    public static void main(String [] args){
        
        //File output= new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\taxID.txt");
       // File inputFolder= new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiExtracted");
        
        File output= new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\taxID.txt");
        File inputFolder= new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\MetazoaExtracted");
        
         String extensions[] = {"dat"};
         boolean recursive = true;
         
         FileWriter fw = null;
         
          try{
                fw = new FileWriter(output);
            }
            catch(IOException e){
              e.printStackTrace();
            }
                 
        try{
         Collection files = FileUtils.listFiles(inputFolder, extensions, recursive);
         
         BufferedReader reader1;
          int fCount=0;

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                 System.out.println("File = " + input.getName());
                 System.out.println("Num processed: "+fCount);
                 
                     Path path =Paths.get(input.getAbsolutePath());
                     FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      
                      int taxonId=-1;
                       int geneSection=0;
                       int wrong=0; int locOK=0, locEx=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file

                          if(line.contains("db_xref=\"taxon:")){
                              String tn[]=line.split(":");
                              tn[1]=tn[1].replaceAll("\"", "");
                              taxonId=Integer.parseInt(tn[1]);
                              fCount++;
                              break;
                          }
                        }
                      fw.write(input.getName()+" "+taxonId+"\n");
                    }
            fw.close();
            }
          catch(IOException e){
                           e.printStackTrace();
                              }
    }
}
