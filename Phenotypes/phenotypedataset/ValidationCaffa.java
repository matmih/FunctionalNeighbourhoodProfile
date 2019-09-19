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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static jdk.nashorn.internal.objects.NativeArray.map;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *work performed while MM was visiting phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to output CaffaGenes
 */
public class ValidationCaffa {
    
    static public void main(String [] args){
        
        HashMap<String,HashSet<String>> geneCaffaFunction = new HashMap<>();
        HashMap<String,String> CaffaToSwiss = new HashMap<>();
        
        File output = new File("CaffaGene.txt");
        File input = new File("C:\\Users\\mmihelcic\\Downloads\\Supplementary_data\\data\\FilesForMappingAll");
       
        String extension[] = {"tfa"};
        BufferedReader reader;
         
      try{
         Collection files = FileUtils.listFiles(input, extension, true);

          FileWriter  fw = new FileWriter(output);
         
          String line = "";
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
               
                File inputF = (File) iterator.next();
                System.out.println("File = " + inputF.getAbsolutePath());
                 System.out.println("File = " + inputF.getName());
               Path path =Paths.get(inputF.getAbsolutePath());
               System.out.println("Path: "+inputF.getAbsolutePath());
               reader = Files.newBufferedReader(path,ENCODING);
                 
               while ((line = reader.readLine()) != null) {
                    if(line.contains(">")){
                        line = line.replace(">", "");
                        String tmp[] = line.split(" ");
                       fw.write(tmp[0]+"\t"+tmp[1]+"\n");
                    }
                    else continue;
          
                 }
               reader.close();
            }
            fw.close();
                  } catch (Exception e) {
                         e.printStackTrace();
                   }         
  }
}
