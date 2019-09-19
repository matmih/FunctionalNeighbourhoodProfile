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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create gene-COG mapping
 */
public class CreateGeneCOGMapping {
    static public void main(String [] args){
        
        String inputPath = "C:\\Users\\mmihelcic\\Desktop\\Sequence alignement and annotatoin\\SequencesClean\\SequencesAllC.txt";
        File input = new File(inputPath);
        String outputPath = "geneOgMapping.txt";
        File output = new File(outputPath);
        
        HashMap<String,HashSet<String>> mapping = new HashMap<>();
        
         BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
       String validPatternC =   "COG\\d{4}";
       Pattern patternC = Pattern.compile(validPatternC);
    Matcher matcherC = null; 
      while ((line = reader.readLine()) != null) {
          String tmp[] = line.split("\t");
          String gene = tmp[0].trim();
           matcherC=patternC.matcher(line);
           
           if(!mapping.containsKey(tmp[0].trim())){
                    mapping.put(tmp[0].trim(), new HashSet<>());
           }
           
           while(matcherC.find())
            {
                String cog = matcherC.group(0).trim();
                 mapping.get(gene).add(cog);
                  //System.out.println(matcherC.group(1));
                    }
                   
      }
         }
         catch(Exception e){
             e.printStackTrace();
         }
         
        
        try{ 
         
           FileWriter fw = new FileWriter(output); 
         
         Iterator<String> it = mapping.keySet().iterator();
         
         while(it.hasNext()){
             String og = it.next();
             
             String line = "";
             
             HashSet<String> cogs = mapping.get(og);
             
             if(cogs.size()==0)
                 continue;
             
             fw.write(og+":");
             
             for(String s:cogs)
                 line+=s+",";
             
             if(line.trim().equals(""))
                 continue;
             
             line = line.substring(0, line.length()-1);
             fw.write(line+"\n");
         }
         fw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
