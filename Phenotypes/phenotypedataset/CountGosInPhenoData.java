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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description helper class to count the number of GOs in the mapping for phenotype prediction
 */
public class CountGosInPhenoData {
    public static void main(String [] args){
        
        File input = new File("C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\ReducedCogOGMappingPhenotype.txt");
     
        HashMap<String,Integer> goCount = new HashMap<>();
        HashMap<String,Integer> goCountRed = new HashMap<>();
        HashMap<String,HashSet<String>> cogGo = new HashMap<>();
        
         Path path =Paths.get(input.getAbsolutePath());
        
         try{
                 BufferedReader reader = Files.newBufferedReader(path,ENCODING);
                 
                 String line = "";
                 
                 while((line=reader.readLine())!=null){
                     
                     String tmp[] = line.split(":");
                     String gos[] = tmp[1].split(",");
                     
                     String cog = tmp[0].trim();
                     cogGo.put(cog, new HashSet<>());
                     
                     for(String g:gos){
                         cogGo.get(cog).add(g);
                         if(!goCount.containsKey(g))
                             goCount.put(g, 1);
                         else{
                             int count = goCount.get(g);
                             count++;
                             goCount.put(g, count);
                         }
                     }
                     
                 }
                 reader.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }
         
         Iterator<String> it = goCount.keySet().iterator();
         
         while(it.hasNext()){
             String g = it.next();
             
             if(goCount.get(g)>=5)
                 goCountRed.put(g, goCount.get(g));
         }
         
         System.out.println("Num GOs: "+goCount.keySet().size());
         System.out.println("Num red GOs: "+goCountRed.keySet().size());
         
         it = cogGo.keySet().iterator();
         
         while(it.hasNext()){
             String cog = it.next();
             
             HashSet<String> go = cogGo.get(cog);
             HashSet<String> toRemove = new HashSet<>();
             
             for(String g:go)
                 if(!goCountRed.containsKey(g))
                     toRemove.add(g);
             
             go.removeAll(toRemove);
             
         cogGo.put(cog, go);
         }
         
         
         it = cogGo.keySet().iterator();
        int countCog = 0;
         
          while(it.hasNext()){
             String cog = it.next();
             
             HashSet<String> go = cogGo.get(cog);
             
            if(go.size()>5)
                countCog++;
         }
          
          System.out.println("Num COGs reduced: "+countCog);
          
          File output = new File("C:\\Users\\mmihelcic\\Documents\\NetBeansProjects\\EColiStrainsDataset\\ReducedCogOGMappingPhenotypeR1.txt");
          
          
          it = cogGo.keySet().iterator();
        
          try{
          FileWriter  fw = new FileWriter(output);
          
          while(it.hasNext()){
              String cog = it.next();
              
              HashSet<String> gos = cogGo.get(cog);
              if(gos.size()<=5)
                  continue;
              
              fw.write(cog+":");
              
              String gosS = "";
              
              for(String g:gos)
                  gosS+=g+",";
              
              gosS = gosS.substring(0,gosS.length()-1);
               fw.write(gosS+"\n");
              
             }
          fw.close();
          }
          catch(IOException e){
              e.printStackTrace();
          }
                  
         
    }
}
