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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to connect Caffa PIDs to OGs
 */
public class AddOgsToCaffaPIDs {
    public static void main(String [] args){
        HashMap<String,String> caffaToPGID = new HashMap<>();
        HashMap<String,String> PIDToCaffa = new HashMap<>();
        HashMap<String,Pair<String,ArrayList<String>>> caffaToPGIDOG = new HashMap<>();
        
        File input = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\caffaPID.txt");
        File inputpid = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
        File output = new File("caffaOGs");
        
        Path path = Paths.get(input.getAbsolutePath());
        BufferedReader read = null;
        
        try{
            read=Files.newBufferedReader(path, StandardCharsets.UTF_8);
        
        String line="";
        
        while((line = read.readLine())!=null){
            String tmp[] = line.split("\t");
            caffaToPGID.put(tmp[0].trim(),tmp[1].trim());
            PIDToCaffa.put(tmp[1].trim(), tmp[0].trim());
        }
        read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        path = Paths.get(inputpid.getAbsolutePath());
        try{
            read=Files.newBufferedReader(path, StandardCharsets.UTF_8);
        
        String line="";
        
        while((line = read.readLine())!=null){
            String tmp[] = line.split("\t");
            String pid = tmp[0].trim();
            String ogs[] = tmp[2].trim().split(",");
            
            if(!PIDToCaffa.containsKey(pid))
                continue;
            
            Pair<String,ArrayList<String>> p = new Pair(pid,new ArrayList<>());
            
            for(int i=0;i<ogs.length;i++){
                p.getValue1().add(ogs[i].trim());
            }
         
            String caffa = PIDToCaffa.get(pid);
            caffaToPGIDOG.put(caffa, p);
            
        }
        read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        try{
              FileWriter fw = new FileWriter(output);
              
              Iterator<String> it = caffaToPGIDOG.keySet().iterator();
              
              while(it.hasNext()){
                  String caffa = it.next();
                  Pair<String,ArrayList<String>> p = caffaToPGIDOG.get(caffa);
                  fw.write(caffa+"\t"+p.getValue0()+"\t");
                  for(int i=0;i<p.getValue1().size();i++){
                      if((i+1)<p.getValue1().size()){
                          fw.write(p.getValue1().get(i)+",");
                      }
                      else
                          fw.write(p.getValue1().get(i)+"\n"); 
                  }               
              }
              fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
    }
}
