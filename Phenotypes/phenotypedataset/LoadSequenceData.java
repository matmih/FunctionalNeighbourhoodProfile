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
import java.util.HashMap;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to load bacterial DNA sequence.
 */
public class LoadSequenceData {
    public HashMap<String,String> contigSequence;
    
    public LoadSequenceData(){
        contigSequence = new HashMap<>();
    }
    
    public void loadData(File sequence){
        
        
        BufferedReader reader;
        StringBuilder b = new StringBuilder();
         try {
      Path path =Paths.get(sequence.getAbsolutePath());
      System.out.println("Path: "+sequence.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      ArrayList<String> contigs = new ArrayList<>();
      while ((line = reader.readLine()) != null) {
        
          if(line.contains(">")){
              if(contigs.size()>0)
                 contigSequence.put(contigs.get(contigs.size()-1), b.toString());
              b = new StringBuilder();
              String tmp[] = line.split(" ");
              String contid = tmp[0].replace(">", "").trim();
              contigSequence.put(contid, "");
              contigs.add(contid);
          }
          else if(!line.equals("")){
             b.append(line.trim());
          } 
          else if(line.equals(""))
              contigSequence.put(contigs.get(contigs.size()-1), b.toString());
      }
      contigSequence.put(contigs.get(contigs.size()-1), b.toString());
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }      
        
    }
    
    public void shortOutput(){
        Iterator<String> it = contigSequence.keySet().iterator();
        
        while(it.hasNext()){
            String cont = it.next();
            System.out.println(cont+" : "+contigSequence.get(cont).substring(0,10));
        }
        
        it = contigSequence.keySet().iterator();
        String f = it.next();
        System.out.println("Short query: ");
        if(f.length()>400)
        System.out.println(contigSequence.get(f).substring(40, 400));
        else System.out.println(f);
        
    }
    
    
}
