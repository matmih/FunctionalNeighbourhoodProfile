/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to extract Prokaryotes that are classificed in a group 'prokaryotic bacteria'
 */
public class ExtractProkoarioticBact {
   public static void main(String args[]){//go_201401-termdb.obo-xml
       
       File oboFile=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\go_20141215-termdb.obo-xml\\go_daily-termdb.obo-xml");
       HashSet<String> procGO=new HashSet<>();
       
       try (BufferedReader bufRdr = new BufferedReader(new FileReader(oboFile)))
        {
            String line;

            String go = null;
            int termGOFound=0;
            
            
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                
                if(line.startsWith("<term>"))
                    termGOFound=1;
                
                if (line.startsWith("<id>GO:") && termGOFound==1){
                    go=line.replaceAll(":", "");
                    go=go.replaceAll("<id>", "");
                    go=go.replace("</id>", "");
                    System.out.println(go);
                }
                 
                if(line.startsWith("<subset>") && termGOFound==1 && line.contains("gosubset_prok")){
                    procGO.add(go);
                }
                
                if(line.startsWith("</term>"))
                    termGOFound=0;
                            
            }
        }
       catch(Exception e){
           e.printStackTrace();
       }
       
       String output="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\prokaryotic.txt";
       FileWriter fwGnu = null;
       try{
            fwGnu=new FileWriter(output);
       
      
       BufferedWriter fw = new BufferedWriter(fwGnu);
       
       for(String s:procGO)
           fw.write(s+"\n");
       }
       catch(Exception e){
           e.printStackTrace();
       }      
    }
}
