/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
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
  *@mail matmih1@gmail.com
  *@description class to create OGGOMappings
 */
 
public class OGGOMapping {
    /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

     final static Charset ENCODING = StandardCharsets.UTF_8;
    public HashMap<String,ArrayList<String>> CogGOmap=new HashMap<String,ArrayList<String>>();
    public HashMap<String,Integer> GOtoIndex=new HashMap<String,Integer>();
    public HashMap<Integer,String> IndexToGO=new HashMap<>();
    
    public void createCOGGOMapping(File input){//implement new func
        int lineNum=0,index=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          if(lineNum>1){
          String[] st=line.split(" ");
          
          String ogName=st[0];
          ogName=ogName.replace(":", "");
          String OGName=ogName;
          
           if(!CogGOmap.containsKey(OGName)){
                ArrayList<String> Gofunc=new ArrayList<String>();
                CogGOmap.put(OGName, Gofunc);
           }
           
         for(int i=1;i<st.length;i++){
             String go=st[i].replace(":", ""); 
             if(go.equals(" "))
                 continue;
             if(!GOtoIndex.containsKey(go)){
                    GOtoIndex.put(go, index++);
                    IndexToGO.put(index-1, go);
              }
             
              CogGOmap.get(OGName).add(go);
         }
        }
      }
            reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         } 
}
