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
import java.util.HashSet;
import static recursivefilesearch.GOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing helper functions to easily translate between the string and the numeric representation of a GO function
 */
 
public class ReducedGOTranslator {
    
    public HashSet<String> goS=new HashSet<>();
    public HashSet<Integer> goI=new HashSet<>();
    
    public void ReadAndTranslate(File input){
        
        BufferedReader reader;
          try{
            Path path =Paths.get(input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      
      while ((line = reader.readLine()) != null) {
          goS.add(line.trim());
           
        }
      reader.close();
    }
          catch(IOException e){
              e.printStackTrace();
          }
          
          for(String s:goS){
              String tmp=s.replace("GO","");
              int numb=Integer.parseInt(tmp);
              System.out.println("GO: "+numb);
              goI.add(numb);
          }
          
    }
    
    public void translate(HashSet<Integer> in){
        goS.clear(); goI=in;
        int numTmpGo, nDigGo;
        for(int go:in){
            String GOName="GO";
                 int nNumGo=go;
             numTmpGo=nNumGo;  nDigGo=0;
              
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=nNumGo;
              goS.add(GOName);
        }
    }
    
   HashSet<Integer> translateToInt(HashSet<String> in){ 
        HashSet<Integer> retGos=new HashSet<>();
        
        for(String go:in){
             String tmp=go.replace("GO","");
              int numb=Integer.parseInt(tmp);
              retGos.add(numb);
        }
        
        return retGos;
   }
    
   HashSet<String> translateToString(HashSet<Integer> in){
       HashSet<String> retGos=new HashSet<>();
       
        int numTmpGo, nDigGo;
        for(int go:in){
            String GOName="GO";
                 int nNumGo=go;
             numTmpGo=nNumGo;  nDigGo=0;
              
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=nNumGo;
              retGos.add(GOName);
        }
        
        return retGos;
    }
    
   
   
    int translateToIntSingle(String in){ 
        

             String tmp=in.replace("GO","");
              int numb=Integer.parseInt(tmp);
        
        return numb;
   }
    
   String translateToStringSingle(Integer in){
       
        int numTmpGo, nDigGo;
            String GOName="GO";
                 int nNumGo=in;
             numTmpGo=nNumGo;  nDigGo=0;
              
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=nNumGo;
        
        return GOName;
    }
    
}
