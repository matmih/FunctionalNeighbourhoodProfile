/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
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
  *@mail matmih1@gmail.com
  *@description class to read OOB predictions
 */
public class ReadOOBPredictions {
    
    public static void main(String [] args){
    File input = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\BaselineOGsk=4_oob.preds");
    HashMap<String, HashMap<String,Double>> data = new HashMap<>();
    
    Path p = Paths.get(input.getAbsolutePath());
    
    try{
        BufferedReader read = Files.newBufferedReader(p);
        
        String line="";
        while((line=read.readLine())!=null){
            String tmp[] = line.split("\t");
            HashMap<String,Double> goPreds = new HashMap<>();
            String cog = tmp[0].trim();
            String preds[] = tmp[2].trim().split(";");
            String gos[] = tmp[1].trim().split("@");
            int count=0;
            for(String gG:gos){
                String gosub[] = gG.trim().split("/");
                goPreds.put(gosub[gosub.length-1].trim(), Double.parseDouble(preds[count++]));
            }
            
            data.put(cog, goPreds);
            
        }
        read.close();
    }
    catch(IOException e){
    }
    
    try{
        FileWriter fw = new FileWriter("transformedPreds.pred");
        
        Iterator<String> it = data.keySet().iterator();
        int count = 0;
        ArrayList<String> goOrder = new ArrayList<>();
        
        while(it.hasNext()){
            String cog = it.next();
            
            HashMap<String,Double> preds = data.get(cog);
            
            Iterator<String> it1 = preds.keySet().iterator();
            ArrayList<Double> predsN = new ArrayList<>();
          
            if(count==0){
            while(it1.hasNext()){
                String go = it1.next();
                goOrder.add(go);
               predsN.add(preds.get(go));
                
            }
            }
            else{
                for(int i=0;i<goOrder.size();i++){
                    String go=goOrder.get(i);
                    predsN.add(preds.get(go));
                }
            }
            
            if(count==0)
            for(int i=0;i<goOrder.size();i++){
                fw.write(goOrder.get(i)+"\n");
            }
            
              fw.write(cog+" ");
             for(int i=0;i<predsN.size();i++){
                 if((i+1)<predsN.size())
                     fw.write(predsN.get(i)+" ");
                 else
                      fw.write(predsN.get(i)+"\n");
            }
            
            count=1;
            
        }
        fw.close();
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
   }
}
