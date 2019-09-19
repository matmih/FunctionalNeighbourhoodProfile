/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about LOR
 */
public class LoadOR {
    public HashMap<String, ArrayList<Double>> cont=new HashMap<>();
    
    public void loadContValue(File input){
         BufferedReader bufRdr=null;
           try 
        {
            bufRdr = new BufferedReader(new FileReader(input));
            String line;
            String cog="";
            int count=0;
            while ((line = bufRdr.readLine()) != null)
            {
                 count++;
                if(line.contains("Table: ")){
                    ArrayList<Double> tmpArr=new ArrayList<>();
                    String tmp[]=line.split(" ");
                    //System.out.println("tmp size: "+tmp.length);
                     cog=tmp[1].split("-")[1];
                     cont.put(cog, tmpArr);
                }
                if(count==3){
                    String tmp[]=line.split(" ");
                    int tp=Integer.parseInt(tmp[0]);
                    int fp=Integer.parseInt(tmp[1]);
                    cont.get(cog).add((double)tp);
                    cont.get(cog).add((double)fp);
                }
                if(count==4){
                    String tmp[]=line.split(" ");
                    int fn=Integer.parseInt(tmp[0]);
                    int tn=Integer.parseInt(tmp[1]);
                    cont.get(cog).add((double)fn);
                    cont.get(cog).add((double)tn);
                    count=-1;
                }     
        }
          }
           catch(Exception e){
               e.printStackTrace();
           }
    }
}
