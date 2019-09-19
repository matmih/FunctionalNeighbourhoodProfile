/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to convert GO number to a GO short name
 */
public class AddNamesToHierarchyAttribute {
    public static void main(String args[]){
        
        File input= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\header.txt");
        
        String outputS="";
        
        try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(input)))
        {
            String line;
            int count=0;
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                String tmp[]=line.split(",");
                
                for(int i=0;i<tmp.length;i++){
                    String pair[]=tmp[i].split("\\/");
                    
                    System.out.println(pair[0]+" "+pair[1]);
                  
                 for(int j=0;j<2;j++){
                    String GOName="";
                  
                    if(pair[j].equals("root")){
                       if(j==0)
                        outputS+="root"+"/";
                       else
                           outputS+="root";
                        continue;
                    }
                        
                    
                    GOName+="GO";
                    int ogNumber=Integer.parseInt(pair[j]);
                    int numTmpGo=ogNumber, nDigGo=0;
                    
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=ogNumber;
                   
                    if(j==0)
                        outputS+=(GOName+"/");
                    else outputS+=GOName;
              
                }
                 if(i<(tmp.length-1))
                 outputS+=",";
              }  
                
            }
            bufRdr1.close();
        }
       catch(Exception e){
           e.printStackTrace();
       }
         
         File output=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\headerGONew.txt");
         
         try
            {
                FileWriter fw = new FileWriter(output);
                
                
                    fw.write(outputS);
                
                   
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            }  
    }
}
