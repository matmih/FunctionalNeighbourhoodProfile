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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import org.javatuples.Pair;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create data for comparative bar plots comparing performance between different predictive models
 */
public class CreateComparativeBarPlotMeasures {
    static public void main(String [] args){
        
        HashMap<String,ArrayList<Double>> result = new HashMap<>();
        
        String paths[]={"C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genome5LogK=1.out","C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genome5LogK=3.out","C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\genome5LogK=10.out","C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\GFPResultsLab.txt","C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\BaselineOGsk=4_600_.oob"};
        int fileType[]={0,0,0,1,2};
        File input;
        int useAUC=0;
        
        String methods[] = {"1-NN","3-NN","10-NN","GFP","GFN"};
        
        ArrayList<String> functions = new ArrayList<>();
        
        functions.add("GO0005975"); functions.add("GO0006260"); functions.add("GO0006974");
        functions.add("GO0016051"); functions.add("GO0006457"); functions.add("GO0046700");
        functions.add("GO0006310"); functions.add("GO0046903"); functions.add("GO0008610");
        functions.add("GO0006865"); functions.add("GO0006508");
        
        for(int i=0;i<paths.length;i++){
            input = new File(paths[i]);
            
            if(fileType[i] == 0){
              BufferedReader reader;  
                 Path path =Paths.get(input.getAbsolutePath());
            String NN = "";
              
              if(i==0)
                  NN="1-NN";
              else if(i==1)
                  NN="3-NN";
              else if(i==2)
                  NN="10-NN";
              
              result.put(NN, new ArrayList<Double>());
     try{
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int summary=1, dataSection=0, test=0;
      ArrayList<Double> tmpRes = new ArrayList();
      for(int k=0;k<functions.size();k++)
          tmpRes.add(0.0);
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Summary")){
              summary=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(summary==1){
              if(line.contains("Testing error")){
                  test=1;
                  summary=0;
                  continue;
              }
          }
          
          if(test==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  test=0;
                  continue;
              }
          }
          
          if(dataSection==1){              
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              if(!functions.contains(first.trim()))
                  continue;
              System.out.println("firstBase: "+first);
              String AUPRC=tmp[2];
              String AUC = tmp[1];
              AUPRC=tmp[2].split(": ")[1];
              AUC=tmp[1].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              int index=-1;
              
              for(int k=0;k<functions.size();k++){
                  if(functions.get(k).equals(first.trim())){
                      index = k;
                      break;
                  }              
              }
              
               if(useAUC==0)
                   tmpRes.set(index, Double.parseDouble(AUPRC));
              else
                   tmpRes.set(index, Double.parseDouble(AUC));
          }
 
      }
       for(int k=0;k<tmpRes.size();k++){
              result.get(NN).add(tmpRes.get(k));
          }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }               
        }
            else if(fileType[i] == 1){
                 BufferedReader reader;
                 try{
                     result.put("GFP", new ArrayList<Double>());
                     ArrayList<Double> tmpRes = new ArrayList();
                     for(int k=0;k<functions.size();k++)
                             tmpRes.add(0.0);
                        Path path =Paths.get(input.getAbsolutePath());

                        reader = Files.newBufferedReader(path,ENCODING);
                        String line = null;
                         while ((line = reader.readLine()) != null) {
                             String s[] = line.split("\t");
                             String first = s[0].trim();
                             if(!functions.contains(first))
                                 continue;
                             
                              int index=-1;
              
                         for(int k=0;k<functions.size();k++){
                                 if(functions.get(k).equals(first.trim())){
                                        index = k;
                                        break;
                                    }              
                            }
                             
                             if(useAUC==0)
                                 tmpRes.set(index, Double.parseDouble(s[2].trim()));
                             else
                                 tmpRes.set(index, Double.parseDouble(s[1].trim()));
                         }
                         
                         for(int k=0;k<tmpRes.size();k++){
                                result.get("GFP").add(tmpRes.get(k));
                           }
                   reader.close();     
                 }
                 catch(IOException e){
                     e.printStackTrace();
                 }
            }
            else if(fileType[i]==2){
                 BufferedReader reader;
                  ArrayList<Double> tmpRes = new ArrayList();
                     for(int k=0;k<functions.size();k++)
                             tmpRes.add(0.0);
         try {
             
             result.put("GFN", new ArrayList<Double>());
             
      Path path =Paths.get(input.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              if(!functions.contains(first.trim()))
                  continue;
              
              int index=-1;
              
                         for(int k=0;k<functions.size();k++){
                                 if(functions.get(k).equals(first.trim())){
                                        index = k;
                                        break;
                                    }              
                            }
              
              String AUPRC=tmp[2];
              String AUC = tmp[1];
              AUPRC=tmp[2].split(": ")[1];
               AUC=tmp[1].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              if(useAUC==0)
                   tmpRes.set(index, Double.parseDouble(AUPRC));
              else
                   tmpRes.set(index, Double.parseDouble(AUC));
          }        
      }
        for(int k=0;k<tmpRes.size();k++){
                                result.get("GFN").add(tmpRes.get(k));
        }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
            }
            
        }

           String out="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\BarPlotComparativeData.txt";
          
           if(useAUC==1)
          out="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\BarPlotComparativeDataAUC.txt";
          
         try{
         FileWriter fw = new FileWriter(out);
         
         for(int i=0;i<methods.length;i++){
             String method = methods[i];
             
             ArrayList<Double> res = result.get(method);
             
             for(int j=0;j<res.size();j++)
                 fw.write(res.get(j)+" ");
             
             fw.write("\n");
             
         }

         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}         
     }   
}
