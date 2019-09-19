/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to create Eukaryotic datasets
 */
public class DatasetWriter {
    
          void createAssociationGraphHeader(String trainSet,OGGOMapping cgmap,LocationSimilarity features ,int isTrain){
           try
{
    FileWriter fw = new FileWriter(trainSet,true); //the true will append the new data

     Iterator<String> keySetIterator = features.locationN.keySet().iterator();//cgmap.CogGOmap.keySet().iterator();

     fw.write("COG\t");

     
      for(int i=0;i</*features.IndexCOG.keySet().size()*/features.fullSize;i++){
         if(features.IndexCOG.get(i) == null || features.IndexCOG.get(i).equals("null"))
             continue;
         fw.write("DT"+features.IndexCOG.get(i)+"\t");
     }
     
    fw.write("\n");
    fw.close();
}
catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
       }
          
          
    void appendWeightsToGraph(String trainSet,OGGOMapping map,LocationSimilarity features, int isTrain){
        
         FileWriter fw;
                try
                {
           fw = new FileWriter(trainSet,true); //the true will append the new data

           Iterator<String> keySetIterator = features.locationN.keySet().iterator();
           String row="";
           
           System.out.println("Number of COGs in append: "+features.locationN.keySet().size());
          while(keySetIterator.hasNext()){
              String COG=keySetIterator.next();
              
              System.out.println("Writing COG: "+COG);
             ArrayList<String> cgm=map.CogGOmap.get(COG);
              int count=0;
              if(map.CogGOmap.containsKey(COG))
                  count=cgm.size();
              else continue;
              
              row+=COG+"\t";
              Pair<ArrayList<Double>,ArrayList<Double>> pt=features.locationN.get(COG);
              ArrayList<Double> tmp=pt.getValue0();
              for(int i=0;i<tmp.size();i++)
                  if(i+1<tmp.size())
                  row+=tmp.get(i)+"\t";
                  else row+=tmp.get(i)+"\n";

              fw.write(row);
              row="";
             }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }       
    
    void createLocationTrainHeader(String trainSet, OGGOMapping cgmap, LocationSimilarity features, File headerInput, int isTrain){
        try
    {
    FileWriter fw = new FileWriter(trainSet,true); //the true will append the new data

     fw.write("@RELATION LocationFeatures\n\n");
     fw.write("@ATTRIBUTE ID         string\n");
     
     for(int i=0;i</*features.IndexCOG.keySet().size()*/features.fullSize;i++){
         if(features.IndexCOG.get(i) == null || features.IndexCOG.get(i).equals("null"))
             continue;
         fw.write("@ATTRIBUTE "+"DT"+features.IndexCOG.get(i)+"     numeric\n");
     }

     String hier= "";
     BufferedReader reader;
     
       try{
                     Path path =Paths.get(headerInput.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            hier+=line;
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 
       
     if(isTrain==1)
         fw.write("@ATTRIBUTE class hierarchical "+hier+"\n\n");
    fw.write("@DATA\n");//appends the string to the file
    
    fw.close();
}
catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }
    
    void appendLocationRowsToSet(String trainSet,OGGOMapping map,LocationSimilarity features, int isTrain){
        
         FileWriter fw;
                try
                {
           fw = new FileWriter(trainSet,true); //the true will append the new data

           Iterator<String> keySetIterator = features.locationN.keySet().iterator();
           String row="";
           StringBuilder sb = new StringBuilder();
           System.out.println("Number of COGs in append: "+features.locationN.keySet().size());
          
           int maxInd=-1;
           
            for(int v:features.IndexCOG.keySet())
                if(v>maxInd)
                    maxInd=v;
                    
           
           for(int i1=0;i1<maxInd+1/*features.IndexCOG.keySet().size()*/;i1++){
        // fw.write("@ATTRIBUTE "+"DT"+features.IndexCOG.get(i)+"     numeric\n");
                
               if(!features.IndexCOG.containsKey(i1))
                   continue;
               
               System.out.println("COGs feat: "+features.IndexCOG.get(i1));
                System.out.println("Index feat: "+i1);
               System.out.println("Diag: "+features.locationN.get(features.IndexCOG.get(i1)).getValue0().get(i1));
               
               int test=1;//remove to write to file
               /*if(test==1)
               continue;*/
          // while(keySetIterator.hasNext()){
              String COG=features.IndexCOG.get(i1);//keySetIterator.next();
              if(COG==null)
                  continue;
              
              System.out.println("Writing COG: "+COG);
             ArrayList<String> cgm=map.CogGOmap.get(COG);
              int count=0;
              if(map.CogGOmap.containsKey(COG))
                  count=cgm.size();
              else continue;
              
              //row+="\""+COG+"\", ";
              sb.append(("\""+COG+"\", "));
              Pair<ArrayList<Double>,ArrayList<Double>> pt=features.locationN.get(COG);
              ArrayList<Double> tmp=pt.getValue0();
              
              int countNonMissing=0;
              
              for(int i=0;i<tmp.size();i++){
                  if(Math.abs(tmp.get(i))<=Double.MAX_VALUE)//mora biti infinity u datoteci
                      countNonMissing++;
                  //if(Math.abs(tmp.get(i))<=Double.MAX_VALUE)
                    //row+=tmp.get(i)+", ";
                    sb.append((tmp.get(i)+", "));
                 // else
                   //   row+="?, ";
              }
              
              System.out.println("cnm: "+countNonMissing);

                if(countNonMissing==1){
                    System.out.println("OG: "+COG+" removed...");
                    //row="";
                    sb.delete(0, sb.length());
                    continue;
                }

              if(count>0 && isTrain==1){
                  StringBuilder target = new StringBuilder();
                 // String target="";
                  for (int j=0;j<cgm.size();j++)
                      if(j<cgm.size()-1)
                        // target+=cgm.get(j)+"@";
                          target.append((cgm.get(j)+"@"));
                      else
                          target.append(cgm.get(j));
                          //target+=cgm.get(j);
                // row+=target;
                 sb.append(target.toString());
                 sb.append("\n");
                 //row+="\n";
                 //fw.write(row);
                 fw.write(sb.toString());
              }
              else if(count==0 || isTrain==0){
                  System.out.println("COG: "+COG+" has no function, error!");
                 // row="";
                  sb.delete(0, sb.length());
                  continue;
                  /*row=row.substring(0, row.length()-2);
                  row+="\n";
                  fw.write(row);*/
              }
              //row="";
              sb.delete(0, sb.length());
             }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }

    
      
    void createTrainHeader(String trainSet, OGGOMapping cgmap, File headerInput, int isTrain){
        try
{
    FileWriter fw = new FileWriter(trainSet,true); //the true will append the new data

     fw.write("@RELATION FungyGenomes\n\n");
     fw.write("@ATTRIBUTE ID	string\n");
     for(int i=0;i<cgmap.GOtoIndex.keySet().size();i++){
          Iterator<String> keySetIterator = cgmap.GOtoIndex.keySet().iterator();
                while(keySetIterator.hasNext()){
                    String go=keySetIterator.next();
                    int index=cgmap.GOtoIndex.get(go);
                    if(index==i){
                        fw.write("@ATTRIBUTE "+go+"			numeric\n");
                        break;
                    }
                   }
     }
     
     String hier= "";
     BufferedReader reader;
     
       try{
                     Path path =Paths.get(headerInput.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            hier+=line;
                    }
                reader.close();
                }
                catch(Exception e){e.printStackTrace();} 

     if(isTrain==1)
         fw.write("@ATTRIBUTE class hierarchical "+hier+"\n\n");
    fw.write("@DATA\n");//appends the string to the file
    fw.close();
}
catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }
    
          void appendBaselineRowsToSet(String trainSet,OGGOMapping map,OGFeatures features, int isTrain){
        
         FileWriter fw;
                try
                {
           fw = new FileWriter(trainSet,true); //the true will append the new data

           Iterator<String> keySetIterator = features.COGbaselinefeaturesmap.keySet().iterator();
           String row="";
          while(keySetIterator.hasNext()){
              row="";
              String COG=keySetIterator.next();
              
              
             ArrayList<String> cgm=map.CogGOmap.get(COG);
              int count=0;
              if(map.CogGOmap.containsKey(COG))
                  count=cgm.size();
              else continue;
              
              if(count==0)
                  continue;
              
              row+="\""+COG+"\", ";
              Pair<ArrayList<Double>,Double> pt=features.COGbaselinefeaturesmap.get(COG);
              ArrayList<Double> tmp=pt.getValue0();
              
              int cnz=0; 
              for(int cnzK = 0;cnzK<tmp.size();cnzK++){
                  if(tmp.get(cnzK)!=0.0)
                      cnz++;
              }
              
              if(tmp.size()==0 || cnz==0)
                  continue;
              for(int i=0;i<tmp.size();i++)
                  row+=tmp.get(i)+", ";


              if(count>0 && isTrain==1){
                  String target="";
                  for (int j=0;j<cgm.size();j++)
                      if(j<cgm.size()-1)
                         target+=cgm.get(j)+"@";
                      else
                          target+=cgm.get(j);
                 row+=target;
                 row+="\n";
                 fw.write(row);
              }
              else if(count==0 && isTrain==0){
                  row=row.substring(0, row.length()-2);
                  row+="\n";
                  fw.write(row);
              }
              row="";
             }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }  
          
    void appendBaselineRowsToSetTrain(String trainSet,OGGOMapping map,OGFeatures features, HashSet<String> leaveoutOgs , int isTrain){
        
         FileWriter fw;
                try
                {
           fw = new FileWriter(trainSet,true); //the true will append the new data

           Iterator<String> keySetIterator = features.COGbaselinefeaturesmap.keySet().iterator();
           String row="";
          while(keySetIterator.hasNext()){
              row="";
              String COG=keySetIterator.next();
              
              if(leaveoutOgs.contains(COG))
                  continue;
              
             ArrayList<String> cgm=map.CogGOmap.get(COG);
              int count=0;
              if(map.CogGOmap.containsKey(COG))
                  count=cgm.size();
              else continue;
              
              if(count==0)
                  continue;
              
              row+="\""+COG+"\", ";
              Pair<ArrayList<Double>,Double> pt=features.COGbaselinefeaturesmap.get(COG);
              ArrayList<Double> tmp=pt.getValue0();
              
              int cnz=0; 
              for(int cnzK = 0;cnzK<tmp.size();cnzK++){
                  if(tmp.get(cnzK)!=0.0)
                      cnz++;
              }
              
              if(tmp.size()==0 || cnz==0)
                  continue;
              for(int i=0;i<tmp.size();i++)
                  row+=tmp.get(i)+", ";


              if(count>0 && isTrain==1){
                  String target="";
                  for (int j=0;j<cgm.size();j++)
                      if(j<cgm.size()-1)
                         target+=cgm.get(j)+"@";
                      else
                          target+=cgm.get(j);
                 row+=target;
                 row+="\n";
                 fw.write(row);
              }
              else if(count==0 && isTrain==0){
                  row=row.substring(0, row.length()-2);
                  row+="\n";
                  fw.write(row);
              }
              row="";
             }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    }  
    
    
    void appendBaselineRowsToSetTest(String trainSet,OGGOMapping map,OGFeatures features, HashSet<String> leaveoutOgs ,int isTrain){
        
         FileWriter fw;
                try
                {
           fw = new FileWriter(trainSet,true); //the true will append the new data

           Iterator<String> keySetIterator = features.COGbaselinefeaturesmap.keySet().iterator();
           String row="";
          while(keySetIterator.hasNext()){
              row="";
              String COG=keySetIterator.next();
              
              if(!leaveoutOgs.contains(COG))
                  continue;
              
             ArrayList<String> cgm=map.CogGOmap.get(COG);
              int count=0;
              if(map.CogGOmap.containsKey(COG))
                  count=cgm.size();
              else continue;
              
              if(count==0)
                  continue;
              
              row+="\""+COG+"\", ";
              Pair<ArrayList<Double>,Double> pt=features.COGbaselinefeaturesmap.get(COG);
              ArrayList<Double> tmp=pt.getValue0();
              
              int cnz=0; 
              for(int cnzK = 0;cnzK<tmp.size();cnzK++){
                  if(tmp.get(cnzK)!=0.0)
                      cnz++;
              }
              
              if(tmp.size()==0 || cnz==0)
                  continue;
              for(int i=0;i<tmp.size();i++)
                  row+=tmp.get(i)+", ";


              if(count>0 && isTrain==1){
                  String target="";
                  for (int j=0;j<cgm.size();j++)
                      if(j<cgm.size()-1)
                         target+=cgm.get(j)+"@";
                      else
                          target+=cgm.get(j);
                 row+=target;
                 row+="\n";
                 fw.write(row);
              }
              else if(count==0 && isTrain==0){
                  row=row.substring(0, row.length()-2);
                  row+="\n";
                  fw.write(row);
              }
              row="";
             }
           fw.close();
                }
          catch(IOException ioe)
{
    System.err.println("IOException: " + ioe.getMessage());
}
    } 
          
    
}