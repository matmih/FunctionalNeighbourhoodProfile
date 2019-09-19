/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypedataset;

import phenotypedataset.Mappings;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import org.javatuples.*;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store OGs from the prokaryotic genome
 */
public class COG {//this function contains an arrayList of all the COGs occuring in the dataset, multiple occurances of the same COG will also be recorded
    final static Charset ENCODING = StandardCharsets.UTF_8;
    ArrayList<Triplet<Integer,Integer,String>> cogs=new ArrayList<Triplet<Integer,Integer,String>>();//all genes
    ArrayList<Triplet<Integer,Integer,String>> ancogs=new ArrayList<Triplet<Integer,Integer,String>>();//genes with corresponding COG groups
    int minCoordinate=0,maxCoordinate=0;//we save the maximum coordinate and take maxCoordinate=0 to compute min distance in both direction

    HashMap<String,Integer> cogToIndex=new HashMap<String,Integer>();//contains COG name COG index mapping

void createIndex(){
    int index=0;
    for(int i=0;i<ancogs.size();i++)
        if(!cogToIndex.containsKey(ancogs.get(i).getValue2()))
            cogToIndex.put(ancogs.get(i).getValue2(), index++);
}

    void findCogs(File input, int mode, Mappings map){ //mode=0, COG mapping from files, mode==1, separate ->gene mapping, one gene mapped to multiple OGs
        int lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          
          if(lineNum==1){
              String[] gensize=line.split(" ");
              String c=gensize[gensize.length-1];
              System.out.println("String c COG: "+c);
              String [] coord=c.split("\\.\\.");
              System.out.println("coord size: "+coord.length);
              minCoordinate=Integer.parseInt(coord[0]);
              maxCoordinate=Integer.parseInt(coord[1]);
          }
          
          if(lineNum>3){
          String[] st=line.split("\t");
          String[] coord=st[0].split("\\..");
          int x=Integer.parseInt(coord[0]), y=Integer.parseInt(coord[1]);

          String cogName="";
         if(mode==0){
          String pattern = "(^COG)(\\d+)(\\D+)";
       
          if(st[7].contains(","))
              st[7]=st[7].split(",")[0];
          cogName=st[7].replaceAll(pattern, "$1$2");
         }
         else{
             cogName=st[3];
         }

          Triplet<Integer,Integer,String> t=new Triplet<Integer,Integer,String>(x,y,cogName);
          cogs.add(t);

          if(map.geneOgMapping.containsKey(cogName))
              ancogs.add(t);
          }
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    void randomize(Random randGen){
       ArrayList<Triplet<Integer,Integer,String>> ShuffledCogs=new ArrayList<Triplet<Integer,Integer,String>>(); 
       HashSet<Integer> usedIndex=new HashSet<>();
       
       for(int i=0;i<cogs.size();i++){
           int ind=randGen.nextInt(cogs.size()-usedIndex.size());
           int count=0;
           
           for(int j=0;j<cogs.size();j++)
               if(!usedIndex.contains(j)){
                   if(count==ind){
                       usedIndex.add(j);
                       ShuffledCogs.add(cogs.get(j));
                       break;
                   }
                   count++;
               }           
       }  
       cogs=ShuffledCogs;
    }
    
    void randomizeLocation(Random randGen, int numShufflingSteps){
       ArrayList<Triplet<Integer,Integer,String>> ShuffledCogs=new ArrayList<Triplet<Integer,Integer,String>>(); 
       HashSet<Integer> usedIndex=new HashSet<>();
       
    for(int k=0;k<numShufflingSteps;k++){
       for(int i=0;i<cogs.size();i++){
           int ind=randGen.nextInt(cogs.size()-i);
           int count=0;
                       
                       int a, b;
                       String c;
                       
                       a=cogs.get(ind+i).getValue0();
                       b=cogs.get(ind+i).getValue1();
                       c=cogs.get(ind+i).getValue2();
                       
                       Triplet tmp=cogs.get(ind+i);
                       tmp=tmp.setAt0(cogs.get(i).getValue0());
                       tmp=tmp.setAt1(cogs.get(i).getValue1());
                       
                       Triplet tmp1=cogs.get(i);
                       tmp1=tmp1.setAt0(a);
                       tmp1=tmp1.setAt1(b);
                       
                       cogs.set(ind+i,tmp);
                       cogs.set(i, tmp1);
                       
                       ShuffledCogs.add(tmp);

                   count++;
       }  
    }
   }
}
