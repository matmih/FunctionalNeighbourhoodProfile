/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;

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
  *@description class to create COG GO mapping from the given arff file
 */
 
public class COGGOMap {
     final static Charset ENCODING = StandardCharsets.UTF_8;
    public HashMap<String,ArrayList<String>> CogGOmap=new HashMap<String,ArrayList<String>>();
    public HashMap<String,Integer> GOtoIndex=new HashMap<String,Integer>();
    public HashMap<Integer,String> IndexToGO=new HashMap<>();

    public void reduceMap(ReducedOgsGOs red){
        
        Iterator<String> it=CogGOmap.keySet().iterator();
        
        while(it.hasNext()){
            String COG=it.next();
           
            if(!red.redCogs.contains(COG))
                it.remove();
      }
        
        it=CogGOmap.keySet().iterator();
        int countRoots=0;
        while(it.hasNext()){
            countRoots=0;
            String COG=it.next();
           
            
            ArrayList<String> GOs=CogGOmap.get(COG);
            
            for(int k=GOs.size()-1;k>=0;k--){
                if(!red.redGOs.contains(GOs.get(k))){
                    GOs.remove(k);
                     continue;
                }
                                
                if(GOs.get(k).equals("GO0008150") || GOs.get(k).equals("GO0003674") || GOs.get(k).equals("GO0005575"))
                        countRoots++;
                
            }
                
                if(countRoots<GOs.size())
                CogGOmap.put(COG, GOs);
        }
        
        
         it=CogGOmap.keySet().iterator();
         
          while(it.hasNext()){
            String COG=it.next();
           
            ArrayList<String> GOs=CogGOmap.get(COG);
            
           if(GOs.size()<1)
               it.remove();
        }
        
          Iterator<String> it1=GOtoIndex.keySet().iterator();
          
          while(it1.hasNext()){
              String GO=it1.next();
              
              if(!red.redGOs.contains(GO)){
                  int goID=GOtoIndex.get(GO);
                  IndexToGO.remove(goID);
              }
          }
          
          int count=0;
          HashMap<Integer,String> tmpMap=new HashMap<>();
          HashSet<Integer> usedIndex=new HashSet<>();
          
          for(int i=0;i<GOtoIndex.keySet().size();i++){
              if(IndexToGO.containsKey(i)){
                   tmpMap.put(count++, IndexToGO.get(i));
              }
          }
          
          IndexToGO=tmpMap;

          GOtoIndex.clear();
          
          Iterator<Integer> it2=IndexToGO.keySet().iterator();
          
          while(it2.hasNext()){
              int id=it2.next();
              String go=IndexToGO.get(id);
              GOtoIndex.put(go, id);
          }
  }
    
    public void createCOGGOMapping(File input){
int lineNum=0,index=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          if(lineNum>1){
          String[] st=line.split("\t");
          
          String ogName=st[0];
          String OGName="";
          int ogNumber=Integer.parseInt(ogName);
          
          int numTmpGo=ogNumber, nDigGo=0;
          
          if(numTmpGo<0){
              numTmpGo=-numTmpGo;
              OGName+="NOG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<6){
                  OGName+="0";
                  nDigGo++;
              }
              
              OGName+=(-ogNumber);
          }
          else{
              OGName+="COG";
              while(numTmpGo>0){
                  numTmpGo=numTmpGo/10;
                  nDigGo++;
              }
              
              while(nDigGo<4){
                  OGName+="0";
                  nDigGo++;
              }
              
              OGName+=ogNumber;
          }
           if(!CogGOmap.containsKey(OGName)){
                ArrayList<String> Gofunc=new ArrayList<String>();
                CogGOmap.put(OGName, Gofunc);
           }
          String go=st[1];

          String[] goSubcategories=go.split("@");
          for(int i=0;i<goSubcategories.length;i++){
     
              
               String GOName="GO";
                 int nNumGo=Integer.parseInt(goSubcategories[i]);
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
              
              CogGOmap.get(OGName).add(GOName);

              if(!GOtoIndex.containsKey(GOName)){
            GOtoIndex.put(GOName, index++);
            IndexToGO.put(index-1, GOName);
              }
          }
      }
       }
            reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         } 
    
    
    public void createCOGGOMappingOldOK(File input){
      int lineNum=0,index=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          if(lineNum>1309){
          String[] st=line.split("\\|");

          String COGname=st[0].substring(1, st[0].length());
           if(!CogGOmap.containsKey(COGname)){
           ArrayList<String> Gofunc=new ArrayList<String>();
           CogGOmap.put(COGname, Gofunc);
           }
          String[] go=st[st.length-1].split(" ");

          String[] goSubcategories=go[go.length-1].split("@");
          for(int i=0;i<goSubcategories.length;i++){
              
              CogGOmap.get(COGname).add(goSubcategories[i]);

              if(!GOtoIndex.containsKey(goSubcategories[i])){
            GOtoIndex.put(goSubcategories[i], index++);
            IndexToGO.put(index-1, goSubcategories[i]);
              }
          }
       }
      }
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         //add E.coli additional GO annotations
         File anEnrFile=new File("Mappings-E.coli.txt");
       
         int test=0;
         
         if(test==0){
         
         HashMap<String,ArrayList<String>> anEnr=new HashMap<>();
         
         try {
      Path path =Paths.get(anEnrFile.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
     lineNum=0;
      while ((line = reader.readLine()) != null) {
          if(lineNum==0){
              lineNum=1;
              continue;
          }
          if(line.contains("cog")){
              String cName="COG";
              String tmp[]=line.split("\t");
              String n=tmp[0].replace("cog", "");
              int nNum=Integer.parseInt(n);
              int numTmp=nNum, nDig=0;
  
              while(numTmp>0){
                  numTmp=numTmp/10;
                  nDig++;
              }
              
              while(nDig<4){
                  cName+="0";
                  nDig++;
              }
              
             cName+=nNum;
              
             String GOs[]=tmp[1].split(", ");
             ArrayList<String> goS=new ArrayList<>();
             for(String g:GOs){
                 String GOName="GO";
                 int nNumGo=Integer.parseInt(g);
              int numTmpGo=nNumGo, nDigGo=0;
              
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
            anEnr.put(cName, goS);     
          }
          
        for(String s:CogGOmap.keySet()){
            if(!anEnr.containsKey(s))
                continue;
            
            ArrayList<String> gos=CogGOmap.get(s);
            ArrayList<String> enGos=anEnr.get(s);
            
            for(int i=0;i<enGos.size();i++)
                if(!gos.contains(enGos.get(i)) && GOtoIndex.containsKey(enGos.get(i)))
                    gos.add(enGos.get(i));
            
        }    
      } 
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         } 
         
    }
    
    
     public void createCOGGOMapping1(File input, int test){
int lineNum=0,index=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          if(lineNum>1309){
          String[] st=line.split("\\|");

          String COGname=st[0].substring(1, st[0].length());
           if(!CogGOmap.containsKey(COGname)){
           ArrayList<String> Gofunc=new ArrayList<String>();
           CogGOmap.put(COGname, Gofunc);
           }
          String[] go=st[st.length-1].split(" ");

          String[] goSubcategories=go[go.length-1].split("@");
          for(int i=0;i<goSubcategories.length;i++){
     
              CogGOmap.get(COGname).add(goSubcategories[i]);

              if(!GOtoIndex.containsKey(goSubcategories[i])){
            GOtoIndex.put(goSubcategories[i], index++);
            IndexToGO.put(index-1, goSubcategories[i]);
              }
          }
       }
      }
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         //add E.coli additional GO annotations
         File anEnrFile=new File("Mappings-E.coli.txt");
                
         if(test==0){
         
         HashMap<String,ArrayList<String>> anEnr=new HashMap<>();
         
         try {
      Path path =Paths.get(anEnrFile.getAbsolutePath());

      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
     lineNum=0;
      while ((line = reader.readLine()) != null) {
          if(lineNum==0){
              lineNum=1;
              continue;
          }
          if(line.contains("cog")){
              String cName="COG";
              String tmp[]=line.split("\t");
              String n=tmp[0].replace("cog", "");
              int nNum=Integer.parseInt(n);
              int numTmp=nNum, nDig=0;
  
              while(numTmp>0){
                  numTmp=numTmp/10;
                  nDig++;
              }
              
              while(nDig<4){
                  cName+="0";
                  nDig++;
              }
              
             cName+=nNum;
              
             String GOs[]=tmp[1].split(", ");
             ArrayList<String> goS=new ArrayList<>();
             for(String g:GOs){
                 String GOName="GO";
                 int nNumGo=Integer.parseInt(g);
              int numTmpGo=nNumGo, nDigGo=0;
              
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
            anEnr.put(cName, goS);     
          }
          
        for(String s:CogGOmap.keySet()){
            if(!anEnr.containsKey(s))
                continue;
            
            ArrayList<String> gos=CogGOmap.get(s);
            ArrayList<String> enGos=anEnr.get(s);
            
            for(int i=0;i<enGos.size();i++)
                if(!gos.contains(enGos.get(i)) && GOtoIndex.containsKey(enGos.get(i)))
                    gos.add(enGos.get(i));
            
        }    
      } 
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         } 
         
    }
    
    void printGI(String pString){
        Iterator<String> keySetIterator = GOtoIndex.keySet().iterator();
        Path path= Paths.get(pString);
        try{
        BufferedWriter writer = Files.newBufferedWriter(path,ENCODING);
        while(keySetIterator.hasNext()){
            String go=keySetIterator.next();
            int index=GOtoIndex.get(go);
            writer.write(go+" "+index+"\n");
        }
        writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    public void printCOGMappings(){
        Iterator<String> it=CogGOmap.keySet().iterator();
        
        while(it.hasNext()){
            String cog=it.next();
            ArrayList<String> gos=CogGOmap.get(cog);
            
            System.out.println("COG: "+cog+"\n");
            for(int i=0;i<gos.size();i++)
                   System.out.print(gos.get(i)+" ");
            System.out.println("\n");
        }
    }
    
    void printCOG(String pString){
        Iterator<String> keySetIterator = CogGOmap.keySet().iterator();
        Path path= Paths.get(pString);
        try{
        BufferedWriter writer = Files.newBufferedWriter(path,ENCODING);
        while(keySetIterator.hasNext()){
            String cog=keySetIterator.next();
            writer.write(cog+"\n");
        }
        writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    public void saveCogGOMapping(String outFile){
        
        Iterator<String> keySetIterator = CogGOmap.keySet().iterator();
        Path path= Paths.get(outFile);
        try{
        BufferedWriter writer = Files.newBufferedWriter(path,ENCODING);
        while(keySetIterator.hasNext()){
            String cog=keySetIterator.next();
            ArrayList<String> sTmp= CogGOmap.get(cog);
            
            writer.write(cog+" ");
            
            for(String s:sTmp)
                writer.write(s+" ");;
                writer.write("\n");
        }
        writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

}
