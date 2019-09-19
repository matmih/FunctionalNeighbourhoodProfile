/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing mappings
 */
public class Mappings {
    
    public HashMap<String,HashSet<Pair<String,String>>> mappings=new HashMap<>(); //all genes mappings from file eggnog proteinid conversion
    public HashMap<String,HashSet<String>> uniprotMappings = new HashMap<>();
    
    public void loadMappings(File input){
        
        try{
        FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     BufferedReader reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      int count=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                          String tmp[]=line.split("\t");
                          String id1=tmp[1].trim();
                          String id2=tmp[2].trim();
                          String tax=tmp[0].trim();
                                 
                            if(id1.contains(".") && !tmp[3].contains("KEGG_ec")){
                                String tmpt[]=id1.split("\\.");
                                 if(tmpt.length==3 && id1.trim().matches("SP[^.]+\\.[^.]+\\.\\d+")){ 
                                     id1=tmpt[0].trim()+"."+tmpt[1].trim();//id1.trim();
                                     System.out.println("mapping: ");
                                     System.out.println(id1);
                                  }
                                 else{
                                   
                                 if(tmpt.length<2)
                                    continue;
                                      id1=tmpt[0].trim();
                                 }
                            }
                            
                            if(id1.contains("ncr:") && id1.contains("NCU") && !id1.contains(".")){
                                id1=id1.replace("ncr:", "");
                            }
                            
                            if(id2.contains("ncr:") && id2.contains("NCU") && !id2.contains(".")){
                                id2=id2.replace("ncr:", "");
                            }
                            
                             if(id2.contains(".") && !tmp[3].contains("KEGG_ec")){
                                 
                                 String tmpt[]=id2.split("\\.");
                                 if(tmpt.length==3 && id2.trim().matches("SP[^.]+\\.[^.]+\\.\\d+")){ 
                                     id2=tmpt[0].trim()+"."+tmpt[1].trim();//id2.trim();
                                       System.out.println("mapping: ");
                                     System.out.println(id2);
                                  }
                                 else{
                                
                                 if(tmpt.length<2)
                                    continue;
                                id2=tmpt[0].trim();
                                         }
                            }
                          
                          if(!mappings.containsKey(id1)){
                              mappings.put(id1, new HashSet<Pair<String,String>>());
                              mappings.get(id1).add(new Pair(id2,tax));
                               mappings.get(id1).add(new Pair(id1,tax));
                          }
                          else if(mappings.containsKey(id1)){
                              HashSet<Pair<String,String>> s=mappings.get(id1);
                              s.add(new Pair(id2,tax));
                              mappings.put(id1, s);
                          }
                          
                          if(!mappings.containsKey(id2)){
                              mappings.put(id2, new HashSet<Pair<String,String>>());
                              mappings.get(id2).add(new Pair(id1,tax));
                               mappings.get(id2).add(new Pair(id2,tax));
                          }
                          else if(mappings.containsKey(id2)){
                              HashSet<Pair<String,String>> s=mappings.get(id2);
                              s.add(new Pair(id1,tax));
                              mappings.put(id2, s);
                          }
                          count++;
                          
                          if(((double)count/5200000)*100%5==0)
                          System.out.println(((double)count/5200000)*100+"%");
                          
                    }
                      reader1.close();
        }
        catch(IOException e){
            e.printStackTrace();
        } 
        System.out.println("Translation mappings loaded!");
    }
    
    
       public void loadMappingsMetazoa(File input){
        
        try{
        FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     BufferedReader reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      int count=0;
                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                          String tmp[]=line.split("\t");
                          String id1=tmp[1].trim();
                          String id2=tmp[2].trim();
                          String tax=tmp[0].trim();
                          int taxI = Integer.parseInt(tax);
                                 
                            if(id1.contains(".") && taxI!=6239){
                                String tmpt[]=id1.split("\\.");
                                 
                                      id1=tmpt[0].trim();
                                 
                            }
                            
                            if(id1.contains("-PA")){
                                id1=id1.replace("-PA", "");
                            }
                            
                            if(id1.contains("__mRNA")){
                                id1=id1.replace("__mRNA", "");
                            }
                            
                            if(id1.contains("-tr")){
                                id1=id1.replace("-tr", "");
                            }
                            
                            if(id1.contains("ncr:") && id1.contains("NCU") && !id1.contains(".")){
                                id1=id1.replace("ncr:", "");
                            }
                            
                            if(id2.contains("ncr:") && id2.contains("NCU") && !id2.contains(".")){
                                id2=id2.replace("ncr:", "");
                            }
                            
                             if(id2.contains(".") && taxI!=6239){
                                 
                                 String tmpt[]=id2.split("\\.");
                                 
                                id2=tmpt[0].trim();
                             }
                            
                          if(!mappings.containsKey(id1)){
                              mappings.put(id1, new HashSet<Pair<String,String>>());
                              mappings.get(id1).add(new Pair(id2,tax));
                               mappings.get(id1).add(new Pair(id1,tax));
                          }
                          else if(mappings.containsKey(id1)){
                              HashSet<Pair<String,String>> s=mappings.get(id1);
                              s.add(new Pair(id2,tax));
                              mappings.put(id1, s);
                          }
                          
                          if(!mappings.containsKey(id2)){
                              mappings.put(id2, new HashSet<Pair<String,String>>());
                              mappings.get(id2).add(new Pair(id1,tax));
                               mappings.get(id2).add(new Pair(id2,tax));
                          }
                          else if(mappings.containsKey(id2)){
                              HashSet<Pair<String,String>> s=mappings.get(id2);
                              s.add(new Pair(id1,tax));
                              mappings.put(id2, s);
                          }
                          count++;
                          
                         // if(((double)count/21665229)*100%5==0)
                          if(count%1000==0)
                          System.out.println(((double)count/21665229)*100+"%");
                          
                    }
                      reader1.close();
        }
        catch(IOException e){
            e.printStackTrace();
        } 
        System.out.println("Translation mappings loaded!");
    }
    
    
            public void loadUniprotMappings(File input){
                    
                 try{
        FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     BufferedReader reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;

                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                          String tmp[]=line.split("\t");
                          String id1=tmp[1].trim();
                          String go=tmp[4].trim();
                           
                          if(tmp[3].toLowerCase().contains("not"))
                              continue;
                          
                          /*System.out.println(line);
                          System.out.println(id1);
                          System.out.println(go);*/
                          
                          if(!uniprotMappings.containsKey(id1)){
                              uniprotMappings.put(id1, new HashSet<String>());
                              uniprotMappings.get(id1).add(go);
                          }
                          else if(uniprotMappings.containsKey(id1)){
                              HashSet<String> s=uniprotMappings.get(id1);
                              s.add(go);
                             uniprotMappings.put(id1, s);
                          }
                    }
                      reader1.close();
        }
        catch(IOException e){
            e.printStackTrace();
        } 
            
              }  
}
