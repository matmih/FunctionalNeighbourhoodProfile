/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
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
  *@description class to compute the average co-expressions for significantly enriched pairs of GO functions divided by semantic similarity and insigifnicantly
  *enriched pairs of functions having LOR values in the interval [-0.5,0.5]
 */
public class ComputeAverageCorrelationsEColi {
    public static void main(String args[]){
        
        int significant = 1;//use insignificant
        
        File inputLORS = null;
                
        if(significant == 1)
                inputLORS= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\SupplementaryMaterials\\SignificantAll01Prokaryot.txt");
        else 
               inputLORS= new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\SupplementaryMaterials\\InSignificantAll01Prokaryot.txt");
      
        File inputpid2OGMap = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\pid2ogs.txt");
        File inputpid2gene = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\pid2geneMapping.txt");
        File inputOG2GO = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
        File inputGeneCorr = new File("C:\\Users\\matej\\Downloads\\ecoli_compendium_data\\correlationsFiltered.txt");
        
        HashMap<String,HashMap<String,Double>> logORsClose = new HashMap<>();
         HashMap<String,HashMap<String,Double>> logORsDist = new HashMap<>();
        HashMap<String,String> pid2gene = new HashMap<>();
        HashMap<String, HashMap<String,Double>> geneCorr = new HashMap<>();
        HashMap<String, HashSet<String>> go2OG = new HashMap<>();
        
        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(inputOG2GO);
       
        Iterator<String> it= cgmap.CogGOmap.keySet().iterator();
        
        while(it.hasNext()){
            String og = it.next();
            ArrayList<String> gos = cgmap.CogGOmap.get(og);
            
            for(int i=0;i<gos.size();i++){
                if(!go2OG.containsKey(gos.get(i)))
                    go2OG.put(gos.get(i), new HashSet<>());
                go2OG.get(gos.get(i)).add(og);
            }
        }
        
        GeneOGMapping pidOGMap=new GeneOGMapping();
        
        pidOGMap.loadGOGMapping(inputpid2OGMap);
        BufferedReader read = null;
        
        try{
            Path p = Paths.get(inputLORS.getAbsolutePath());
            read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            
            String line = "";
            int lineCount = 0;
            
            while((line = read.readLine())!=null){
               
                if(lineCount == 0){
                    lineCount = 1;
                    continue;
                }
                    
                
                String tmp[] = line.split("\t");
                
                String go1 = tmp[0].trim();
                String go2 = tmp[5].trim();
                
                if(go1.equals(go2)){                   
                    continue;
                }
                
                Double lor = Double.parseDouble(tmp[7].trim().replaceAll(",", "."));
                Double rs = Double.parseDouble(tmp[12].trim().replaceAll(",", "."));
                
               if(rs>=6 && rs<8){
                  // System.out.println(go1+" "+go2+" "+lor+" "+rs);
                if(!logORsClose.containsKey(go1))
                    logORsClose.put(go1, new HashMap<>());
                logORsClose.get(go1).put(go2, lor);
               }
               else{
                   int t= 1;
                   if(t==1)
                   continue;
                  if(!logORsDist.containsKey(go1))
                    logORsDist.put(go1, new HashMap<>());
                logORsDist.get(go1).put(go2, lor); 
               }                
            }
            
            read.close();
            
             p = Paths.get(inputpid2gene.getAbsolutePath());
              read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
              
              while((line=read.readLine())!=null){
                  String tmp[] = line.split(" ");
                  String pid = tmp[0].trim().replace(":", "");
                  String gene = tmp[1].trim();
                  pid2gene.put(pid, gene);
              }
              
              read.close();
              
               p = Paths.get(inputGeneCorr.getAbsolutePath());
              read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
            
              while((line = read.readLine())!=null){
                  String tmp[] = line.split("\t");
                  String g1 =  tmp[0].trim();
                  String g2 = tmp[1].trim();
                  Double corr = Double.parseDouble(tmp[2].trim().replaceAll(",", "."));
                  if(!geneCorr.containsKey(g1))
                      geneCorr.put(g1, new HashMap<>());
                  geneCorr.get(g1).put(g2, corr);
              }
              
              read.close();
              
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        Double averageCorrelationClose = 0.0, averageCorrelationDist = 0.0;
        it =   logORsClose.keySet().iterator();
        double cnt1 = 0.0;
        
        while(it.hasNext()){
            String go1 = it.next();
            HashMap<String,Double> mt = logORsClose.get(go1);
            Iterator<String> it2 = mt.keySet().iterator();
            
            while(it2.hasNext()){
                String go2 = it2.next();
                double lor = mt.get(go2);
                
                if(lor>2 || significant == 0){
                    
                    HashSet<String> og2 = go2OG.get(go2);
                    HashSet<String> og1 = go2OG.get(go1);
                    
                    HashSet<String> pid1 = new HashSet<>();
                    HashSet<String> pid2 = new HashSet<>();
                    
                    
                    for(String o:og1){
                        if(pidOGMap.OGGeneMap.containsKey(o))
                                 pid1.addAll(pidOGMap.OGGeneMap.get(o));
                    }
                    
                     for(String o:og2){
                        if(pidOGMap.OGGeneMap.containsKey(o))
                                 pid2.addAll(pidOGMap.OGGeneMap.get(o));
                    }
                    
                     HashSet<String> gene1 = new HashSet<>();
                     HashSet<String> gene2 = new HashSet<>();
                     
                     for(String pd:pid1){
                         if(pid2gene.containsKey(pd))
                             gene1.add(pid2gene.get(pd));
                     }
                     
                     for(String pd:pid2){
                         if(pid2gene.containsKey(pd))
                             gene2.add(pid2gene.get(pd));
                     }
                     for(String g:gene1){
                           if(!geneCorr.containsKey(g))
                               continue;
                           HashMap<String,Double> t = geneCorr.get(g);
                          for(String g1:gene2){
                              if(t.containsKey(g1)){
                                  averageCorrelationClose+=t.get(g1);
                                  cnt1+=1.0;
                              }
                         }
                     }
                }
            }
        }
        
        
        it =   logORsDist.keySet().iterator();
        
        double cnt = 0.0;
        
        while(it.hasNext()){
            String go1 = it.next();
            HashMap<String,Double> mt = logORsDist.get(go1);
            Iterator<String> it2 = mt.keySet().iterator();
            
            while(it2.hasNext()){
                String go2 = it2.next();
                double lor = mt.get(go2);
                
                if(lor>2 || significant == 0){
                    
                    HashSet<String> og2 = go2OG.get(go2);
                    HashSet<String> og1 = go2OG.get(go1);
                    
                    HashSet<String> pid1 = new HashSet<>();
                    HashSet<String> pid2 = new HashSet<>();
                    
                    
                    for(String o:og1){
                        if(pidOGMap.OGGeneMap.containsKey(o))
                                 pid1.addAll(pidOGMap.OGGeneMap.get(o));
                    }
                    
                     for(String o:og2){
                        if(pidOGMap.OGGeneMap.containsKey(o))
                                 pid2.addAll(pidOGMap.OGGeneMap.get(o));
                    }
                    
                     HashSet<String> gene1 = new HashSet<>();
                     HashSet<String> gene2 = new HashSet<>();
                     
                     for(String pd:pid1){
                         if(pid2gene.containsKey(pd))
                             gene1.add(pid2gene.get(pd));
                     }
                     
                     for(String pd:pid2){
                         if(pid2gene.containsKey(pd))
                             gene2.add(pid2gene.get(pd));
                     }
                     for(String g:gene1){
                           if(!geneCorr.containsKey(g))
                               continue;
                           HashMap<String,Double> t = geneCorr.get(g);
                          for(String g1:gene2){
                              if(t.containsKey(g1)){
                                  averageCorrelationDist+=t.get(g1);
                                  cnt+=1.0;
                              }
                         }
                     }
                }
            }
        }
       
        System.out.println("AverageCorrelationClose: "+averageCorrelationClose/cnt1);
        System.out.println("AverageCorrelationDist: "+averageCorrelationDist/cnt);
        
    }
}


