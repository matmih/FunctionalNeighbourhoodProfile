/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package semanticsimilarity;

import eukaryoticdatasetcreator.OGGOMapping;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to compute different measures of semantic similarity between all pairs of GO terms in eukaryotes
 */
public class MainEucariots {
    static public void main(String args[]){
		
        int type = Integer.parseInt(args[0]);
        
         File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
         if(type==1)
             ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
         OGGOMapping oggomap=new OGGOMapping();
        
        oggomap.createCOGGOMapping(ogFuncFile);
               
        GOMap map1=new GOMap();
        map1.CreateGOMap(oggomap);
        map1.printGOMapNumeric();
        File output=null;
        if(type==0)
            output=new File("GosimFungi.txt");
        else output = new File("GosimMetazoa.txt");
        
        map1.saveGO(output);
        File fileWithUniprotFrequencyCounts=null;
        
        if(type ==0)
             fileWithUniprotFrequencyCounts=new File("GOFrequencyFungyN.txt");
        else
           fileWithUniprotFrequencyCounts= new File("GOFrequencyMetazoaN.txt");
        
        map1.loadFrequencies1(fileWithUniprotFrequencyCounts);
        map1.printFrequencies();
        
        try{
           File input = null;
           if(type==0)
               input = new File("hierarchyFungy.txt");
           else input = new File("hierarchyMetazoa.txt");
           
           java.nio.file.Path p = Paths.get(input.getAbsolutePath());
           BufferedReader buf = Files.newBufferedReader(p, StandardCharsets.UTF_8);
           
           String line = "";
           
           while((line = buf.readLine())!=null){
               map1.hierarchy = line;
           }
           
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        Graph g=new Graph(map1.GOmap.keySet().size());
        g.createAdjacency1(map1);
        File outputM=null;
        if(type==0)
             outputM=new File("AdjMatFungi.txt");
        else  outputM=new File("AdjMatMetazoa.txt");
        g.saveAdjacency(outputM, map1.GOmap.keySet().size());
        
        map1.printGOMap();

        ArrayList<Path> Paths=new ArrayList<Path>();
        Path visited=new Path();
         
          Paths.clear();
         visited.path.clear();
        
        Wang w=new Wang();
     
        GOSimilarity sim=new GOSimilarity((map1.GOmap.keySet().size()));//-1
              
        Object GO[]=map1.GOmap.keySet().toArray();
                
        int missing=0;
        for(int i=0;i<GO.length;i++){

            if(map1.GOmap.containsKey(GO[i]) && !GO[i].equals("root")){
           sim.addSimilarity(((int)map1.GOmap.get(GO[i]))+1, ((int)map1.GOmap.get(GO[i]))+1, g.computeSimilarity(map1.GOmap.get(GO[i]), map1.GOmap.get(GO[i]), map1));
            }
            else{
                missing++;
              }
        }
               
        for(int i=0;i<GO.length-1;i++){
            for(int j=i+1;j<GO.length;j++){
                if(map1.GOmap.containsKey(GO[i]) && map1.GOmap.containsKey(GO[j]) && !GO[i].equals("root") && !GO[j].equals("root")){
        sim.addSimilarity(((int)map1.GOmap.get(GO[i])+1), ((int)map1.GOmap.get(GO[j])+1), g.computeSimilarity(map1.GOmap.get(GO[i]), map1.GOmap.get(GO[j]), map1));  
                }
            }
            if(i%100==0)
                 System.out.println("index: "+i);
        }
           
        for(int i=0;i<GO.length;i++)
            if(map1.GOmap.containsKey(GO[i]) && !GO[i].equals("root")){
                sim.addWangSimilarity(((int)map1.GOmap.get(GO[i]))+1, ((int)map1.GOmap.get(GO[i]))+1, w.sGO(map1.GOmap.get(GO[i]), map1.GOmap.get(GO[i]), 1.0, g));
            }
        
        for(int i=0;i<GO.length;i++){
            for(int j=i+1;j<GO.length;j++){
        if(map1.GOmap.containsKey(GO[i]) && map1.GOmap.containsKey(GO[j]) && !GO[i].equals("root") && !GO[j].equals("root"))
          sim.addWangSimilarity(((int)map1.GOmap.get(GO[i]))+1, ((int)map1.GOmap.get(GO[j]))+1, w.sGO(map1.GOmap.get(GO[i]), map1.GOmap.get(GO[j]), 1.0, g)); 
            }
             if(i%100==0)
                 System.out.println("index: "+i);
        }
        
        File []outSim=new File[4];
        if(type==0){
            outSim[0]=new File("PairwiseGOSimFungi.txt");
         outSim[1]=new File("PairwiseGOSimLinFungi.txt");
         outSim[2]=new File("PairwiseGOSimRelFungi.txt");
         outSim[3]=new File("PairwiseGOSimWangFungi.txt");
        }
        else{
            outSim[0]=new File("PairwiseGOSimMetazoa.txt");
         outSim[1]=new File("PairwiseGOSimLinMetazoa.txt");
         outSim[2]=new File("PairwiseGOSimRelMetazoa.txt");
         outSim[3]=new File("PairwiseGOSimWangMetazoa.txt");
        }
        sim.writeSimilarity(outSim);
    }
}
