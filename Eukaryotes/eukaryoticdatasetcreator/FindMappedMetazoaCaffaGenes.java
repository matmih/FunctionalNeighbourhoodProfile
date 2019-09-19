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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to find Caffa genes that are contained in the metazoa mapping
 */
public class FindMappedMetazoaCaffaGenes {
   static public void main(String [] args){
       
        HashMap<String,HashSet<Pair<String,String>>> geneOGMap = new HashMap<>();
     
            String geneOGMapFile = "geneOGMappingFinalOKNew.txt";
       
       BufferedReader reader=null;
            
             try {
                     Path path =Paths.get(new File(geneOGMapFile).getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split(" ");
                            String gene = tmp[0].trim();
                            
                            if(!geneOGMap.containsKey(gene)){
                                geneOGMap.put(gene, new HashSet<Pair<String,String>>());
                            }
                            
                            HashSet<Pair<String,String>> ogs=geneOGMap.get(gene);
                            
                            for(int i=1;i<tmp.length;i++){
                                String og[]=tmp[i].split(":");
                               Pair<String, String> np = new Pair(og[0].trim(),og[1].trim());
                                
                               ogs.add(np);
                            }
                            geneOGMap.put(gene, ogs);
                    }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}
       
               File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
               
        OGGOMapping oggomap=new OGGOMapping();       
        oggomap.createCOGGOMapping(ogFuncFile);
       
        File caffaGenes = new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\metazoaCaffaGenes.txt");
        
        HashMap<String,Triplet<String,Integer, String>> caffa = new HashMap<>();
        HashSet<String> nogs = new HashSet<>();
        
        try{
            Path path = Paths.get(caffaGenes.getAbsolutePath());
            reader = Files.newBufferedReader(path,ENCODING);
            String line = "";
            while((line = reader.readLine())!=null){
                String tmp[] = line.split(" ");
                String caf = tmp[0].trim();
                String encaf = tmp[1].trim();
                String t1[] = encaf.split("\\.");
                int org = Integer.parseInt(t1[0]);
                String prot = t1[1].trim();
                Pair<String,Integer> p = new Pair(prot,org);
                if(geneOGMap.containsKey(prot)){
                    HashSet<Pair<String,String>> ogs = geneOGMap.get(prot);
                    for(Pair<String,String> ps:ogs){
                        if(Integer.parseInt(ps.getValue1()) == org){
                            if(oggomap.CogGOmap.containsKey(ps.getValue0())){
                                nogs.add(ps.getValue0());
                                Triplet<String,Integer,String> t = new Triplet(p.getValue0(),p.getValue1(),ps.getValue0()); 
                                     caffa.put(caf, t);
                            }
                        }
                    }
                }
                
            }
            reader.close();
            System.out.println("Different ogs: "+nogs.size());
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        File output = new File("caffaAnnotatedMetazoa.txt");
        
        try{
            FileWriter fw = new FileWriter(output);
            
            Iterator<String> cafIt = caffa.keySet().iterator();
            
            while(cafIt.hasNext()){
                String cafID = cafIt.next();
                
                Triplet<String,Integer,String> t = caffa.get(cafID);
                
                fw.write(cafID+" "+t.getValue0()+" "+t.getValue1()+" "+t.getValue2()+"\n");
                
            }
          fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
   }
}
