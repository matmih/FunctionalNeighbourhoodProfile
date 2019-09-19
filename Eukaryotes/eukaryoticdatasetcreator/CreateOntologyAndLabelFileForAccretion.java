/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create files for accrition computation
  * type = 0 -> fungi, type = 1 -> metazoa
 */
public class CreateOntologyAndLabelFileForAccretion{
    static public void main(String args[]){
        int type = 1;//Integer.parseInt(args[0]);//0-fungi, 1-metazoa
       
        String geneFunctionMapping = "C:\\Users\\matej\\Downloads\\Eukaryot data\\GeneFunctionsNew1.txt";//
        
        if(type == 1)
             geneFunctionMapping = "C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\GeneFunctionsNew1.txt";
        
        HashMap<String,HashSet<String>> genFunctions = new HashMap<>();
        
        Path p = Paths.get(geneFunctionMapping);
        
        try{
        BufferedReader read = Files.newBufferedReader(p,StandardCharsets.UTF_8);
        
        String line = "";
        
        while((line = read.readLine())!=null){
            String tmp[] = line.split(" ");
            String gen = tmp[0].trim();
            /*if(!genFunctions.containsKey(gen)){
                genFunctions.put(gen,new HashSet<String>());
            }*/
            
            for(int i=1;i<tmp.length;i++){
                String subtmp[] = tmp[i].split("::");
                String gen1 = gen+subtmp[1];
                if(!genFunctions.containsKey(gen1))
                    genFunctions.put(gen1, new HashSet<String>());
                genFunctions.get(gen1).add(subtmp[0].trim());
            }
        }
        
        read.close();
        
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        File output = null;
        
       
        
        String ontologyFile = "C:\\Users\\matej\\Downloads\\Eukaryot data\\go_daily-termdb.obo-xml\\go_daily-termdb.obo-xml";
        
        p = Paths.get(ontologyFile);
        
        HashMap<String,GOTerm> goTermMapping = new HashMap<>();
        HashSet<String> allgos = new HashSet<>();
        
           Iterator<String> Itk = genFunctions.keySet().iterator();
           
           while(Itk.hasNext()){
               HashSet<String> s =genFunctions.get(Itk.next());
               allgos.addAll(s);
           }
        
           
                if(type == 0)
            output = new File("labelsFungiIA.txt");
        else output = new File("labelsMetazoaIA.txt");
           
             try{
        GeneOntology go = new GeneOntology(ontologyFile);
        LinkedHashSet<GOTerm> allTerms = (LinkedHashSet)go.getAllTerms(false, false);
        FileWriter fw = new FileWriter("ontologyIA.txt");
        int count = 0;
       
        System.out.println("All terms size: "+allTerms.size());
        int count1=0;
        for(GOTerm g:allTerms){
            
             for(String ggo:allgos){
                 if(ggo.equals(g.getFormattedId())){
                     goTermMapping.put(ggo, g);
                     break;
                 }
                   
             }
     
             count1++;  
             if(count1%1000==0)
             System.out.println((((double)count1)/allTerms.size())+"%");
             
            
            Set<GOTerm> ch =  g.getChildren();
            
            for(GOTerm go2:ch)
                fw.write(g.getFormattedId()+"\t"+go2.getFormattedId()+"\n");
            
          
            
        }
      
            fw.close();
            
            
          fw = new FileWriter(output);
            
            Itk=genFunctions.keySet().iterator();
            int numC = genFunctions.keySet().size(), iter = 0;
            System.out.println("Num genes: "+numC);
            HashSet<String> gosN = null;
               while(Itk.hasNext()){
                String f = Itk.next();
                
                HashSet<String> gos = genFunctions.get(f);
                    gosN = new HashSet<>(1000);
                for(String ggg:gos){
                    if(!goTermMapping.containsKey(ggg))
                        continue;
                    GOTerm gc = goTermMapping.get(ggg);
                    Set<GOTerm> par = gc.getAllParents();
                    
                    for(GOTerm d:par)
                        gosN.add(d.getFormattedId());    
                }
                
                count+=gosN.size();
                gosN.addAll(gos);
                
                 for(String gogo: gosN)
                     fw.write(f+"\t"+gogo+"\t"+1+"\n"); 
                
                iter++;
                
              //  genFunctions.put(f, gos);
                if(iter%1000==0){
                System.out.println((((double)iter/numC)*100)+"%");
                System.out.println("Num new gens: "+count);
                }
            }
            
                System.out.println("Newly added: "+count); 
               fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
             
         
        
       /* try{
             FileWriter fw = new FileWriter(output);
             
             Iterator<String> it = genFunctions.keySet().iterator();
             
             while(it.hasNext()){
                 String g = it.next();
                 
                 HashSet<String> gos = genFunctions.get(g);
      
                 for(String go: gos)
                     fw.write(g+"\t"+go+"\t"+1+"\n");                
             }         
            fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }*/
        
             
        
    }
}
