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

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the Jaccard index at the gene level between pairs of functions on the E.Coli bacteria
 */
public class ComputeJSGeneLevelEColi {
    static public void main(String [] args){
        
        HashMap<String,HashSet<String>> goGenes = new HashMap<>();
        HashSet<String> genes = new HashSet<>();
        
        File input = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\2014-11-26_6.B_subtilis_168+18.E_coli_MG1655.goa");
        
        Path p = Paths.get(input.getAbsolutePath());
        try{
             BufferedReader read = Files.newBufferedReader(p, StandardCharsets.UTF_8);
             String line="";
             int count=0;
             while((line=read.readLine())!=null){
                 String tmp[] = line.split("\t");
                 String gene = tmp[1].trim();
                 String go = tmp[4].trim().replace(":", "");
                 
                 genes.add(gene);
                 
                 if(!goGenes.containsKey(go)){
                     goGenes.put(go, new HashSet<String>());
                     goGenes.get(go).add(gene);
                 }
                 else{
                     HashSet<String> tmpS = goGenes.get(go);
                     tmpS.add(gene);
                     goGenes.put(go, tmpS);
                 }
                 
                //// System.out.println(gene+" "+go);
                // count++;
                 
                 //if(count==6)
                   //  break;
             }
             read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        ArrayList<String> go1 = new ArrayList<>();
        go1.add("GO0005975"); go1.add("GO0006260"); go1.add("GO0006974");
        go1.add("GO0016051"); go1.add("GO0006457"); go1.add("GO0046700");
        go1.add("GO0006310"); go1.add("GO0046903"); go1.add("GO0008610");
        go1.add("GO0006865"); go1.add("GO0006974"); go1.add("GO0006508");
        
        ArrayList<String> go2 = new ArrayList<>();
        go2.add("GO0008643"); go2.add("GO0032506"); go2.add("GO0006265");
        go2.add("GO0043163"); go2.add("GO0016226"); go2.add("GO0051180");
        go2.add("GO0006952"); go2.add("GO0006935"); go2.add("GO0051668");
        go2.add("GO0009310"); go2.add("GO0033866"); go2.add("GO0019682");
       
        for(int i=0;i<go1.size();i++){
            String g1 = go1.get(i);
            String g2 = go2.get(i);
            
            HashSet<String> gen1 = goGenes.get(g1);
            HashSet<String> gen2 = goGenes.get(g2);
            
            if(gen1==null)
                gen1=new HashSet<>();
            if(gen2 ==null)
                gen2=new HashSet<>();
            
            HashSet<String> intersection = new HashSet();
            
            for(String s: gen1)
                if(gen2.contains(s))
                    intersection.add(s);
            
            System.out.println(g1+" "+g2+" "+gen1.size()+" "+gen2.size()+" "+intersection.size());
            
            if(gen1.size()>0 && gen2.size()>0)
              System.out.println("GO1: "+g1+" - GO2: "+g2+" JS: "+(intersection.size()/((double)gen1.size()+gen2.size()-intersection.size())));
            else
               System.out.println("NC! GO1: "+g1+" - GO2: "+g2+" JS: "+0.0); 
            
        }
        
    }
}
