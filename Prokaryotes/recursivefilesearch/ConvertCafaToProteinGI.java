/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class create mapping between protein ids and caffa ids
 */
public class ConvertCafaToProteinGI {
    final static Charset ENCODING = StandardCharsets.UTF_8;
    static public void main(String [] args){
        HashMap<String,Pair<String,String>> caffaToNCBI = new HashMap<>();
        HashMap<String,HashSet<Pair<String,String>>> NCBIToCaffa = new HashMap<>();
        HashMap<String,String> geneToAccession = new HashMap<>();
        HashMap<String,Pair<String,String>> caffaToProteinGI = new HashMap<>();
        
        File inputCToNCBIF = new File("C:\\Users\\matej\\Downloads\\gene2accession\\CaffaToNCBITax.txt");
        File geneToAccessionF = new File("C:\\Users\\matej\\Downloads\\gene2accession\\gene2accession");
        File cafaToProteinGIF = new File("C:\\Users\\matej\\Downloads\\gene2accession\\caffaPIDT.txt");
        
        Path path = Paths.get(inputCToNCBIF.getAbsolutePath());
        
        BufferedReader read = null;
        try{
            read=Files.newBufferedReader(path,ENCODING);
            
            String line = "";
            while((line = read.readLine())!=null){
                String tmp[] = line.split("\t");
                Pair<String,String> p = new Pair(tmp[1].trim(),tmp[2].trim());
                caffaToNCBI.put(tmp[0].trim(), p);
                if(!NCBIToCaffa.containsKey(tmp[1].trim())){
                    HashSet<Pair<String,String>> tmpS = new HashSet<>();            
                    Pair<String,String> p1 = new Pair(tmp[0].trim(),tmp[2].trim());
                     tmpS.add(p1);
                    NCBIToCaffa.put(tmp[1].trim(), tmpS);
                }
                else{
                    HashSet<Pair<String,String>> tmpS = NCBIToCaffa.get(tmp[1].trim());
                    tmpS.add(new Pair(tmp[0].trim(),tmp[2].trim()));
                NCBIToCaffa.put(tmp[1].trim(), tmpS);
                }
            }
           read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        System.out.println("MSize: "+caffaToNCBI.keySet().size()+" "+NCBIToCaffa.keySet().size());
        
         path = Paths.get(geneToAccessionF.getAbsolutePath());
         try{
            read=Files.newBufferedReader(path,ENCODING);
            
            int count=0;
            String line = "";
            while((line = read.readLine())!=null){
                if(count==0){
                    count=1;
                    continue;
                }
                        
                String tmp[] = line.split("\t");
                String accession = tmp[5].trim();
                accession = tmp[5].split("\\.")[0].trim();
                String pid = tmp[6].trim();
                if(NCBIToCaffa.containsKey(accession)){
                    Pair<String,String> p = new Pair(pid, tmp[0].trim());
                    HashSet<Pair<String,String>> prots = NCBIToCaffa.get(accession);
                    for(Pair<String,String> pT: prots){
                         Pair<String,String> p2 = new Pair(pid,pT.getValue1());
                    caffaToProteinGI.put(pT.getValue0(), p2);
                    }
                }
               // caffaToNCBI.put(tmp[0].trim(), tmp[1].trim());
            }
            read.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
         
         System.out.println("Initial size: "+caffaToNCBI.keySet().size());
         System.out.println("New size: "+caffaToProteinGI.keySet().size());
        
        try{ 
         FileWriter fw = new FileWriter(cafaToProteinGIF);
         
         Iterator<String> it = caffaToProteinGI.keySet().iterator();
         
         while(it.hasNext()){
             String caffa = it.next();
             fw.write(caffa+"\t"+caffaToProteinGI.get(caffa).getValue0()+"\t"+caffaToProteinGI.get(caffa).getValue1()+"\n");
         }
         fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
         
    }
}
