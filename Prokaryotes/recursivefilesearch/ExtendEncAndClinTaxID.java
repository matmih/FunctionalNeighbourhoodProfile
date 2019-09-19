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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to extend the bacterial selection to all strends of a given species
 */
public class ExtendEncAndClinTaxID {
      public static void main(String []args){
        
        File input = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\TaxIDs.txt");
        File inputClin = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDFreeLiving.txt");
        File inputEnv = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDPathogenicMammals.txt");
        File outputClinInSet = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDFreeLivingInSet.txt");
        File outputEnvInSet =  new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDPathogenicMammalsInSet.txt");
        BufferedReader reader = null;
        
         HashMap<Integer,Integer> taxIDTranslationMap= new HashMap<>();
         String taxIDTranslation = "TaxIDTranslationFile.txt";
        
        File translationFile=new File(taxIDTranslation);
        
             try {
                     Path path =Paths.get(translationFile.getAbsolutePath());
                     reader = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader.readLine()) != null) {
                            String tmp[] = line.split("\t");
                            int speciesTax = Integer.parseInt(tmp[1]);
                            int strainTax = Integer.parseInt(tmp[2]);
                            taxIDTranslationMap.put(strainTax, speciesTax);
                            }
      reader.close();
       }
         catch(Exception e){e.printStackTrace();}  
        
        
        Path p = Paths.get(input.getAbsolutePath());
        HashSet<Integer> clinTax = new HashSet<>();
        HashSet<Integer> envTax = new HashSet<>();
        HashSet<Integer> AllInTax = new HashSet<>();
         HashSet<Integer> clinInTax = new HashSet<>();
        HashSet<Integer> envInTax = new HashSet<>();
        
        try{
             reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
            
             String s = "";
             
             
             while((s=reader.readLine())!=null){
                    String tmp[] = s.split(" ");
                    AllInTax.add(Integer.parseInt(tmp[1].trim()));   
             }
             
             reader.close();
             
              p = Paths.get(inputClin.getAbsolutePath());
              reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
              
              while((s=reader.readLine())!=null){
                   
                    clinTax.add(Integer.parseInt(s.trim()));   
             } 
             
              reader.close();
              
               p = Paths.get(inputEnv.getAbsolutePath());
              reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
              
              while((s=reader.readLine())!=null){
                   
                    envTax.add(Integer.parseInt(s.trim()));   
             } 
             
              reader.close();
            
             for(int t:AllInTax){
                 if(taxIDTranslationMap.containsKey(t)){
                         
                 int spec = taxIDTranslationMap.get(t);
                 
                 for(int t1:clinTax){
                     if(!taxIDTranslationMap.containsKey(t1))
                         continue;
                     if(taxIDTranslationMap.get(t1) == spec){
                         clinInTax.add(t);
                         break;
                     }
                 }
             
                 for(int t1:envTax){
                     if(!taxIDTranslationMap.containsKey(t1))
                         continue;
                     if(taxIDTranslationMap.get(t1) == spec){
                         envInTax.add(t);
                         break;
                     }
                 }
                }
                 else{
                     
                     if(clinTax.contains(t))
                         clinInTax.add(t);
                     
                     if(envTax.contains(t))
                         envInTax.add(t);
                     
                 }
             } 
              
             FileWriter fw = new FileWriter(outputClinInSet.getAbsoluteFile());
             
             for(int t:clinInTax)
                 fw.write(t+"\n");
             
             fw.close();
             
             fw = new FileWriter(outputEnvInSet.getAbsoluteFile());
             
             for(int t:envInTax)
                 fw.write(t+"\n");
             
             fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        
    }
}
