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
import java.util.HashSet;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to extract taxIDs belonging to bacteria classified in a group of free-living or pathogenic in mammals
 */
public class ExtractEnvAndClinTaxID {
    public static void main(String []args){
        
        File input = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\ProTraits_precisionScores.txt");
        File outputClin = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDFreeLiving.txt");
        File outputEnv = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search New\\taxIDPathogenicMammals.txt");
        
        BufferedReader reader = null;
        
        Path p = Paths.get(input.getAbsolutePath());
        HashSet<Integer> clinTax = new HashSet<>();
        HashSet<Integer> envTax = new HashSet<>();
        
        try{
             reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
            
             String s = "";
             
             
             while((s=reader.readLine())!=null){
                    String tmp[] = s.split(";");
                    
                    if(tmp[2].equals("pathogenic_in_mammals")){
                        if(Double.parseDouble(tmp[28])>=0.95)
                            envTax.add(Integer.parseInt(tmp[0]));
                    }
                    else if(tmp[2].equals("habitat=freeliving")){
                        if(Double.parseDouble(tmp[28])>=0.95)
                            clinTax.add(Integer.parseInt(tmp[0]));
                    }
                    
                    
             }
             
             reader.close();
             
             FileWriter fw = new FileWriter(outputClin.getAbsoluteFile());
             
             for(int t:clinTax)
                 fw.write(t+"\n");
             
             fw.close();
             
             fw = new FileWriter(outputEnv.getAbsoluteFile());
             
             for(int t:envTax)
                 fw.write(t+"\n");
             
             fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        
    }
}
