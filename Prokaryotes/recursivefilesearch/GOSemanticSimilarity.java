/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about GO semantic similarity
 */
public class GOSemanticSimilarity {
    
    double similarity[][];
    int numGOs=0;
    
    GOSemanticSimilarity(int num){
        numGOs=num;
        similarity=new double[num][num];
    }
    
    void readSimilarity(File input){
        BufferedReader reader;
         
        try {
            Path path =Paths.get(input.getAbsolutePath());

            reader = Files.newBufferedReader(path,ENCODING);
            String line = null;
            int count=0;
             while ((line = reader.readLine()) != null) {
                 String tmp[]=line.split(" ");
                 for(int i=0;i<tmp.length;i++)
                     similarity[count][i]=Double.parseDouble(tmp[i]);
                 count++;
      }
             reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
    }
    
    void printSimilarity(String GO, COGGOMap mapping){
        int index=mapping.GOtoIndex.get(GO);
        
        System.out.println("Printing semantic similarity vector for GO: "+GO);
        System.out.println(); System.out.println();
        
        for(int i=0;i<this.numGOs;i++)
            System.out.print(this.similarity[index][i]+" ");
        
    }
    
}
