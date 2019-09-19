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
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class compute gene co-expressions from colombos data
 */
 
public class ComputePairwiseGeneCorrelation {
    public static void main(String args[]){
        
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.REMOVED);
        SpearmansCorrelation sc = new SpearmansCorrelation(ranking);
       
        double input[][] = new double[4080][4321];
        String inGenes[] = new String[4321];
        int filter = 1;
        
        BufferedReader reader = null;
        File inF = new File("C:\\Users\\matej\\Downloads\\ecoli_compendium_data\\colombos_ecoli_exprdata_20151029Tmp.txt");
        Path p= Paths.get(inF.getAbsolutePath());
       
        try{
            reader = Files.newBufferedReader(p, StandardCharsets.UTF_8);
            
            String line = "";
            int count = 0, countNans = 0;
            
            while((line=reader.readLine())!=null){
                countNans=0;
                 String tmp[]= line.split("\t");
                 inGenes[count] = tmp[1].trim();
                 
                 for(int i=3;i<tmp.length;i++)
                     if(!tmp[i].equals("NaN"))
                         countNans++;
                
                 if(filter == 1)
                 if((countNans/(double)(tmp.length-3)<0.5))
                     continue;
                 
                 for(int i=3;i<tmp.length;i++){
                     if(!tmp[i].equals("NaN"))
                         input[i-3][count] = Double.parseDouble(tmp[i].trim());
                     else input[i-3][count] = Double.NaN;
                 }
                 count++;

            }         
            reader.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        RealMatrix res = null;
        res = sc.computeCorrelationMatrix(input);
        
        try{
                FileWriter fw = null;
                
                if(filter == 0)
                    fw = new FileWriter("C:\\Users\\matej\\Downloads\\ecoli_compendium_data\\correlations.txt");
                else 
                  fw = new FileWriter("C:\\Users\\matej\\Downloads\\ecoli_compendium_data\\correlationsFiltered.txt");
                
                for(int i=0;i<res.getColumnDimension();i++){
                    double tmp[] = res.getColumn(i);
                    
                    for(int j=0;j<tmp.length;j++){
                        fw.write(inGenes[i]+"\t"+inGenes[j]+"\t"+tmp[j]+"\n");
                    }
                }
                    
                fw.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    
}
