/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create files to draw comparative voilin plots
 */
public class CreateFileForViolinPlot {
    static public void main(String args[]){
        
        ArrayList<ArrayList<Double>> methods=new ArrayList<>();
        
        try (BufferedReader bufRdr = new BufferedReader(new FileReader("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\Biological_process_data.dat")))
        {
            String line;
            int gcount=0;
            HashSet<Integer> gos=new HashSet<>();
            while ((line = bufRdr.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");
                
                int go = Integer.parseInt(parts[0]);
                if(!gos.contains(go)){
                    gos.add(go);
                    if(methods.size()==0){
                        ArrayList<Double> tmp=new ArrayList<>();
                        methods.add(tmp);
                    }
                    
                }
                else {
                    gos.clear();
                    ArrayList<Double> tmp=new ArrayList<>();
                    methods.add(tmp);
                    gos.add(go);
                    gcount++;
                }
                double auprc = Double.parseDouble(parts[2]);
                methods.get(gcount).add(auprc);
            }
            bufRdr.close();
        }
        catch(Exception e){}
        
        
         try (BufferedReader bufRdr = new BufferedReader(new FileReader("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\Cellular_component_data.dat")))
        {
            String line;
            int gcount=0;
            HashSet<Integer> gos=new HashSet<>();
            while ((line = bufRdr.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");
                
                int go = Integer.parseInt(parts[0]);
                if(!gos.contains(go)){
                    gos.add(go);
                    
                }
                else {
                    gos.clear();
                    ArrayList<Double> tmp=new ArrayList<>();
                    gos.add(go);
                    gcount++;
                }
                double auprc = Double.parseDouble(parts[2]);
                methods.get(gcount).add(auprc);
            }
            bufRdr.close();
        }
        catch(Exception e){}
        
         try (BufferedReader bufRdr = new BufferedReader(new FileReader("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\Molecular_function_data.dat")))
        {
            String line;
            int gcount=0;
            HashSet<Integer> gos=new HashSet<>();
            while ((line = bufRdr.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");
                
                int go = Integer.parseInt(parts[0]);
                if(!gos.contains(go)){
                    gos.add(go);
                    
                }
                else {
                    gos.clear();
                    ArrayList<Double> tmp=new ArrayList<>();
                    gos.add(go);
                    gcount++;
                }
                double auprc = Double.parseDouble(parts[2]);
                methods.get(gcount).add(auprc);
            }
            bufRdr.close();
        }
        catch(Exception e){}
         
        String out="C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\MethodViolinComparisonNew.txt";
         try{
         FileWriter fw = new FileWriter(out);
         System.out.println("metods size: "+methods.size());
         for(int i=0;i<methods.get(0).size();i++){
             for(int j=0;j<methods.size();j++){
                 if((j+1)<methods.size())
                    fw.write(methods.get(j).get(i)+"\t");
                 else
                     fw.write(methods.get(j).get(i)+"\n");
             }
         }
         
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}    
        
        
    }
}
