/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import static eukaryoticdatasetcreator.OGGOMapping.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to compute information needed to compute Caffa validation scores (eukaryotes)
 */
 
public class CafaValidation {
        public static void main(String [] args){
        
        File inputPredictions = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\BaselineOGsk=4_oob.preds");
        File inputPredictionsBaseline = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\genomeMetazoaBB1Cafa.test.pred.arff");
        //File inputCaffaGeneOg = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\fungiInfo.txt");
        File inputCaffaGeneOg = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\metazoaInfo.txt");
        File inputCaffaGeneGO = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\EukaryoticDatasetCreator\\AllTargetAnnotations.txt");
        
        HashMap<String,ArrayList<Double>> cogPredictions = new HashMap<>(); 
        HashMap<String,HashMap<String,Double>> CaffaGenePredictions = new HashMap<>();
        HashMap<String,HashSet<String>> CaffaGeneTrueLabels = new HashMap<>();
        HashMap<String,HashSet<String>> CaffaInputGeneTrueLabels = new HashMap<>();
        HashMap<String,HashSet<String>> CaffaGeneCOG = new HashMap<>();

        HashMap<Integer,HashMap<Integer, ArrayList<Integer>>> outputResults=new HashMap<>();
        
         File ogFuncFile=new File("NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
           // ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");
        //File ogFuncFile=new File("NogFunctionsProp3_10_30NewF.txt");//metazoa
       /* OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);*/
         

             HashMap<String, HashMap<String,Double>> data = new HashMap<>();
    
    Path p = Paths.get(inputPredictions.getAbsolutePath());
    
    int mode = 1; // 0 - GFN, 1-baseline
    
    //if GFN
    if(mode==0){
    try{
        BufferedReader read = new BufferedReader(new FileReader(inputPredictions));
        
        String line="";
        while((line=read.readLine())!=null){
            String tmp[] = line.split("\t");
            HashMap<String,Double> goPreds = new HashMap<>();
            String cog = tmp[0].replaceAll("\"", "").trim();
            String preds[] = tmp[2].trim().split(";");
            String gos[] = tmp[1].trim().split("@");
            int count=0;
            for(String gG:gos){
                String gosub[] = gG.trim().split("/");
                String goM =gosub[gosub.length-1].trim(); 
                String goM1 = goM.replace("GO", "");
                goM = "GO:"+goM1;
                goPreds.put(goM.trim(), Double.parseDouble(preds[count++]));
            }
            
            data.put(cog, goPreds);
        }
        read.close();
    }
    catch(IOException e){
        e.printStackTrace();
    }
    }
    else if(mode == 1){
        HashMap<String,Double> goPreds = new HashMap<>();
        ///implement baseline COG->prediction reading
        try{
        BufferedReader read = new BufferedReader(new FileReader(inputPredictionsBaseline));
               // Files.newBufferedReader(p);
        ArrayList<String> gos = new ArrayList<>();
        String line="";
        int dataSection = 0, count =0;
        while((line=read.readLine())!=null){
            
            if(line.contains("@ATTRIBUTE")){
                if(!line.contains("-a-GO"))
                    continue;
                String tmp[] = line.split(" ");
                String go = tmp[1].replace("class-a-","");
                go = go.trim();
                gos.add(go);
            }
            
            if(line.contains("DATA")){
                dataSection = 1;
                continue;
            }
            
            if(dataSection ==1){
                 goPreds = new HashMap<>();
                String tmp[] = line.split(",");
                String cog =tmp[0].replaceAll("\"", "");
                cog = cog.trim();
                for(int i=2+gos.size();i<2+2*gos.size();i++){
                    double value = Double.parseDouble(tmp[i]);
                    goPreds.put(gos.get(count++), value);
                }
                count=0;
                data.put(cog, goPreds);
            }
        }
        read.close();
    }
    catch(IOException e){
        e.printStackTrace();
    }
    }
    
    try{
         p = Paths.get(inputCaffaGeneOg.getAbsolutePath());
       BufferedReader read = Files.newBufferedReader(p, ENCODING);
        
        String line = "";
        
        while((line = read.readLine())!=null){
            String tmp[] = line.split("\t");
            
            String caffaGene = tmp[0].trim();
            String cog = tmp[2].trim();
            
            if(!CaffaGeneCOG.containsKey(caffaGene)){
                CaffaGeneCOG.put(caffaGene, new HashSet<String>());
                CaffaGeneCOG.get(caffaGene).add(cog);
            }
            else{
                CaffaGeneCOG.get(caffaGene).add(cog);
            }
        }
        
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
    Iterator<String> it = CaffaGeneCOG.keySet().iterator();
    
    while(it.hasNext()){
        String gene = it.next();
        
        if(!CaffaGenePredictions.containsKey(gene)){
            CaffaGenePredictions.put(gene, new HashMap<String,Double>());
        }
        
        HashSet<String> cogs = CaffaGeneCOG.get(gene);
        
        for(String s:cogs){
            if(!data.containsKey(s))
                continue;
            HashMap<String,Double> goPreds = data.get(s);

            HashMap<String,Double> geneGOPreds = CaffaGenePredictions.get(gene);
            
            Iterator<String> it1 = goPreds.keySet().iterator();
            
            while(it1.hasNext()){
                String go = it1.next();
                
                if(!geneGOPreds.containsKey(go)){
                    geneGOPreds.put(go, goPreds.get(go)/cogs.size());
                }
                else{
                    double pred = geneGOPreds.get(go);
                    pred+=goPreds.get(go)/cogs.size();
                    geneGOPreds.put(go, goPreds.get(go));
                }     
            } 
        }      
    }
    
    //load original labels
    
    
    try{
        p = Paths.get(inputCaffaGeneGO.getAbsolutePath());
        BufferedReader read = Files.newBufferedReader(p, ENCODING);
        
        String line = "";
        
        while((line = read.readLine())!=null){
            String tmp[] = line.split("\t");
            String gene = tmp[0].trim();
            String go = tmp[1].replace(":", "").trim();
            
            if(!CaffaGeneTrueLabels.containsKey(gene)){
                CaffaGeneTrueLabels.put(gene, new HashSet<String>());
                CaffaGeneTrueLabels.get(gene).add(go);
            }
            else CaffaGeneTrueLabels.get(gene).add(go);
        }
        
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
    try{
    
        FileWriter fw = new FileWriter("geneCaffaPreds.txt");
    //convert data to go -> gene scores (to evaluate GO performance measures)
    HashMap<String,Integer> geneToIndex = new HashMap<>();
    HashMap<Integer,String> indexToGene = new HashMap<>();
    
    HashMap<String,ArrayList<Double>> classifierScores = new HashMap<>();
    HashMap<String,ArrayList<Integer>> originalLabels = new HashMap<>();
        
         it = CaffaGenePredictions.keySet().iterator();
        int count = 0;
        ArrayList<String> goOrder = new ArrayList<>();
        int countC = 0;
        while(it.hasNext()){//change to geneGO
            String gene = it.next();//.replaceAll("\"", "");
            System.out.println("gene read: "+gene);
            geneToIndex.put(gene, countC++);
            indexToGene.put(countC-1, gene);
            HashMap<String,Double> preds = CaffaGenePredictions.get(gene);
            if(preds.keySet().isEmpty())
                continue;
            
            Iterator<String> it1 = preds.keySet().iterator();
            ArrayList<Double> predsN = new ArrayList<>();
          
            if(count==0){
            while(it1.hasNext()){
                String go = it1.next();
                goOrder.add(go);
               predsN.add(preds.get(go));
                
            }
            
                for(int i=0;i<goOrder.size();i++){
                    classifierScores.put(goOrder.get(i), new ArrayList<>(Collections.nCopies(CaffaGenePredictions.keySet().size(), 0.0)));
                    originalLabels.put(goOrder.get(i), new ArrayList<>(Collections.nCopies(CaffaGenePredictions.keySet().size(), 0)));
                }
            }
            else{
                for(int i=0;i<goOrder.size();i++){
                    String go=goOrder.get(i);
                    predsN.add(preds.get(go));
                }
            }
            
            if(!CaffaGeneTrueLabels.containsKey(gene))
                continue;
            
              fw.write(gene+" ");
             for(int i=0;i<predsN.size();i++){
                  classifierScores.get(goOrder.get(i)).set(geneToIndex.get(gene), predsN.get(i));
                  String goPreds = goOrder.get(i).replace(":", "");
                  if(CaffaGeneTrueLabels.get(gene).contains(goPreds))
                      originalLabels.get(goOrder.get(i)).set(geneToIndex.get(gene), 1);
                  else  originalLabels.get(goOrder.get(i)).set(geneToIndex.get(gene), 0);
                 
                 if((i+1)<predsN.size())
                     fw.write(predsN.get(i)+" ");
                 else
                      fw.write(predsN.get(i)+"\n");
            }
            
            count=1;
            
        }
        fw.close();
        
        
        fw = new FileWriter("ClassifierScoresCaffa.txt");
        
        for(int i=0;i<goOrder.size();i++){
            fw.write(goOrder.get(i)+"\t");
            ArrayList<Double> t = classifierScores.get(goOrder.get(i));
            
            for(int j=0;j<t.size();j++){
                if((j+1)<t.size())
                    fw.write(t.get(j)+"\t");
                else fw.write(t.get(j)+"\n");
            }
        }
        
        fw.close();
        
        fw = new FileWriter("TrueLabelsCaffa.txt");
        
         for(int i=0;i<goOrder.size();i++){
            fw.write(goOrder.get(i)+"\t");
            ArrayList<Integer> t = originalLabels.get(goOrder.get(i));
            
            for(int j=0;j<t.size();j++){
                if((j+1)<t.size())
                    fw.write(t.get(j)+"\t");
                else fw.write(t.get(j)+"\n");
            }
        }
        
        fw.close();
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
    }
}
