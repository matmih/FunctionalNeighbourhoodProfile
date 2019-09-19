/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about predictions made by the classification algorithm.
 */
public class ClusPredictionInfo {
    
    public HashMap<String,ArrayList<Double>> classifierScores=new HashMap<>();//index in the array corresopnds to COG index in cgmap
    public HashMap<String,ArrayList<Double>> originalLabels=new HashMap<>();
    HashMap<String,Integer> cogToIndex=new HashMap<>();
    HashMap<Integer,String> indexToCog=new HashMap<>();
    ArrayList<String> gos = new ArrayList<>();
    
    public void readPredictionInfoSCXval(COGGOMap cgmap, String GO, File input){
        BufferedReader bufRdr1=null; 
         
         try{
         bufRdr1= new BufferedReader(new FileReader(input));
         }
         catch(IOException e){
             e.printStackTrace();
         }
        
          String line;
            int data=0;
            
            int cogCount=0;
            ArrayList<String> GOs=new ArrayList();
            GOs.add(GO);
         try{
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    continue;
                }
                else if(line.contains("@ATTRIBUTE"))
                    continue;
                else if(line.contains("@DATA")){
                    data=1;
                    continue;
                }
                else if(line.contains("RELATION"))
                    continue;
                
                if(line.equals(""))
                    continue;

                if(data==1){
                    
                    if(line.contains("Fold = "))
                        continue;
                    
                String tmp[]=line.split(",");
                String cog=tmp[0].replace("\"", "");

                if(!cogToIndex.containsKey(cog)){
                    cogToIndex.put(cog, cogCount++);
                    indexToCog.put(cogCount-1, cog);
                }

                ArrayList<Double> sc=new ArrayList<>(); 
                    
                    if(!classifierScores.containsKey(GOs.get(0))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            classifierScores.put(GOs.get(0), tmpA);
                    }
                    
                    if(!originalLabels.containsKey(GOs.get(0))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            originalLabels.put(GOs.get(0), tmpA);
                    }
                    
                    double score=Double.parseDouble(tmp[3]);
                    double score1=Double.parseDouble(tmp[1]);
                    classifierScores.get(GOs.get(0)).set(cogToIndex.get(cog), score);
                    originalLabels.get(GOs.get(0)).set(cogToIndex.get(cog), score1);
                }
              }
            
            bufRdr1.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }    
    }
    
    public void readPredictionInfo(COGGOMap cgmap, File input){
         BufferedReader bufRdr1=null; 
         
         try{
         bufRdr1= new BufferedReader(new FileReader(input));
         }
         catch(IOException e){
             e.printStackTrace();
         }
        
          String line;
            int data=0;
            
            int cogCount=0;
            ArrayList<String> GOs=new ArrayList();
         try{
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    func=func.replace("Original-p-", "");
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    System.out.println("func: "+func);
                    GOs.add(func);
                    gos.add(func);
                    continue;
                }
                else if(line.contains("@ATTRIBUTE"))
                    continue;
                else if(line.contains("@DATA")){
                    data=1;
                    continue;
                }
                else if(line.contains("RELATION"))
                    continue;
                
                if(line.equals(""))
                    continue;
                if(data==1){
                String tmp[]=line.split(",");
                String cog=tmp[0].replace("\"", "");
                if(!cogToIndex.containsKey(cog)){
                    cogToIndex.put(cog, cogCount++);
                    indexToCog.put(cogCount-1, cog);
                }
                ArrayList<Double> sc=new ArrayList<>(); 
                for(int i=cgmap.GOtoIndex.keySet().size()+2;i<2*cgmap.GOtoIndex.keySet().size()+2;i++){//1488-2975
                  
                    if(!classifierScores.containsKey(GOs.get(i-2-cgmap.GOtoIndex.keySet().size()))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            classifierScores.put(GOs.get(i-2-cgmap.GOtoIndex.keySet().size()), tmpA);
                    }
                    
                    if(!originalLabels.containsKey(GOs.get(i-2-cgmap.GOtoIndex.keySet().size()))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            originalLabels.put(GOs.get(i-2-cgmap.GOtoIndex.keySet().size()), tmpA);
                    }
                    
                    double score=Double.parseDouble(tmp[i]);
                    double score1=Double.parseDouble(tmp[i-cgmap.GOtoIndex.size()]);
                    classifierScores.get(GOs.get(i-2-cgmap.GOtoIndex.keySet().size())).set(cogToIndex.get(cog), score);
                    originalLabels.get(GOs.get(i-2-cgmap.GOtoIndex.keySet().size())).set(cogToIndex.get(cog), score1);
                }
               }
            }
            bufRdr1.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }       
    }
    
    public ArrayList<Double> sortedScores(ArrayList<Double> scores){
        
        ArrayList<Double> tmp=new ArrayList(scores);
        
        Collections.sort(tmp);
        
        return tmp;  
    }
    
    public void loadGOOrder(File input){
        BufferedReader bufRdr1 = null;
        
        try{
            bufRdr1 = new BufferedReader(new FileReader(input));
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        String line = "";
        
        try{
            while((line = bufRdr1.readLine())!=null){
                line = line.trim();
                System.out.println("line loadGOOrder: "+line+"\n");
                 if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    func=func.replace("Original-p-", "");
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    System.out.println("func: "+func);
                    gos.add(func);
                    continue;
                }
                else continue;
                
            }
            bufRdr1.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
    }
    
    public void readPredictionInfoOOB(COGGOMap cgmap, ARFF infile, File input){         
          HashMap<String, HashMap<String,Double>> data = new HashMap<>();
    
    Path p = Paths.get(input.getAbsolutePath());
    
    try{
        BufferedReader read = new BufferedReader(new FileReader(input));
        
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
    
        Iterator<String> it = data.keySet().iterator();
        int count = 0;
        ArrayList<String> goOrder = new ArrayList<>();
        int countC = 0;
        while(it.hasNext()){
            String cog = it.next();//.replaceAll("\"", "");
            System.out.println("cog read: "+cog);
            cogToIndex.put(cog, countC++);
            indexToCog.put(countC-1, cog);
            HashMap<String,Double> preds = data.get(cog);
            
            Iterator<String> it1 = preds.keySet().iterator();
            ArrayList<Double> predsN = new ArrayList<>();
          
            if(count==0){
            while(it1.hasNext()){
                String go = it1.next();
                goOrder.add(go);
               predsN.add(preds.get(go));
                
            }
            
                for(int i=0;i<goOrder.size();i++){
                    classifierScores.put(goOrder.get(i), new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0)));
                    originalLabels.put(goOrder.get(i), new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0)));
                }
            }
            else{
                for(int i=0;i<goOrder.size();i++){
                    String go=goOrder.get(i);
                    predsN.add(preds.get(go));
                }
            }

            int ind = infile.findCogIndex(cog);
            System.out.println("cog index: "+ind);
            if(ind==-1)
                continue;
            
           String tar = infile.Targets.get(ind);

             for(int i=0;i<predsN.size();i++){
                  classifierScores.get(goOrder.get(i)).set(cogToIndex.get(cog), predsN.get(i));
                  String goPreds = goOrder.get(i).replace(":", "");
                  if(tar.contains(goPreds)/*goOrder.get(i))*/)
                         originalLabels.get(goOrder.get(i)).set(cogToIndex.get(cog), 1.0);
                  else  originalLabels.get(goOrder.get(i)).set(cogToIndex.get(cog), 0.0);
            }           
            count=1;            
        }  
     }
    
    public ArrayList<Integer> numberOfNewGenes(ArrayList<Double> predicted, ArrayList<Double> trueClass, GeneOGMapping geneOGMap){
        int oGenes=0, nGenes=0, difference=0;
        HashSet<String> trueGenes=new HashSet<>();
        HashSet<String> predictedGenes=new HashSet<>();
        HashSet<String> newlypredictedGenes=new HashSet<>();     
        
        for(int i=0;i<trueClass.size();i++){
            if(trueClass.get(i)==1.0){
                HashSet<String> st=geneOGMap.getGenes(indexToCog.get(i));
                trueGenes.addAll(st);
            }
            if(predicted.get(i)==1.0){
                HashSet<String> st=geneOGMap.getGenes(indexToCog.get(i));
                predictedGenes.addAll(st);
            }
            
            if(trueClass.get(i)==0.0 && predicted.get(i)==1.0){
                HashSet<String> st=geneOGMap.getGenes(indexToCog.get(i));
                newlypredictedGenes.addAll(st);
            }
            
        }
        difference=predictedGenes.size()-trueGenes.size();
        ArrayList<Integer> res=new ArrayList<>();
        System.out.println("trueGenes size: "+trueGenes.size());
        res.add(trueGenes.size()); res.add(predictedGenes.size()); res.add(difference); res.add(newlypredictedGenes.size());
        
        return res;
    }
    
     public ArrayList<Integer> numberOfNewGenesOrganism(ArrayList<Double> predicted, ArrayList<Double> trueClass, GenesAndOgsOrganism organismMappings, HashSet<String> npgenes){//geneOGMap for only those genes in a organism, set containing only ogs contained in a organism
        int oGenes=0, nGenes=0, difference=0;
        HashSet<String> trueGenes=new HashSet<>();
        HashSet<String> predictedGenes=new HashSet<>();
        HashSet<String> newlypredictedGenes=new HashSet<>();
        HashSet<String> ogsOrg = organismMappings.organismOgs;
        
        for(int i=0;i<trueClass.size();i++){
            if(!ogsOrg.contains(indexToCog.get(i)))
                continue;
            if(trueClass.get(i)==1.0){
                HashSet<String> st= organismMappings.OGGeneMappingOrganism.get(indexToCog.get(i));
                trueGenes.addAll(st);
            }
            if(predicted.get(i)==1.0){
                HashSet<String> st=organismMappings.OGGeneMappingOrganism.get(indexToCog.get(i));
                predictedGenes.addAll(st);
            }
            
            if(trueClass.get(i)==0.0 && predicted.get(i)==1.0){
                HashSet<String> st=organismMappings.OGGeneMappingOrganism.get(indexToCog.get(i));
                newlypredictedGenes.addAll(st);
            }
            
        }
        npgenes.addAll(newlypredictedGenes);
        difference=predictedGenes.size()-trueGenes.size();
        ArrayList<Integer> res=new ArrayList<>();
        System.out.println("trueGenes size: "+trueGenes.size());
        res.add(trueGenes.size()); res.add(predictedGenes.size()); res.add(difference); res.add(newlypredictedGenes.size());
        
        return res;
    }
    
    public ArrayList<Double> predictedClass(ArrayList<Double> scores, double threshold){
        
        ArrayList<Double> tmp=new ArrayList<>();
        
        for(int i=0;i<scores.size();i++)
            if(scores.get(i)>=threshold)
                tmp.add(1.0);
            else tmp.add(0.0);
        
        return tmp;
    }
    
    public ArrayList<Double> precRec(ArrayList<Double> scores, ArrayList<Double> trueLabel){
        
        ArrayList<Double> precRec=new ArrayList<>();
        double tp=0.0,fp=0.0,tn=0.0,fn=0.0;
        
        for(int i=0;i<scores.size();i++){
            if(scores.get(i)==1 && trueLabel.get(i)==1)
                tp=tp+1.0;
            else if(scores.get(i)==0 && trueLabel.get(i)==0)
                tn=tn+1.0;
            else if(scores.get(i)==1.0 && trueLabel.get(i)==0.0)
                fp=fp+1.0;
            else if(scores.get(i)==0.0 && trueLabel.get(i)==1.0)
                fn=fn+1.0;
        }
        
        precRec.add(tp/(tp+fp));
        precRec.add(tp/(tp+fn));
        
        return precRec;
    }
    
}
