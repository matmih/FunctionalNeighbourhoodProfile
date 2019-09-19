/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing prediction information (eukaryotic organisms)
 */
 
public class ClusPredictionInfo {
    
    public HashMap<String,ArrayList<Double>> classifierScores=new HashMap<>();//index in the array corresopnds to COG index in cgmap
    public HashMap<String,ArrayList<Double>> originalLabels=new HashMap<>();
    HashMap<String,Integer> cogToIndex=new HashMap<>();
    HashMap<Integer,String> indexToCog=new HashMap<>();
    ArrayList<String> gos = new ArrayList<>();
    
    public ClusPredictionInfo(){
        
    }
    
    public ClusPredictionInfo(ClusPredictionInfo c){
        
        Iterator<String> it = c.classifierScores.keySet().iterator();
        
        while(it.hasNext()){
            String k = it.next();
            
            ArrayList<Double> t = c.classifierScores.get(k);
            
            ArrayList<Double> tN = new ArrayList(t.size());
            
            for(Double v:t)
                tN.add(v);
            
            classifierScores.put(k, tN);          
        }
        
        it = c.originalLabels.keySet().iterator();
        
         while(it.hasNext()){
            String k = it.next();
            
            ArrayList<Double> t = c.originalLabels.get(k);
            
            ArrayList<Double> tN = new ArrayList(t.size());
            
            for(Double v:t)
                tN.add(v);
            
            originalLabels.put(k, tN);          
        }
         
         it = c.cogToIndex.keySet().iterator();
         
         while(it.hasNext()){
            String s = it.next();
            int v = c.cogToIndex.get(s);
            cogToIndex.put(s, v);
         }
         
         Iterator<Integer> it1 = c.indexToCog.keySet().iterator();
         while(it1.hasNext()){
                 int k = it1.next();
                 String v = c.indexToCog.get(k);
                 indexToCog.put(k, v);
         }
         
         this.gos.addAll(c.gos);
         
    }
    
     public void readPredictionInfoSCXval(OGGOMapping cgmap, String GO, File input){
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
     
     public void readPredictionInfoOOB(OGGOMapping cgmap, ARFF infile, File input){
          //classifierScores.get(GOs.get(i-2-cgmap.GOtoIndex.keySet().size())).set(cogToIndex.get(cog), score);
          //originalLabels.get(GOs.get(i-2-cgmap.GOtoIndex.keySet().size())).set(cogToIndex.get(cog), score1);
         
          HashMap<String, HashMap<String,Double>> data = new HashMap<>();
    
    Path p = Paths.get(input.getAbsolutePath());
    
    try{
        BufferedReader read = new BufferedReader(new FileReader(input));
               // Files.newBufferedReader(p);
        
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
           // System.out.println("data size: "+data.keySet().size());
        }
        read.close();
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
     //ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
    
   
       // FileWriter fw = new FileWriter("transformedPreds.pred");
        int numAnnots = 0;
        Iterator<String> it = data.keySet().iterator();
        int count = 0;
        ArrayList<String> goOrder = new ArrayList<>();
        int countC = 0;
        while(it.hasNext()){
            String cog = it.next();//.replaceAll("\"", "");
            //System.out.println("cog read: "+cog);
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
           // set(cogToIndex.get(cog), score1);
            }
            else{
                for(int i=0;i<goOrder.size();i++){
                    String go=goOrder.get(i);
                    predsN.add(preds.get(go));
                }
            }
            
           /* if(count==0){
                 for(int i=0;i<goOrder.size();i++){
                    fw.write(goOrder.get(i)+"\n");   
                }
            }*/
            
            int ind = infile.findCogIndex(cog);
           // System.out.println("cog index: "+ind);
            if(ind==-1)
                continue;
            
           String tar = infile.Targets.get(ind);
           String tarT[] = tar.split("@");
                  
                  HashSet<String> tarS = new HashSet<>();
                  
                  for(String t:tarT)
                      tarS.add(t.trim());
            
             // fw.write(cog+" ");
             for(int i=0;i<predsN.size();i++){
                 //System.out.println("p: "+predsN.get(i));
                 //System.out.println(tar.contains(goOrder.get(i)));
                  classifierScores.get(goOrder.get(i)).set(cogToIndex.get(cog), predsN.get(i));
                  String goPreds = goOrder.get(i).replace(":", "");
                  
                  
                  if(tarS.contains(goPreds)/*goOrder.get(i))*/){
                         originalLabels.get(goOrder.get(i)).set(cogToIndex.get(cog), 1.0);
                         numAnnots++;
                  }
                  else  originalLabels.get(goOrder.get(i)).set(cogToIndex.get(cog), 0.0);
                 
                 /*if((i+1)<predsN.size())
                     fw.write(predsN.get(i)+" ");
                 else
                      fw.write(predsN.get(i)+"\n");*/
            }
            
            count=1;
            
        }
        
        System.out.println("Number of annotations: "+numAnnots);
         
     }
    
     public void loadGOOrder(File input){
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
            int numAnnots=0, numAnnotsInt = 0;
            ArrayList<String> GOs=new ArrayList();
         try{
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    if(!func.contains("Original-p"))
                        System.out.println(line);
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
            }
            bufRdr1.close();
         }
            catch(IOException e){
                    e.printStackTrace();
                    }
     }
     
    public void readPredictionInfo(OGGOMapping cgmap, File input){
         BufferedReader bufRdr1=null; 
         
         try{
         bufRdr1= new BufferedReader(new FileReader(input));
         }
         catch(IOException e){
             e.printStackTrace();
         }
        
          String line;
            int data=0;
            
            
            System.out.println("map size cog: "+cgmap.CogGOmap.keySet().size());
             System.out.println("map size go: "+cgmap.GOtoIndex.keySet().size());
            
            int cogCount=0;
            int numAnnots=0, numAnnotsInt = 0;
            ArrayList<String> GOs=new ArrayList();
         try{
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    if(!func.contains("Original-p"))
                        System.out.println(line);
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
               // System.out.println(line);
                
                System.out.println("GOs data size: "+GOs.size());
                System.out.println("GOMap size: "+cgmap.GOtoIndex.keySet().size());
                
                if(data==1){
                String tmp[]=line.split(",");
                String cog=tmp[0].replace("\"", "");
                //System.out.println("cog: "+cog);
                //System.out.println("numFunctions: "+cgmap.GOtoIndex.keySet().size());
                if(!cogToIndex.containsKey(cog)){
                    cogToIndex.put(cog, cogCount++);
                    indexToCog.put(cogCount-1, cog);
                }
                //System.out.println("tmp size: "+tmp.length);
                ArrayList<Double> sc=new ArrayList<>(); 
                for(int i=/*cgmap.GOtoIndex.keySet().size()*/GOs.size()+2;i<2*GOs.size()/*2*cgmap.GOtoIndex.keySet().size()*/+2;i++){//1488-2975
                    
                   /* System.out.println("GO size: "+GOs.size());
                    System.out.println("index: "+(i-2-cgmap.GOtoIndex.keySet().size()));
                    System.out.println("GO: "+GOs.get(i-2-cgmap.GOtoIndex.keySet().size()));
                    System.out.println("cog: "+cog);*/
                    
                    if(!classifierScores.containsKey(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            classifierScores.put(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/), tmpA);
                    }
                    
                    if(!originalLabels.containsKey(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/))){
                            ArrayList<Double> tmpA=new ArrayList<>(Collections.nCopies(cgmap.CogGOmap.keySet().size(), 0.0));
                            originalLabels.put(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/), tmpA);
                    }
                    
                    double score=Double.parseDouble(tmp[i]);
                    double score1=Double.parseDouble(tmp[i-GOs.size()/*cgmap.GOtoIndex.size()*/]);
                    //System.out.println("score, score1: "+score+" "+score1);
                    classifierScores.get(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/)).set(cogToIndex.get(cog), score);
                    originalLabels.get(GOs.get(i-2-GOs.size()/*cgmap.GOtoIndex.keySet().size()*/)).set(cogToIndex.get(cog), score1);
                    
                    if(score1 == 1)
                        numAnnots++;
                    
                    int score2 = Integer.parseInt(tmp[i-GOs.size()]);
                    
                    if(score2!=0 && score2!=1){
                        System.out.println("SW index: "+i);
                    }
                    
                    if(score2 == 1)
                        numAnnotsInt++;
                    
                   /* if(score1!= 0 && score!= 1){
                        System.out.println("Index: "+i);
                        System.out.println("Score: "+score1);
                    }*/
                    //sc.add(score);
                }
               /* cogFunc.put(cog, sc);
                cogs.add(cog);*/
                }
            }
            bufRdr1.close();
            
            System.out.println("GO functions read from file: "+GOs.size());
            System.out.println("OGs read from file: "+cogToIndex.keySet().size());
          System.out.println("Number of annotations: "+numAnnots);  
           System.out.println("Number of annotations int: "+numAnnotsInt); 
         }
         catch(IOException e){
             e.printStackTrace();
         }       
    }
    
    public ArrayList<Double> sortedScores(ArrayList<Double> scores){
        HashSet<Double> tmpS = new HashSet<>();
        tmpS.addAll(scores);
        ArrayList<Double> tmp=new ArrayList(tmpS);
        
        Collections.sort(tmp);
        
        return tmp;  
    }
    
    
     public HashSet<String> getGenes(String Og, HashMap<String,HashSet<Pair<String,String>>> geneOGMap){
        
        Iterator<String> geneIt= geneOGMap.keySet().iterator();
        HashSet<String> associatedGenes=new HashSet<>();
        
        while(geneIt.hasNext()){
            String gene=geneIt.next();
            HashSet<Pair<String,String>> ogs=geneOGMap.get(gene);
            for(Pair<String,String> p:ogs){
            if(p.getValue0().equals(Og)){
                      associatedGenes.add(gene);
                      break;
                }
            }
        }
        return associatedGenes;
    }
    
    public ArrayList<Integer> numberOfNewGenes(ArrayList<Double> predicted, ArrayList<Double> trueClass, HashMap<String,HashSet<Pair<String,String>>> geneOGMap){
        int oGenes=0, nGenes=0, difference=0;
        HashSet<String> trueGenes=new HashSet<>();
        HashSet<String> predictedGenes=new HashSet<>();
        HashSet<String> newlypredictedGenes=new HashSet<>();
        /*System.out.println("true class size: "+trueClass.size());
        
        System.out.println("Labs in NOfNewGenes: ");
                            for(int i=0;i<trueClass.size();i++)
                                if(trueClass.get(i)!=0)
                                    System.out.print(trueClass.get(i)+" ");
                            System.out.println();*/
        
        for(int i=0;i<trueClass.size();i++){
            if(trueClass.get(i)==1.0){
                //System.out.println("Before adding");
               // System.out.println("Cog: "+indexToCog.get(i));
               // System.out.println("Set: "+geneOGMap.getGenes(indexToCog.get(i)).toString());
                HashSet<String> st=getGenes(indexToCog.get(i),geneOGMap);//geneOGMap.getGenes(indexToCog.get(i));
               // System.out.println("Adding element");
                trueGenes.addAll(st);
            }
            if(predicted.get(i)==1.0){
                HashSet<String> st=getGenes(indexToCog.get(i),geneOGMap);//geneOGMap.getGenes(indexToCog.get(i));
                predictedGenes.addAll(st);
            }
            
            if(trueClass.get(i)==0.0 && predicted.get(i)==1.0){
                HashSet<String> st=getGenes(indexToCog.get(i),geneOGMap);//geneOGMap.getGenes(indexToCog.get(i));
                newlypredictedGenes.addAll(st);
            }
            
        }
        //System.out.println("True genes: "+trueGenes.size());
        difference=predictedGenes.size()-trueGenes.size();
        ArrayList<Integer> res=new ArrayList<>();
        System.out.println("trueGenes size: "+trueGenes.size());
        res.add(trueGenes.size()); res.add(predictedGenes.size()); res.add(difference); res.add(newlypredictedGenes.size());
        
        return res;
    }
    
    
     public ArrayList<Integer> numberOfNewGenesOrganism(ArrayList<Double> predicted, ArrayList<Double> trueClass, GenesAndOgsOrganism organismMapping, HashSet<String> npgenes){
        int oGenes=0, nGenes=0, difference=0;
        HashSet<String> trueGenes=new HashSet<>();
        HashSet<String> predictedGenes=new HashSet<>();
        HashSet<String> newlypredictedGenes=new HashSet<>();
        /*System.out.println("true class size: "+trueClass.size());
        
        System.out.println("Labs in NOfNewGenes: ");
                            for(int i=0;i<trueClass.size();i++)
                                if(trueClass.get(i)!=0)
                                    System.out.print(trueClass.get(i)+" ");
                            System.out.println();*/
        
        System.out.println("predicted size: "+predicted.size());
        System.out.println("trueClass size: "+trueClass.size());
        
        System.out.println("mapping size og: "+organismMapping.OGGeneMappingOrganism.size());
        System.out.println("mapping size org: "+organismMapping.geneOGMappingOrganism.size());
      
        for(int i=0;i<trueClass.size();i++){
            if(!indexToCog.containsKey(i)){
               // System.out.println("Cog index not contained! Warning!");
                continue;
            }
            
             if(!organismMapping.OGGeneMappingOrganism.containsKey(indexToCog.get(i)))
                    continue;
             
            if(trueClass.get(i)==1.0){
               
                //System.out.println("Before adding");
               // System.out.println("Cog: "+indexToCog.get(i));
               // System.out.println("Set: "+geneOGMap.getGenes(indexToCog.get(i)).toString());
                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));
               // System.out.println(indexToCog.get(i));
               // System.out.println("contained in mapping: "+organismMapping.OGGeneMappingOrganism.containsKey(indexToCog.get(i)));
               // System.out.println("Num genes true: "+st.size());
               // System.out.println("Adding element");
                trueGenes.addAll(st);
            }
            if(predicted.get(i)==1.0){
                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));
               // System.out.println("Num genes pred: "+st.size());
                predictedGenes.addAll(st);
            }
            
            if(trueClass.get(i)==0.0 && predicted.get(i)==1.0){
                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));
               // System.out.println("Num genes new: "+st.size());
                newlypredictedGenes.addAll(st);
            }
            
        }
        //System.out.println("True genes: "+trueGenes.size());
        difference=predictedGenes.size()-trueGenes.size();
        ArrayList<Integer> res=new ArrayList<>();
       // System.out.println("trueGenes size: "+trueGenes.size());
        res.add(trueGenes.size()); res.add(predictedGenes.size()); res.add(difference); res.add(newlypredictedGenes.size());
        npgenes.addAll(newlypredictedGenes);
        
        return res;
    }
     
     public ArrayList<Integer> numberOfNewGenesOrganismDist(int taxID,String go, ArrayList<Double> predicted, ArrayList<Double> trueClass, GenesAndOgsOrganism organismMapping, Mappings map ,HashSet<String> npgenes){
        int oGenes=0, nGenes=0, difference=0;
        HashSet<String> trueGenes=new HashSet<>();
        HashSet<String> predictedGenes=new HashSet<>();
        HashSet<String> newlypredictedGenes=new HashSet<>();

       // System.out.println("predicted size: "+predicted.size());
       // System.out.println("trueClass size: "+trueClass.size());
        
        //System.out.println("mapping size og: "+organismMapping.OGGeneMappingOrganism.size());
       // System.out.println("mapping size org: "+organismMapping.geneOGMappingOrganism.size());
      
        for(int i=0;i<trueClass.size();i++){
            if(!indexToCog.containsKey(i)){
               // System.out.println("IndexToCog not contains index: "+i);
                continue;
            }
            
             if(!organismMapping.OGGeneMappingOrganism.containsKey(indexToCog.get(i))){
                 //System.out.println("Organism mapping does not contain: "+indexToCog.get(i));
                    continue;
             }
             
            if(trueClass.get(i)==1.0){
               

                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));
                
                for(String gp : st){
                    int contained =0;
                   
                    if(map.mappings.containsKey(gp)){
                         HashSet<Pair<String,String>> cand = map.mappings.get(gp);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if(trueGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }      
                    }
                  }
                    
                    for(String t:trueGenes){
                        if(!map.mappings.containsKey(t))
                            continue;
                         HashSet<Pair<String,String>> cand = map.mappings.get(t);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if(trueGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }     
                    }
                    if(contained == 1)
                        break;
                   }
                    
                    if(contained==0)
                        trueGenes.add(gp);
                    
                }
               // System.out.println("true genes: "+trueGenes.size()+" "+go);
               // trueGenes.addAll(st);
            }
            if(predicted.get(i)==1.0){
                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));

                for(String gp : st){
                    int contained =0;
               
                    if(map.mappings.containsKey(gp)){
     HashSet<Pair<String,String>> cand = map.mappings.get(gp);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if( predictedGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }      
                    }
                    }  
                    
                     for(String t:predictedGenes){
                         if(!map.mappings.containsKey(t))
                             continue;
                         HashSet<Pair<String,String>> cand = map.mappings.get(t);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if(predictedGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }     
                    }
                    if(contained == 1)
                        break;
                   }
                    
                    if(contained==0){
                         predictedGenes.add(gp);
                        // System.out.println("Predicted: "+go+" "+gp);
                    }
                    
                }
               // System.out.println("predicted genes: "+predictedGenes.size());
                //predictedGenes.addAll(st);
            }
            
            if(trueClass.get(i)==0.0 && predicted.get(i)==1.0){
                HashSet<String> st=organismMapping.OGGeneMappingOrganism.get(indexToCog.get(i));//geneOGMap.getGenes(indexToCog.get(i));

                for(String gp : st){
                    int contained =0;
                    
                    if(map.mappings.containsKey(gp)){
                        HashSet<Pair<String,String>> cand = map.mappings.get(gp);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if(newlypredictedGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }      
                    }
                    
                     for(String t:newlypredictedGenes){
                         if(!map.mappings.containsKey(t))
                             continue;
                          cand = map.mappings.get(t);
                    for(Pair<String,String> p:cand){
                        if(p.getValue1().equals(gp)){
                            if(newlypredictedGenes.contains(p.getValue1()))
                                contained=1;
                            break;
                        }     
                    }
                    if(contained == 1)
                        break;
                   }
                    
                    }   
                    if(contained==0){
                         newlypredictedGenes.add(gp);
                        // System.out.println("Predicted: ");
                         //System.out.println(gp+" "+organismMapping.geneOGMappingOrganism.get(gp)+" "+go);
                    }                    
                }
                
                //newlypredictedGenes.addAll(st);
            }
            
        }

        difference=predictedGenes.size()-trueGenes.size();
        ArrayList<Integer> res=new ArrayList<>();

        res.add(trueGenes.size()); res.add(predictedGenes.size()); res.add(difference); res.add(newlypredictedGenes.size());
        npgenes.addAll(newlypredictedGenes);
        
        return res;
    }
    
    public ArrayList<Double> predictedClass(ArrayList<Double> scores, double threshold){
        
        ArrayList<Double> tmp=new ArrayList<>(scores.size());
        
        for(int i=0;i<scores.size();i++)
            if(scores.get(i)>=threshold)
                tmp.add(1.0);
            else tmp.add(0.0);
        
        return tmp;
    }
    
     public void predictedClass1(ArrayList<Double> scores, ArrayList<Double> tmp, double threshold){
        
        for(int i=0;i<scores.size();i++)
            if(scores.get(i)>=threshold)
                tmp.set(i,1.0);
            else tmp.set(i,0.0);
        
        //return tmp;
    }
    
    public ArrayList<Double> precRec(ArrayList<Double> scores, ArrayList<Double> trueLabel){
        
        ArrayList<Double> precRec=new ArrayList<>(2);
        double tp=0.0,fp=0.0,tn=0.0,fn=0.0;
       
        try{
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
        }
        catch(Exception e){
            e.printStackTrace();
            System.out.println("Something is wrong!");
        }
        precRec.add(tp/(tp+fp));
        precRec.add(tp/(tp+fn));
        
        return precRec;
    }
    
    public void precRec1(ArrayList<Double> scores, ArrayList<Double> trueLabel, ArrayList<Double> precRec){
        
        double tp=0.0,fp=0.0,tn=0.0,fn=0.0;
       
       // try{
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
        /*}
        catch(Exception e){
            e.printStackTrace();
            System.out.println("Something is wrong!");
        }*/
        precRec.set(0,tp/(tp+fp));
        precRec.set(1,tp/(tp+fn));

    }
    
}
