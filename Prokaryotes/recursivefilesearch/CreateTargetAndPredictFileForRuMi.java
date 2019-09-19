/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create target and prediction files for the Ru-MI curves
 */
public class CreateTargetAndPredictFileForRuMi {
    
    static public void main(String args[]){
        File input=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\genome5k=10.train.1.pred.arff");
         
           File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         
         String redCog="COGRed150.txt"; String redGO="GORed5.txt";

        GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(OgGOs);
        
        ReducedGOTranslator rgt=new ReducedGOTranslator();
        rgt.ReadAndTranslate(new File(redGO));

    OntologyTools.GeneOntology myGO=null;

    try{
    myGO= new GeneOntology("go_201401-termdb.obo-xml");
    }
    catch(IOException e){
        e.printStackTrace();
    }

    HashSet<Integer> extendedGOs=new HashSet<>();

    for(int go:rgt.goI){
        Collection<GOTerm> parents = myGO.get(go).getAllParents();
    for ( GOTerm curPar : parents ) {
         if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
           extendedGOs.add(curPar.getId());
    }

    extendedGOs.add(go);

    }

    rgt.translate(extendedGOs);

        ReducedOgsGOs red=new ReducedOgsGOs();
        red.LoadReducedOGGos(new File(redCog), new File(redGO));
        red.redGOs=rgt.goS;
        System.out.println("Loading complete");
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);
         
         
         ArrayList<String> cogsOrdered=new ArrayList<>();
         HashMap<String,Pair<Double,Double>> scores=new HashMap<>();
         
         ArrayList<String> cogs=new ArrayList<>();
         ArrayList<String> GOs=new ArrayList<>();
         HashMap<String,ArrayList<Double>> cogFunc=new HashMap<>();
         
         int FileType=0;
         
        if(FileType==0){ 
         try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(input)))
        {
            String line;
            int count=0;
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    func=func.replace("Original-p-", "");
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    GOs.add(func);
                    continue;
                }
                else if(line.contains("@ATTRIBUTE"))
                    continue;
                else if(line.contains("@DATA"))
                    continue;
                else if(line.contains("RELATION"))
                    continue;
                
                if(line.equals(""))
                    continue;
                String tmp[]=line.split(",");
                String cog=tmp[0].replace("\"", "");
                System.out.println("tmp size: "+tmp.length);
                ArrayList<Double> sc=new ArrayList<>();
                for(int i=cgmap.GOtoIndex.keySet().size()+2;i<2*cgmap.GOtoIndex.keySet().size()+2;i++){//1488-2975
                    double score=Double.parseDouble(tmp[i]);
                    sc.add(score);
                }
                cogFunc.put(cog, sc);
                cogs.add(cog);
            }
            bufRdr1.close();
        }
       catch(Exception e){
           e.printStackTrace();
       }
     }
 
        if(FileType==1){
        
            try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(input)))
        {
            String line;
            int count=0;
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                if(line.contains("@ATTRIBUTE") && line.contains("numeric")){
                    String tmp[]=line.split(" ");
                    String func=tmp[1];
                    func=func.replace("Original-p-", "");
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    GOs.add(func);
                    continue;
                }
                else if(line.contains("@ATTRIBUTE"))
                    continue;
                else if(line.contains("@DATA")){
                    break;
                }
                else if(line.contains("RELATION"))
                    continue;
                
                if(line.equals(""))
                    continue;
                String tmp[]=line.split(",");
                String cog=tmp[0].replace("\"", "");
                System.out.println("tmp size: "+tmp.length);
                ArrayList<Double> sc=new ArrayList<>(); 
                for(int i=cgmap.GOtoIndex.size();i<2*cgmap.GOtoIndex.size()+1;i++){//1488-2975
                    double score=Double.parseDouble(tmp[i]);
                    sc.add(score);
                }
                cogFunc.put(cog, sc);
                cogs.add(cog);
            }
            bufRdr1.close();
            
            input=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\DistanceLocD_oob.preds");
            BufferedReader bufRdr2 = null;
            
            try{
                bufRdr2=new BufferedReader(new FileReader(input));
            }
            catch(IOException e){
                e.printStackTrace();
            }
            
            while ((line = bufRdr2.readLine()) != null)
            {
                line = line.trim();
                String tmp[]=line.split("\t");
                String cog=tmp[0].replace("\"", "");
                System.out.println("tmp size: "+tmp.length);
                ArrayList<Double> sc=new ArrayList<>();
                String scs[]=tmp[2].split(";");
                for(int i=0;i<cgmap.GOtoIndex.size();i++){//1487
                    double score=Double.parseDouble(scs[i]);
                    sc.add(score);
                }
                cogFunc.put(cog, sc);
                cogs.add(cog);
            }
            bufRdr2.close();
            
        }
       catch(Exception e){
           e.printStackTrace();
       }
        }
        
        
        if(FileType==2){
        
            try{
               File input1=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\AssociationLogDFullAnot.txt");
               BufferedReader bufRdr2 = new BufferedReader(new FileReader(input1));
               String line;
               int count=0;
               
               while((line=bufRdr2.readLine())!=null){
                   if(count==0){
                      count=1;
                       continue;
                   }
                   if(line.equals(""))
                       continue;
                   line.trim();
                   String tmp[]=line.split("\t");
                   cogs.add(tmp[0]);
               }
               bufRdr2.close();
            }
            catch(IOException e){
                e.printStackTrace();
            }
            
            
            try
        {
            String line;
            int count=0;
        
            input=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\GFPScoreMatrix.txt");
            BufferedReader bufRdr2 = null;
            
            try{
                bufRdr2=new BufferedReader(new FileReader(input));
            }
            catch(IOException e){
                e.printStackTrace();
            }
            
            double mat[][]=new double[cgmap.GOtoIndex.keySet().size()][cogs.size()];
            
            while ((line = bufRdr2.readLine()) != null)
            {
                line = line.trim();
                String tmp[]=line.split("\t");
                System.out.println("tmp size: "+tmp.length);
                for(int i=0;i<cogs.size();i++){
                    double score=Double.parseDouble(tmp[i]);
                    mat[count][i]=score;
                }
                count++;
            }
            bufRdr2.close();
            
            for(int i=0;i<cogs.size();i++){
                ArrayList<Double> sc=new ArrayList<>();
                String cog=cogs.get(i);
                for(int j=0;j<cgmap.GOtoIndex.keySet().size();j++){
                    sc.add((mat[j][i]+1)/2);
                }
                cogFunc.put(cog, sc);
            }
            
             for(int i=0;i<cgmap.IndexToGO.keySet().size();i++){
                  String func=cgmap.IndexToGO.get(i);
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    GOs.add(func);
             }
            
        }
       catch(Exception e){
           e.printStackTrace();
       }
        }
        
        if(FileType == 3){
            ClusPredictionInfo predInfo=new ClusPredictionInfo();
             File inputData = new File("BaselineOrganism1k=4.arff");
         File predInput=new File("BaselineOGsk=4_oob.preds");
              ARFF d = new ARFF();
            d.loadARFFTarget(inputData);
               
            predInfo.readPredictionInfoOOB(cgmap, d, predInput);
            
            for(int i=0;i<cgmap.IndexToGO.keySet().size();i++){
                  String func=cgmap.IndexToGO.get(i);
                    func=func.replace("GO", "");
                    func="GO:"+func;
                    GOs.add(func);
             }
            
            System.out.println("PD cog ind size: "+predInfo.cogToIndex.keySet().size());
            
             for(int i=0;i<predInfo.cogToIndex.keySet().size();i++){
                 ArrayList<Double> sc = new ArrayList<>();
                    String cog = predInfo.indexToCog.get(i);
                    for(int j=0;j<cgmap.IndexToGO.keySet().size();j++){
                        String go = cgmap.IndexToGO.get(j);
                        go= go.replace("GO", "");
                        go = "GO:"+go;
                        if(!predInfo.classifierScores.containsKey(go))
                            continue;
                        sc.add(predInfo.classifierScores.get(go).get(predInfo.cogToIndex.get(cog)));
                    }
                    System.out.println("sc: "+i+" "+sc.size());
                    cogFunc.put(cog, sc);
                    
                }
             
             System.out.println("CogFunc size: "+cogFunc.keySet().size());
        }
        
        
        File output=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PredsCurve.txt");
         
         try
            {
                FileWriter fw = new FileWriter(output);
                
                for(String cog:cogFunc.keySet()){
                    ArrayList<Double> tmp=cogFunc.get(cog);
                    
                    for(int i=0;i<tmp.size();i++){
                        fw.write(cog+"\t"+GOs.get(i)+"\t"+tmp.get(i)+"\n");
                    }
                }
                
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            }  
       
         int generateLabels=1;
         
         if(generateLabels==1){
         output=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\LabelsCurve.txt");
  
         try
            {
                FileWriter fw = new FileWriter(output);
                
                for(String cog:cgmap.CogGOmap.keySet()){
                    ArrayList<String> tmp=cgmap.CogGOmap.get(cog);
                    
                    for(int i=0;i<tmp.size();i++){
                        String g="GO:"+tmp.get(i).replace("GO", "");
                        fw.write(cog+"\t"+g+"\t"+1.0+"\n");
                    }
                }
                
                fw.close();
           
            }
           catch(Exception e){
                e.printStackTrace();
            } 
         }
       }     
    }
