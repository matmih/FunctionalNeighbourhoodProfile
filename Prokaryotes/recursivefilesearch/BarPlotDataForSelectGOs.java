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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import semanticsimilarity.Graph;
import semanticsimilarity.Path;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create data to plot a comparative LOR bar plot for selected GOs (displaying share of (in)significantly enriched GO functions with the selected GO)
 */
public class BarPlotDataForSelectGOs {


    static public void main(String arg[]){
        
           File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File OgGOs30=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_30_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File semFile=new File("PairwiseGOSim.txt");
         
         String redCog="COGRed150.txt"; String redGO="GORed5.txt";
        
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
             
         semanticsimilarity.GOMap mapFull=new semanticsimilarity.GOMap();
         mapFull.CreateGOMap(cgmap);
         System.out.println("Num GOs: "+cgmap.GOtoIndex.keySet().size());
         System.out.println("GO on index 0 "+cgmap.IndexToGO.get(0));
         GOSemanticSimilarity goSem=new GOSemanticSimilarity(mapFull.GOmap.keySet().size());
         goSem.readSimilarity(semFile);
        
          ArrayList<HashSet<String>> categories=new ArrayList<>();
        GODevider gd=new GODevider();
        try{
            categories=gd.ExtractCategories(new File("go_201401-termdb.obo-xml.gz"),cgmap);
        }
        catch(Exception e){
            e.printStackTrace();
        }  
         
        semanticsimilarity.GOMap mapMod=new semanticsimilarity.GOMap();
        File fileWithUniprotFrequencyCounts=new File("Uniprot-freqs-from idmapping-2014-06-11.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
        semanticsimilarity.GOMap mapFull1=new semanticsimilarity.GOMap();
        mapFull1.CreateGOMap(cgmap);
        mapFull1.printGOMap();
        mapFull1.loadFrequencies(fileWithUniprotFrequencyCounts);

          String path="";

        ArrayList<String> functions=new ArrayList<>();      
       HashMap<Integer,ArrayList<Double>> pvals=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> pvalsCT=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> pvalsDiff=new HashMap<>();
          
        for(int iter=0;iter<cgmap.GOtoIndex.keySet().size();iter++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            ArrayList<Double> tmp1=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            ArrayList<Double> tmp2=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvals.put(iter, tmp);
            pvalsDiff.put(iter, tmp1);
            pvalsCT.put(iter,tmp2);
        }
       
      File pvalsA = new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsAll.txt");  
       File fdrPvals=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORZeroFDR.txt");
       File fdrPvalsDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORDiffFDR.txt");
         BufferedReader bufRdr=null;
         BufferedReader bufRdrDiff=null;
         
         
            try 
        {
            bufRdr = new BufferedReader(new FileReader(pvalsA));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double pv=Double.parseDouble(line.trim());
                pvalsCT.get(count).set(count1, pv);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
          }
           catch(Exception e){
               e.printStackTrace();
           }
         
         
           try 
        {
            bufRdr = new BufferedReader(new FileReader(fdrPvals));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrCorr=Double.parseDouble(line.trim());
                pvals.get(count).set(count1, fdrCorr);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
          }
           catch(Exception e){
               e.printStackTrace();
           }
           
           
            try 
        {
            bufRdrDiff = new BufferedReader(new FileReader(fdrPvalsDiff));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrCorr=Double.parseDouble(line.trim());
                pvalsDiff.get(count).set(count1, fdrCorr);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
          }
           catch(Exception e){
               e.printStackTrace();
           }
 
        functions.add("GO0005975"); functions.add("GO0006260"); functions.add("GO0006974");
        functions.add("GO0016051"); functions.add("GO0006457"); functions.add("GO0046700");
        functions.add("GO0006310"); functions.add("GO0046903"); functions.add("GO0008610");
        functions.add("GO0006865"); functions.add("GO0006974"); functions.add("GO0006508");
        
          System.out.println("Adjacency graph created!");
          
          HashMap<String,HashSet<String>> GOParents=new HashMap<>();
        
          boolean writePVals=false;
          
        for(int i=0;i<functions.size();i++){
            ArrayList<Path> Paths=new ArrayList<Path>();
        semanticsimilarity.Path visited=new semanticsimilarity.Path();
        visited.path.add(mapFull1.GOmap.get(functions.get(i)));
        HashSet<String> gos=new HashSet<>();
        gos.add(functions.get(i));
        
        for(int i1=0;i1<Paths.size();i1++){
             for(int j=0;j<Paths.get(i1).path.size();j++)
            gos.add(mapFull1.GOmapNumeric.get(Paths.get(i1).path.get(j)));
        }
        System.out.println("Function: "+i+" path computed...");
        GOParents.put(functions.get(i), gos);
      }
        
        File prokaryotInput=new File("prokaryotic.txt");
        HashSet<String> pkGOs=new HashSet<>();
        
          try (BufferedReader bufRdr1 = new BufferedReader(new FileReader(prokaryotInput)))
        {
            String line;
            
            while ((line = bufRdr1.readLine()) != null)
            {
                line = line.trim();
                
                pkGOs.add(line);
                            
            }
        }
       catch(Exception e){
           e.printStackTrace();
       }
          
          HashMap<String,ArrayList<String>> sigPairs=new HashMap<>();
          
          HashMap<String,ArrayList<Double>> medClHistSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHistSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHistSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> medClHistSigSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHistSigSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHistSigSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> medClHistSigNSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHistSigNSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHistSigNSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> medClHistNSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHistNSig=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHistNSig=new HashMap<>();
          
        for(int iter=0;iter<functions.size();iter++){ 
            medClHistSig.put(functions.get(iter),new ArrayList<Double>());
            clParHistSig.put(functions.get(iter), new ArrayList<Double>());
            distHistSig.put(functions.get(iter), new ArrayList<Double>());
             medClHistSigSig.put(functions.get(iter),new ArrayList<Double>());
            clParHistSigSig.put(functions.get(iter), new ArrayList<Double>());
            distHistSigSig.put(functions.get(iter), new ArrayList<Double>());
            medClHistNSig.put(functions.get(iter),new ArrayList<Double>());
            clParHistNSig.put(functions.get(iter), new ArrayList<Double>());
            distHistNSig.put(functions.get(iter), new ArrayList<Double>());
            medClHistSigNSig.put(functions.get(iter),new ArrayList<Double>());
            clParHistSigNSig.put(functions.get(iter), new ArrayList<Double>());
            distHistSigNSig.put(functions.get(iter), new ArrayList<Double>());
            
            if(!pkGOs.contains(functions.get(iter)))
                continue;
            
     LoadOR lor=new LoadOR();
     File ORInp=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+functions.get(iter)+".txt");//cgmap.IndexToGO.get(iter)
     lor.loadContValue(ORInp);
     
    ArrayList<Double> ctOrig=lor.cont.get(cgmap.IndexToGO.get(iter));
        Double orOrig=((ctOrig.get(0)+0.5)/(ctOrig.get(1)+0.5))/((ctOrig.get(2)+0.5)/(ctOrig.get(3)+0.5));
        Double logOrOrig=Math.log10(orOrig)/Math.log10(2); 

        int aInd=mapFull.GOmap.get(functions.get(iter));       
        
    for(int goindex=0;goindex<lor.cont.keySet().size();goindex++){
        String GOname=cgmap.IndexToGO.get(goindex);
        int iGO=mapFull.GOmap.get(GOname);
        if(!pkGOs.contains(GOname) || !categories.get(0).contains(GOname))
            continue;
        
            ArrayList<Double> ct=lor.cont.get(GOname);
        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
        Double logOr=Math.log10(or)/Math.log10(2);

        if(goSem.similarity[iGO-1][aInd-1]>2 && !GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
            if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && pvalsDiff.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>logOrOrig && logOrOrig>0 && pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05){
                if(!(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0)&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05)
                     medClHistSigNSig.get(functions.get(iter)).add(logOr);
                else
                    medClHistSigSig.get(functions.get(iter)).add(logOr);        
            }
            else if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05)
                 medClHistSig.get(functions.get(iter)).add(logOr);
            else  medClHistNSig.get(functions.get(iter)).add(logOr);
        }
        
        if(/*goSem.similarity[iGO-1][aInd-1]>2 && !*/GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
            if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && pvalsDiff.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>logOrOrig && logOrOrig>0&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05){
                if(!(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0)&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05)
                     clParHistSigNSig.get(functions.get(iter)).add(logOr);
                else
                    clParHistSigSig.get(functions.get(iter)).add(logOr);
            }
            else if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0 && pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05)
                 clParHistSig.get(functions.get(iter)).add(logOr);
            else  clParHistNSig.get(functions.get(iter)).add(logOr);
        }
        
        if(goSem.similarity[iGO-1][aInd-1]<=2 && !GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
           if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && pvalsDiff.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>logOrOrig && logOrOrig>0&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05){
               if(!(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0)&& pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05){
                   System.out.println(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)+" "+pvalsDiff.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)+" "+logOr+" "+logOrOrig);
                     distHistSigNSig.get(functions.get(iter)).add(logOr);
               }
               else
                    distHistSigSig.get(functions.get(iter)).add(logOr);
           }
            else if(pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1 && logOr>0 && pvalsCT.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.05)
                 distHistSig.get(functions.get(iter)).add(logOr);
            else  distHistNSig.get(functions.get(iter)).add(logOr);
        }
    }
    System.out.println("Iter complited: "+iter);
     }
       
       String fileName="DensityOverlayBySignif";
       
       for(int i=0;i<functions.size();i++){
           String fileNameTmp=fileName+functions.get(i)+".txt";
           try{
         FileWriter fw = new FileWriter(fileNameTmp);                    
         
         fw.write("logOR \t category \t hist\n");
         
          ArrayList<Double> tmp=medClHistSigSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"close \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LORD\n");
          }
          
          tmp=medClHistSigNSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"close \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LORD\n");
          }
          
           tmp=medClHistSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"close \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LORD\n");
          }
          
          tmp=medClHistNSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"close \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LORD\n");
          }
          
          
          tmp=distHistSigSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"distant \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LORD\n");
          }
          
          tmp=distHistSigNSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"distant \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LORD\n");
          }
          
           tmp=distHistSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"distant \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"significant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LORD\n");
          }
          
          tmp=distHistNSig.get(functions.get(i));
          
          for(int ind=0;ind<tmp.size();ind++){
              fw.write(tmp.get(ind)+" \t "+"distant \t "+"D\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LOR\n");
              fw.write(tmp.get(ind)+" \t "+"unsignificant \t "+"LORD\n");
          }
         
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}
       }         
       }
        
    }
