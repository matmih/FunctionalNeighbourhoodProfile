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
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import semanticsimilarity.*;
/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to assess predictive importance of features derived from functional neighbourhoods of different GO terms
  *divided by semantical similarity to the target GO function. The number of created attributes is equal in all produced datasets
  *attributes corresponding to functions that are not at the target similarity level to the target function are replaced with randomized attributes
 */
public class SingleClassPerformanceTest {
     
    public static void main(String[] args) {
        
         File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
         File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File semFile=new File("PairwiseGOSim.txt");
         
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
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);
        
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
		
        ArrayList<String> functions=new ArrayList<>();
        
        functions.add("GO0005975"); functions.add("GO0006260"); functions.add("GO0006974");
        functions.add("GO0016051"); functions.add("GO0006457"); functions.add("GO0046700");
        functions.add("GO0006310"); functions.add("GO0046903"); functions.add("GO0008610");
        functions.add("GO0006865"); functions.add("GO0006508");
     
        HashMap<String,HashSet<String>> GOParents=new HashMap<>();
  
        for(int i=0;i<functions.size();i++){
        HashSet<String> gos=new HashSet<>();
        gos.add(functions.get(i));
        
        int go1=rgt.translateToIntSingle(functions.get(i));
        HashSet<Integer> pars=new HashSet<>();
        Collection<GOTerm> parents = myGO.get(go1).getAllParents();
    for ( GOTerm curPar : parents ) {
         if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
           pars.add(curPar.getId());
    }
    
    HashSet<String> parsString=rgt.translateToString(pars);
    gos.addAll(parsString);
        
        GOParents.put(functions.get(i), gos);
      }

         String path="";
         
         File input=new File(path+"BaselineOrganism1.arff");
        
        Classify cl=new Classify();
        FileDeleter del=new FileDeleter();
        Random rand = new Random();
        int runInd=0;
        ArrayList<HashMap<Integer,ArrayList<Double>>> results=new ArrayList<>();
        
        for(int i=0;i<functions.size();i++){
            HashMap<Integer,ArrayList<Double>> tmp=new HashMap<>();
            results.add(tmp);
        }            
        
        File FdrPvalsGZeroLOR=new File("pValsLORZeroFDR.txt");
         
          HashMap<Integer,ArrayList<Double>> fdrPvalsGZeroLOR=new HashMap<>();   
       
        for(int iter1=0;iter1<cgmap.GOtoIndex.keySet().size();iter1++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            fdrPvalsGZeroLOR.put(iter1, tmp);
        }
        
        BufferedReader bufRdr=null;
        
          try 
        {
            bufRdr = new BufferedReader(new FileReader(FdrPvalsGZeroLOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrpvalsGZeroLOR=Double.parseDouble(line.trim());
                fdrPvalsGZeroLOR.get(count).set(count1, fdrpvalsGZeroLOR);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
            bufRdr.close();
          }
           catch(Exception e){
               e.printStackTrace();
           }
        
         File pValsOR=new File("PvalsAll.txt");
         
          HashMap<Integer,ArrayList<Double>> pvalsOR=new HashMap<>();   
       
        for(int iter1=0;iter1<cgmap.GOtoIndex.keySet().size();iter1++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvalsOR.put(iter1, tmp);
        }
        
         bufRdr=null;
        
          try 
        {
            bufRdr = new BufferedReader(new FileReader(pValsOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double pOR=Double.parseDouble(line.trim());
                pvalsOR.get(count).set(count1, pOR);
                count1++;
                if((count1)%cgmap.GOtoIndex.keySet().size()==0){
                    count1=0;
                    count++;
                }
                    
             }
            bufRdr.close();
          }
           catch(Exception e){
               e.printStackTrace();
           }
          
          
          LoadOR lor=new LoadOR();
          
    for( runInd=0;runInd<200;runInd++){
         int seed=rand.nextInt((100 - 1) + 1) + 1;
          System.out.println("Iterations... "+runInd);
      for(int iter=0;iter<functions.size();iter++){  
         
           ARFF f=new ARFF();
          String GO=functions.get(iter);
         OddsRatio or=new OddsRatio();
         
         File ORInp=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContAll/"+"ContingencyOut"+GO+".txt");
         lor.loadContValue(ORInp);

         f.loadARFF(input);
         f.assignLabels(GO);
         ArrayList<Double> maxs=new ArrayList<>();
         ArrayList<Double> mins=new ArrayList<>();
         f.attrStats(maxs, mins);
         
         f.writeARFFPO(path+"SingleFunction"+GO+"BPP.arff",categories,false);
         cl.classify(path+"SingleFunction"+GO+"BPP.arff", path+"SingleFunction"+GO+"Ind:"+(runInd+1)+"BPP.out",rand,runInd,iter,1,seed,results);
       
          f.writeARFFPO(path+"SingleFunction"+GO+"ClNClBPPar.arff", GO, GOParents,goSem,mapFull ,cgmap,or,categories,Double.POSITIVE_INFINITY,20, false,false,true,maxs,mins);///CL/Par OK!
          cl.classify(path+"SingleFunction"+GO+"ClNClBPPar.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"ClNClBPPar.out",rand,runInd,iter,3,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"ClNClBPPar.arff"));
          
          f.writeARFFPO(path+"SingleFunction"+GO+"ClNBPParNEnr.arff", GO, GOParents,goSem,mapFull ,cgmap,lor,categories,Double.POSITIVE_INFINITY,25, false,false,true,maxs,mins,fdrPvalsGZeroLOR,pvalsOR);///CL/Par/Enr OK!
          cl.classify(path+"SingleFunction"+GO+"ClNBPParNEnr.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"ClNBPParNEnr.out",rand,runInd,iter,4,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"ClNBPParNEnr.arff"));
          
          f.writeARFFPO(path+"SingleFunction"+GO+"ClBPPar.arff", GO, GOParents,goSem,mapFull ,cgmap,or,categories,Double.POSITIVE_INFINITY,12, false,false,true,maxs,mins);///CLPar OK!
          cl.classify(path+"SingleFunction"+GO+"ClBPPar.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"ClBPPar.out",rand,runInd,iter,11,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"ClBPPar.arff"));
         
         f.writeARFFPO(path+"SingleFunction"+GO+"MedBP.arff",GO,goSem,mapFull,cgmap,or,categories,Double.POSITIVE_INFINITY,1,false,false,true,maxs,mins);///Med
          cl.classify(path+"SingleFunction"+GO+"MedBP.arff", path+"SingleFunction"+GO+"Ind:"+(runInd+1)+"MedBP.out",rand,runInd,iter,5,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"MedBP.arff"));         
          
           f.writeARFFPO(path+"SingleFunction"+GO+"MedNBPPar.arff", GO, GOParents,goSem,mapFull ,cgmap,or,categories,Double.POSITIVE_INFINITY,21, false,false,true,maxs,mins);///Med/Par OK!
          cl.classify(path+"SingleFunction"+GO+"MedNBPPar.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"MedNBPPar.out",rand,runInd,iter,6,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"MedNBPPar.arff"));       
          
          f.writeARFFPO(path+"SingleFunction"+GO+"MedNBPParNEnr.arff", GO, GOParents,goSem,mapFull ,cgmap,lor,categories,Double.POSITIVE_INFINITY,24, false,false,true,maxs,mins,fdrPvalsGZeroLOR,pvalsOR);///Med/Par/Enr OK!
          cl.classify(path+"SingleFunction"+GO+"MedNBPParNEnr.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"MedNBPParNEnr.out",rand,runInd,iter,7,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"MedNBPParNEnr.arff"));
         
         f.writeARFFPO(path+"SingleFunction"+GO+"DistBP.arff",GO,goSem,mapFull,cgmap,or,categories,Double.POSITIVE_INFINITY,0,false,false,true,maxs,mins);//Dist
         cl.classify(path+"SingleFunction"+GO+"DistBP.arff", path+"SingleFunction"+GO+"Ind:"+(runInd+1)+"DistBP.out",rand,runInd,iter,8,seed,results);
         del.deleteFile(new File(path+"SingleFunction"+GO+"DistBP.arff"));
         
         f.writeARFFPO(path+"SingleFunction"+GO+"DistNBPPar.arff", GO, GOParents,goSem,mapFull ,cgmap,or,categories,Double.POSITIVE_INFINITY,22, false,false,true,maxs,mins);///Dist/Par OK!
          cl.classify(path+"SingleFunction"+GO+"DistNBPPar.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"DistNBPPar.out",rand,runInd,iter,9,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"DistNBPPar.arff"));
         
           f.writeARFFPO(path+"SingleFunction"+GO+"DistNBPParNEnr.arff", GO, GOParents,goSem,mapFull ,cgmap,lor,categories,Double.POSITIVE_INFINITY,23, false,false,true,maxs,mins,fdrPvalsGZeroLOR,pvalsOR);///Dist/Par/Enr OK!
          cl.classify(path+"SingleFunction"+GO+"DistNBPParNEnr.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"DistNBPParNEnr.out",rand,runInd,iter,10,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"DistNBPParNEnr.arff"));
         
          f.writeARFFPO(path+"SingleFunction"+GO+"ClBP.arff", GO, goSem,mapFull ,cgmap,or,categories,Double.POSITIVE_INFINITY,2, false,false,true,maxs,mins);
          cl.classify(path+"SingleFunction"+GO+"ClBP.arff", path+"SingleFunction"+GO+"Ind"+(runInd+1)+"ClBP.out",rand,runInd,iter,2,seed,results);
          del.deleteFile(new File(path+"SingleFunction"+GO+"ClBP.arff"));
       }     
     }
         
     System.out.println("Writing summary: ");
      try{
         FileWriter fw = new FileWriter(path+"Summary.out");
         
               for(int iter=0;iter<functions.size();iter++){
                            fw.write("Function: "+functions.get(iter)+"\n\n");
                            HashMap<Integer,ArrayList<Double>> values=results.get(iter);
                            
                            Iterator<Integer> it=values.keySet().iterator();
                            
                            while(it.hasNext()){
                                int key=it.next();
                                ArrayList<Double> tmp=values.get(key);
                                //6 - CL/CLPar
                                if(key==1)
                                fw.write("BPP: \n");
                                else if(key==2)
                                    fw.write("Cl: \n");
                                else if(key==3)
                                    fw.write("CL/CLPar: \n");
                                 else if(key==4)
                                    fw.write("CL/CLPar/Enr: \n");
                                else if(key==5)
                                    fw.write("Med: \n");
                                else if(key==6)
                                    fw.write("Med/Par: \n");
                                 else if(key==7)
                                    fw.write("Med/Par/Enr: \n");
                                else if(key==8)
                                    fw.write("Dist: \n");
                                else if(key==9)
                                    fw.write("Dist/Par: \n");
                                else if(key==10)
                                    fw.write("Dist/Par/Enr: \n");
                                else if(key==11)
                                    fw.write("ClPar: \n");
                                else if(key==12)
                                    fw.write("DO: \n");
                                else if(key==13)
                                    fw.write("BPOR>0Sig: \n");
                                for(int i=0;i<tmp.size();i++)
                                    fw.write(tmp.get(i)+" ");
                                fw.write("\n\n");
                            }  
                   }
                
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}    
      }
    }      

