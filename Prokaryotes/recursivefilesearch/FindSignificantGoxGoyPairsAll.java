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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.javatuples.Pair;
import static recursivefilesearch.COGGOMap.ENCODING;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to store all pairs of GOs with significant LOR (with respect to one of the two criteria of significance)
 */
public class FindSignificantGoxGoyPairsAll {
     static public void main(String arg[]){

         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File OgGOs30=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_30_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         File semFile=new File("PairwiseGOSim.txt");
         
         String redCog="COGRed150.txt"; String redGO="GORed5.txt";

        GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

        COGGOMap cgmap=new COGGOMap();
        cgmap.createCOGGOMapping(OgGOs);

        
        COGGOMap cgmap30=new COGGOMap();
        cgmap30.createCOGGOMapping(OgGOs30);
        
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
        cgmap30.reduceMap(red);
        geneOGMap.reduceMap(red);

        
        HashMap<String,String> goDescription=new HashMap<>();
       Iterator<String> itDesc=cgmap.GOtoIndex.keySet().iterator();
       
       while(itDesc.hasNext()){
           String go=itDesc.next();
           int goNum=Integer.parseInt(go.replaceAll("GO", ""));
           String desc=myGO.get(goNum).getName();
           
           goDescription.put(go, desc);
       }
       
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
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
     
          String path="";

        FisherExactTest myFisher = new FisherExactTest();   
            
        ArrayList<String> functions=new ArrayList<>();      
       HashMap<Integer,ArrayList<Double>> pvalsOR=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> pvalsLORDiff=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> fdrPvalsOR=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> fdrPvalsLORDiff=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> pvalsGZeroLOR=new HashMap<>();
       HashMap<Integer,ArrayList<Double>> fdrPvalsGZeroLOR=new HashMap<>();   
       
        for(int iter=0;iter<cgmap.GOtoIndex.keySet().size();iter++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvalsOR.put(iter, tmp);
            tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvalsLORDiff.put(iter, tmp);
            tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            fdrPvalsOR.put(iter, tmp);
            tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            fdrPvalsLORDiff.put(iter, tmp);
            tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvalsGZeroLOR.put(iter, tmp);
            tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            fdrPvalsGZeroLOR.put(iter, tmp);
        }
        
        int rand = 0;//0-non random, 1-random
       
       File PvalsOR = null;  File PvalsLORDiff = null;
       File FdrPvalsOR = null; File FdrPvalsLORDiff = null;
       File PvalsGZeroLOR = null; File FdrPvalsGZeroLOR = null;
       
       if(rand == 0){//pvals on original data
              PvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsAll.txt");
              PvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORDiffAll.txt");
              FdrPvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsFDR.txt");
              FdrPvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORDiffFDR.txt");
              PvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORZero.txt");
              FdrPvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORZeroFDR.txt");
       }
       else{//pvals after randomization
           PvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsAllRand.txt");
              PvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORDiffAllRand.txt");
              FdrPvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsFDRRand.txt");
              FdrPvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORDiffFDRRand.txt");
              PvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORZeroRand.txt");
              FdrPvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORZeroFDRRand.txt");
       }
       
         BufferedReader bufRdr=null;
           try 
        {
            bufRdr = new BufferedReader(new FileReader(PvalsOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double pvalOR=Double.parseDouble(line.trim());
                pvalsOR.get(count).set(count1, pvalOR);
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
           
           try 
        {
            bufRdr = new BufferedReader(new FileReader(PvalsLORDiff));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double pvalsLORdiff=Double.parseDouble(line.trim());
                pvalsLORDiff.get(count).set(count1, pvalsLORdiff);
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
       
           
           
           try 
        {
            bufRdr = new BufferedReader(new FileReader(FdrPvalsOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrpvalsor=Double.parseDouble(line.trim());
                fdrPvalsOR.get(count).set(count1, fdrpvalsor);
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
       
           
           
           try 
        {
            bufRdr = new BufferedReader(new FileReader(FdrPvalsLORDiff));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double fdrpvalslordiff=Double.parseDouble(line.trim());
                fdrPvalsLORDiff.get(count).set(count1, fdrpvalslordiff);
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
           
           
           try 
        {
            bufRdr = new BufferedReader(new FileReader(PvalsGZeroLOR));
            String line;
            String cog="";
            int count=0,count1=0,lineNum=0;
            while ((line = bufRdr.readLine()) != null)
            {
                lineNum++;
                Double pvalsGZerolor=Double.parseDouble(line.trim());
                pvalsGZeroLOR.get(count).set(count1, pvalsGZerolor);
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

           
        HashMap<String,Pair<Double,Double>> average=new HashMap<>();
        HashMap<String,Pair<Double,Double>> baseline=new HashMap<>();
        HashMap<String,Pair<Double,Double>> locNeigh=new HashMap<>();
        File functionFile= new File("BaselineOGs_600_.oob");
        File locNeighFile= new File("DistanceLocD_600_.oob");
        File baselineFile = new File("genome5k=5logDist.out");
        HashMap<String,Double> realFreq=new HashMap<>();
        
        try (BufferedReader bufRdr1 = new BufferedReader(new FileReader("Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt")))
        {
            String line;
            
            while ((line = bufRdr1.readLine()) != null)
            {
                if (line.startsWith("#")) //header
                    continue;
                
                String[] parts = line.split("\t");

                int go = Integer.parseInt(parts[0]);
                String gString=go+"";
                while(gString.length()<7){
                    gString="0"+gString;
                }
                gString="GO"+gString;
                double generality = Double.parseDouble(parts[1]);
                realFreq.put(gString, generality);
                System.out.println("GO: "+gString);
            }
        }
        catch(Exception e){}

        BufferedReader reader;
         try {
             
      Path path1 =Paths.get(functionFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              average.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
         
         
      try {
             
      Path path1 =Paths.get(locNeighFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int oob=0, dataSection=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Out-Of-Bag")){
              oob=1;
              continue;
          }
          
           if(line.contains("Default") || line.equals(""))
              dataSection=0;
          
          if(oob==1){
              if(line.contains("P100R")){
                  dataSection=1;
                  oob=0;
                  continue;
              }
          }

          if(dataSection==1){
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              locNeigh.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)/*Double.parseDouble(freq)*/));
          }
          
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }     

         try {
      Path path1 =Paths.get(baselineFile.getAbsolutePath());

      reader = Files.newBufferedReader(path1,ENCODING);
      String line = null;
      int summary=0, dataSection=0, test=0;
      while ((line = reader.readLine()) != null) {
          
          if(line.contains("Training error")){
              summary=1;
              continue;
          }
          
           if(line.contains("P100R") && summary==1){
               System.out.println("Data section reached baseline");
                  dataSection=1;
                 summary=0;
                  continue;
              }
          
           if(line.contains("Testing error") || line.equals(""))
              dataSection=0;

          if(dataSection==1){
              System.out.println("Entered data section");
              line=line.trim();
              String tmp[]=line.split(", ");
              System.out.println(tmp[0]);
              String first=(tmp[0].split(": "))[1];
              first=first.split("\\[")[0];
              System.out.println("first baseline: "+first);
              String AUPRC=tmp[2];
              AUPRC=tmp[2].split(": ")[1];
              String freq=tmp[3].split(": ")[1];
              System.out.println("AUPRC: "+AUPRC);
              System.out.println("freq: "+freq);
              baseline.put(first, new Pair(Double.parseDouble(AUPRC),realFreq.get(first)));
          }
      }
      reader.close();
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            } 
        
         System.out.println("Average: "+average.keySet().size());
         System.out.println("DistLoc: "+locNeigh.keySet().size());
         System.out.println("Baseline: "+baseline.keySet().size());
         DecimalFormat df = new DecimalFormat("#.###");
         HashMap<String,Integer> GOCount=new HashMap<>();
        
         HashSet<String> functionsCond = new HashSet<>();
         
        for(String s:cgmap.CogGOmap.keySet()){
            System.out.println(s+" \n");
            if(!cgmap.CogGOmap.containsKey(s))
                continue;
            for(String go:cgmap.CogGOmap.get(s))
                if(!GOCount.containsKey(go))
                    GOCount.put(go, 1);
                else{
                    int count=GOCount.get(go);
                    count=count+1;
                    GOCount.put(go, count);
                }       
       }
               
         Iterator<String> itFunc=average.keySet().iterator();
         
         while(itFunc.hasNext()){
             String go=itFunc.next();

             double ICTest=mapMod.frequency.get(mapMod.GOmap.get(go));
             
                     functions.add(go);
             
         }
         
         System.out.println("Broj funkcija: "+functions.size());
         int test=1;
     
        DataGOFunctionInfo datInfo=new DataGOFunctionInfo();
        datInfo.computeOccurence(cgmap);
        DataGOFunctionInfo datInfo1=new DataGOFunctionInfo();
        datInfo1.computeOccurence(cgmap30);
        
        double maxJS=0.0; String maxGO1="",maxGO2="";
        HashMap<String,ArrayList<String>> sigPairs=new HashMap<>();
          
        int criteria=1;//0-strict (logOR(GOx,GOy)>logOR(GOx,GOx), sigg.), 1-medium logOR(GOx,GOy)>0 sigg.
        int countPairs=0;
        
        
        try{
        
            FileWriter fw = null;
            
            if(rand == 0)
                fw = new FileWriter("SignificantAll01Prokaryot.txt");
            else fw = new FileWriter("SignificantPairsBPDistRand01.txt");
            
            String header = "GO1\tDescription1\tLogOR GO1-GO1 CI1\tGO2\tDescription2\tLOGOR GO1-GO2 CI2\tPvalOR>0\tResnik SemDist\tJaccI_50\tF1_50\tJaccI_30\tF1_30\n";
            
            if(rand==1)
                 fw.write(header);
            else 
                fw.write("GO1\tDescription1\tLogORGO1-GO1\tsign\tCI1\tGO2\tDescription2\tLOGORGO1-GO2\tsign\tCI2\tLOGORGO1-GO2Rand\tPvalOR>0\tResnik SemDist\tJaccI_50\tF1_50\tJaccI_30\tF1_30\n");
            
            HashSet<String> bpOntol = new HashSet<>();
            
        for(int iter=0;iter<functions.size();iter++){ 
            ArrayList<String> tmp=new ArrayList<>();
            sigPairs.put(functions.get(iter), tmp);
            
            if(categories.get(0).contains(functions.get(iter)))
                bpOntol.add(functions.get(iter));
            
        int aInd=mapFull.GOmap.get(functions.get(iter));
        int aInd1=cgmap.GOtoIndex.get(functions.get(iter));
        int aInd2=cgmap30.GOtoIndex.get(functions.get(iter));
        double pvalOROrig=pvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(aInd);
        double pValFDROROrig=fdrPvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(aInd);
        
        double IC=0.0;
        IC=mapMod.frequency.get(mapMod.GOmap.get(functions.get(iter)));
        
         LoadOR lor=new LoadOR();
         LoadOR lorRand=new LoadOR();
     File ORInp=null;
     File ORInpRand=null;
     
     if(rand == 0){
        ORInp = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
        ORInpRand = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesRandomNew\\"+"RandomContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
     }
     else
          ORInp = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesRandomNew\\"+"RandomContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
     lor.loadContValue(ORInp);
     if(rand ==0)
         lorRand.loadContValue(ORInpRand);
     
      ArrayList<Double> ctSelf=lor.cont.get(functions.get(iter));
    Double orSelf=((ctSelf.get(0)+0.5)/(ctSelf.get(1)+0.5))/((ctSelf.get(2)+0.5)/(ctSelf.get(3)+0.5));
    Double logOrSelf=Math.log10(orSelf)/Math.log10(2);
    
    double SESelf=Math.sqrt(1/(ctSelf.get(0)+0.5)+1/(ctSelf.get(1)+0.5)+1/(ctSelf.get(2)+0.5)+1/(ctSelf.get(3)+0.5));
    
    double confUSelf=logOrSelf+1.96*SESelf;
    double confLSelf=logOrSelf-1.96*SESelf;
    
    for(int goindex=0;goindex<mapFull.GOmap.keySet().size();goindex++){
        String GOname=mapFull.GOmapNumeric.get(goindex);
        int iGO=mapFull.GOmap.get(GOname);
        int iGO1=cgmap.GOtoIndex.get(GOname);
        int iGO2=cgmap30.GOtoIndex.get(GOname);    
        
        ArrayList<Double> ct=lor.cont.get(GOname);
         ArrayList<Double> ctRand=null;
         
         if(rand ==0)
            ctRand = lorRand.cont.get(GOname);
         
        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
        Double logOr=Math.log10(or)/Math.log10(2);
        
        Double orRand, logOrRand=-1.0 ;
        
        if(rand ==0){
            orRand=((ctRand.get(0)+0.5)/(ctRand.get(1)+0.5))/((ctRand.get(2)+0.5)/(ctRand.get(3)+0.5));
            logOrRand=Math.log10(orRand)/Math.log10(2);
        }
                
         double pvalORComp=pvalsLORDiff.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         double pvalLORDiffFDR=fdrPvalsLORDiff.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
          
         double pvalORSelf=pvalsOR.get(mapFull.GOmap.get(GOname)).get(goindex);
         double pvalORFDRORSelf=fdrPvalsOR.get(goindex).get(goindex);
         
         double pvalOR = pvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         double pvalORFDR = fdrPvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         
          double pvalLORGZero=fdrPvalsGZeroLOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
          double pvalGZero=pvalsGZeroLOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         
          if(pvalOR>0.05 || pvalORFDR>0.2)
              continue;
          
          double ICSelf=0.0;
        ICSelf=mapMod.frequency.get(mapMod.GOmap.get(GOname));

        double SE=Math.sqrt(1/(ct.get(0)+0.5)+1/(ct.get(1)+0.5)+1/(ct.get(2)+0.5)+1/(ct.get(3)+0.5));
    
        double confU=logOr+1.96*SE;
        double confL=logOr-1.96*SE;
        
        if(logOr>0.0 && logOr>logOrSelf && pvalORComp<0.05 && pvalLORDiffFDR<0.2 && criteria==0){//remove to <1 //very strict criteria, LOR(GOx,GOy) must be > LOR(Gox,GOx) and significant
           sigPairs.get(functions.get(iter)).add(GOname);
            System.out.println(functions.get(iter)+"\t"+" ("+goDescription.get(functions.get(iter))+") "+"\t"+logOrSelf+" +-"+(1.96*SESelf)+" "+"\t"+GOname+"\t"+" ("+goDescription.get(GOname)+") "+"\t"+logOr+" +-"+(1.96*SE)+"\t"+goSem.similarity[iGO][aInd]+"\t"+datInfo.computeJaccard(aInd1, iGO1)+"\t"+datInfo.computeF1(aInd1, iGO1)+"\t"+datInfo1.computeJaccard(aInd2, iGO2)+"\t"+datInfo1.computeF1(aInd2, iGO2));
        }
        else if((logOr>0.0 && pvalGZero<0.05 && pvalLORGZero<0.1 &&  criteria==1)){//less strict criteria, LOR(GOx,GOy) significantly >0
               sigPairs.get(functions.get(iter)).add(GOname);
          
            double jstmp=datInfo.computeJaccard(aInd1, iGO1);
            if(jstmp>maxJS){
                maxJS=jstmp;
                maxGO1=functions.get(iter);
                maxGO2=GOname;
            }
            fw.write(functions.get(iter)+"\t("+goDescription.get(functions.get(iter))+")\t"+df.format(logOrSelf)+"\t+-\t"+df.format((1.96*SESelf))+"\t"+GOname+"\t("+goDescription.get(GOname)+")\t"+df.format(logOr)+"\t+-\t"+df.format((1.96*SE))+"\t"+df.format(logOrRand)+"\t"+pvalGZero+"\t"+df.format(goSem.similarity[iGO][aInd])+"\t"+df.format(datInfo.computeJaccard(aInd1, iGO1))+"\t"+df.format(datInfo.computeF1(aInd1, iGO1))+"\t"+df.format(datInfo1.computeJaccard(aInd2, iGO2))+"\t"+df.format(datInfo1.computeF1(aInd2, iGO2))+"\n");
        }

            countPairs++;
    }

     }

        fw.close();
        System.out.println();
        System.out.println("Biological process ontology num nodes: "+bpOntol.size());
        System.out.println("Num functions satisfying condition: "+functionsCond.size());
        }
        catch(IOException e){
            e.printStackTrace();         
        }
    }
}
