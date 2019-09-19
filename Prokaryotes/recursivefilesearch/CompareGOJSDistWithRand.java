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
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class compare GOs of pairs of functions obtained on original and randomized dataset for segnificantly enriched, distant pairs of functions
 */
public class CompareGOJSDistWithRand {
    static public void main(String [] args ){
        
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
        File input2=new File("PP.arff");
        semanticsimilarity.GOMap mapFull1=new semanticsimilarity.GOMap();
        mapFull1.CreateGOMap(input2);
        mapFull1.printGOMap();
        mapFull1.loadFrequencies(fileWithUniprotFrequencyCounts);
     
          String path="";

        FisherExactTest myFisher = new FisherExactTest();   
       double myProb;
           HashSet<String> pkGOs=new HashSet<>();
            File prokaryotInput=new File("prokaryotic.txt");
            
            try (BufferedReader bufRdr = new BufferedReader(new FileReader(prokaryotInput)))
        {
            String line;
            
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                
                pkGOs.add(line);
                            
            }
        }
       catch(Exception e){
           e.printStackTrace();
       }
            
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
       
       File PvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsAll.txt");
       File PvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORDiffAll.txt");
       File FdrPvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsFDR.txt");
       File FdrPvalsLORDiff=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORDiffFDR.txt");
       File PvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsLORZero.txt");
       File FdrPvalsGZeroLOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pValsLORZeroFDR.txt");
       
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
             
             if(GOCount.get(go)<20)
                 continue;
             
              if(!pkGOs.contains(go) || !categories.get(0).contains(go)){
                    continue;
            }

             double ICTest=mapMod.frequency.get(mapMod.GOmap.get(go));
             
        if(Math.abs(Math.log10(ICTest)/Math.log10(2))<4)
            continue;
             
             if(average.get(go).getValue0()>locNeigh.get(go).getValue0() && average.get(go).getValue0()>=1.5*baseline.get(go).getValue0())
                     functions.add(go);             
         }
         
         System.out.println("Broj funkcija: "+functions.size());
          
        DataGOFunctionInfo datInfo=new DataGOFunctionInfo();
        DataGOFunctionInfo datInfoRand=new DataGOFunctionInfo();
        datInfo.computeOccurence(cgmap);
        datInfoRand.randomizedOccurence(cgmap, datInfo);
        DataGOFunctionInfo datInfo1=new DataGOFunctionInfo();
        DataGOFunctionInfo datInfo1Rand=new DataGOFunctionInfo();
        datInfo1.computeOccurence(cgmap30);//create a function in computeOccurence that randomizes COG-function relation
        datInfo1Rand.randomizedOccurence(cgmap30, datInfo1);

             double maxJS=0.0; String maxGO1="",maxGO2="";
        HashMap<String,ArrayList<String>> sigPairs=new HashMap<>();
          
        
        ArrayList<Double> regJaccard = new ArrayList<>();
        ArrayList<Double> randJaccard = new ArrayList<>();
        ArrayList<Double> regJaccard30 = new ArrayList<>();
        ArrayList<Double> randJaccard30 = new ArrayList<>();
        
        int criteria=1;//0-strict (logOR(GOx,GOy)>logOR(GOx,GOx), sigg.), 1-medium logOR(GOx,GOy)>0 sigg.
        int countPairs=0;
        for(int iter=0;iter<functions.size();iter++){ 
            ArrayList<String> tmp=new ArrayList<>();
            sigPairs.put(functions.get(iter), tmp);
            if(!pkGOs.contains(functions.get(iter))){
                System.out.println("Not contained in prokaryotic: "+iter);
                    continue;
            }
            
        int aInd=mapFull.GOmap.get(functions.get(iter));
        int aInd1=cgmap.GOtoIndex.get(functions.get(iter));
        int aInd2=cgmap30.GOtoIndex.get(functions.get(iter));
        double pvalOROrig=pvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(aInd);
        double pValFDROROrig=fdrPvalsOR.get(mapFull.GOmap.get(functions.get(iter))).get(aInd);
        
        if(pvalOROrig>0.05 || pValFDROROrig>0.2 )
            continue;
        
        double IC=0.0;
        IC=mapMod.frequency.get(mapMod.GOmap.get(functions.get(iter)));
        if(Math.abs(Math.log10(IC)/Math.log10(2))<4){
            System.out.println("To general: "+iter);
            continue;
        }
        
         LoadOR lor=new LoadOR();
     File ORInp=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+mapFull.GOmapNumeric.get(aInd)+".txt");
     lor.loadContValue(ORInp);
     
      ArrayList<Double> ctSelf=lor.cont.get(functions.get(iter));
    Double orSelf=((ctSelf.get(0)+0.5)/(ctSelf.get(1)+0.5))/((ctSelf.get(2)+0.5)/(ctSelf.get(3)+0.5));
    Double logOrSelf=Math.log10(orSelf)/Math.log10(2);
    
    if(logOrSelf<0.0){
        System.out.println("Lor <0: "+iter);
        continue;
    }
    
    double SESelf=Math.sqrt(1/(ctSelf.get(0)+0.5)+1/(ctSelf.get(1)+0.5)+1/(ctSelf.get(2)+0.5)+1/(ctSelf.get(3)+0.5));
    
    double confUSelf=logOrSelf+1.96*SESelf;
    double confLSelf=logOrSelf-1.96*SESelf;
    
    for(int goindex=0;goindex<mapFull.GOmap.keySet().size();goindex++){
        String GOname=mapFull.GOmapNumeric.get(goindex);
        int iGO=mapFull.GOmap.get(GOname);
        int iGO1=cgmap.GOtoIndex.get(GOname);
        int iGO2=cgmap30.GOtoIndex.get(GOname);
        if(!pkGOs.contains(GOname) || !categories.get(0).contains(GOname))
            continue;
        
        
        ArrayList<Double> ct=lor.cont.get(GOname);
        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
        Double logOr=Math.log10(or)/Math.log10(2);
        
         double pvalORComp=pvalsLORDiff.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         double pvalLORDiffFDR=fdrPvalsLORDiff.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
          
         double pvalORSelf=pvalsOR.get(mapFull.GOmap.get(GOname)).get(goindex);
         double pvalORFDRORSelf=fdrPvalsOR.get(goindex).get(goindex);
         
          double pvalLORGZero=fdrPvalsGZeroLOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
          double pvalGZero=pvalsGZeroLOR.get(mapFull.GOmap.get(functions.get(iter))).get(goindex);
         
          if(pvalORSelf>0.05 || pvalORFDRORSelf>0.2)
              continue;
          
          double ICSelf=0.0;
        ICSelf=mapMod.frequency.get(mapMod.GOmap.get(GOname));
        if(Math.abs(Math.log10(ICSelf)/Math.log10(2))<4)
            continue;
        
        
        double SE=Math.sqrt(1/(ct.get(0)+0.5)+1/(ct.get(1)+0.5)+1/(ct.get(2)+0.5)+1/(ct.get(3)+0.5));
    
        double confU=logOr+1.96*SE;
        double confL=logOr-1.96*SE;
    
        if(logOr>0.0 && logOr>logOrSelf && pvalORComp<0.05 && pvalLORDiffFDR<0.2 && goSem.similarity[iGO][aInd]<1 && criteria==0){//remove to <1
           sigPairs.get(functions.get(iter)).add(GOname);
            System.out.println(functions.get(iter)+"\t"+" ("+goDescription.get(functions.get(iter))+") "+"\t"+logOrSelf+" +-"+(1.96*SESelf)+" "+"\t"+GOname+"\t"+" ("+goDescription.get(GOname)+") "+"\t"+logOr+" +-"+(1.96*SE)+"\t"+goSem.similarity[iGO][aInd]+"\t"+datInfo.computeJaccard(aInd1, iGO1)+"\t"+datInfo.computeF1(aInd1, iGO1)+"\t"+datInfo1.computeJaccard(aInd2, iGO2)+"\t"+datInfo1.computeF1(aInd2, iGO2));
            regJaccard.add(datInfo.computeJaccard(aInd1, iGO1));
            regJaccard30.add(datInfo1.computeJaccard(aInd2, iGO2));
            randJaccard.add(datInfoRand.computeJaccard(aInd1, iGO1));
            randJaccard30.add(datInfo1Rand.computeJaccard(aInd2, iGO2));
        }
        else if(logOr>0.0 && pvalGZero<0.05 && pvalLORGZero<0.1 && logOrSelf>0.0 /*&& pvalORComp<0.05 */&& goSem.similarity[iGO][aInd]<1 && criteria==1){
            sigPairs.get(functions.get(iter)).add(GOname);
            double jstmp=datInfo.computeJaccard(aInd1, iGO1);
            if(jstmp>maxJS){
                maxJS=jstmp;
                maxGO1=functions.get(iter);
                maxGO2=GOname;
            }
             System.out.println(functions.get(iter)+"\t ("+goDescription.get(functions.get(iter))+") \t"+df.format(logOrSelf)+" +- "+df.format((1.96*SESelf))+"\t  "+GOname+"\t ("+goDescription.get(GOname)+")\t "+df.format(logOr)+" +- "+df.format((1.96*SE))+"\t "+df.format(goSem.similarity[iGO][aInd])+"\t "+df.format(datInfo.computeJaccard(aInd1, iGO1))+"\t "+df.format(datInfo.computeF1(aInd1, iGO1))+"\t "+df.format(datInfo1.computeJaccard(aInd2, iGO2))+"\t "+df.format(datInfo1.computeF1(aInd2, iGO2)));
         regJaccard.add(datInfo.computeJaccard(aInd1, iGO1));
            regJaccard30.add(datInfo1.computeJaccard(aInd2, iGO2));
            randJaccard.add(datInfoRand.computeJaccard(aInd1, iGO1));
            randJaccard30.add(datInfo1Rand.computeJaccard(aInd2, iGO2));
        }
        
        if(goSem.similarity[iGO][aInd]<1 )
            countPairs++;
    }
     }
        
        //write to file norm,rand,norm30,rand30
       try{
           String output = "JSDistWithRand.txt";
        FileWriter fw = new FileWriter(output); 
        for(int i=0;i<regJaccard.size();i++){
                          
                fw.write(regJaccard.get(i)+","+randJaccard.get(i)+","+regJaccard30.get(i)+","+randJaccard30.get(i)+"\n");
        }
       fw.close();
       }
            catch(Exception e)
                {e.printStackTrace();}  
       }        
        
    }
