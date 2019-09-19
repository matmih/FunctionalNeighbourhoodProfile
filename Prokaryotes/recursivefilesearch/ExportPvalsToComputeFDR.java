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
import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;
import org.apache.commons.math3.distribution.NormalDistribution;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class to compute p-values for each computed contingecy table
 */
public class ExportPvalsToComputeFDR {
    
    static public void main(String args[]){
        
        
         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");

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
        System.out.println("Loading complete");
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
     
        int random = 1;//Integer.parseInt(args[0]);
        
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
            
            FileWriter fw=null;
            FileWriter fw1=null;
            FileWriter fw2=null;
            FileWriter fw3=null;

   if(random == 0){         
    try{
        fw = new FileWriter(path+"PvalsAll"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }

    try{
        fw1= new FileWriter(path+"PvalsLORDiffAll"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
    try{
        fw2= new FileWriter(path+"PvalsLORZero"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }
    
    try{
        fw3 = new FileWriter(path+"PvalsLOROne"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }
   }
   else if(random == 1){
       try{
        fw = new FileWriter(path+"PvalsAllRand"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }

    try{
        fw1= new FileWriter(path+"PvalsLORDiffAllRand"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }

    try{
        fw2= new FileWriter(path+"PvalsLORZeroRand"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    }

    try{
        fw3 = new FileWriter(path+"PvalsLOROneRand"+".txt");   
    }
    catch(IOException e){
        e.printStackTrace();
    } 
   }
   
        HashMap<String,HashSet<String>> sigFunctions=new HashMap<>();
        NormalDistribution norm=new NormalDistribution();
        
        for(int iter=0;iter<mapFull.GOmap.keySet().size();iter++){  //main iteration starts..........
            String goName=mapFull.GOmapNumeric.get(iter);
            if(goName.contains("root")){
               System.out.println("Found root!");
                continue;
            }
     LoadOR lor=new LoadOR();
     File ORInp=null;
     if(random == 0)
        ORInp = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+mapFull.GOmapNumeric.get(iter)+".txt");
     else if(random == 1)
        ORInp = new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesRandomNew\\"+"RandomContingencyOut"+mapFull.GOmapNumeric.get(iter)+".txt");
     lor.loadContValue(ORInp);
     HashMap<String,Double> significant=new HashMap<>();
     HashMap<String,Double> significantLORDiff=new HashMap<>();
     HashMap<String,Double> significantLORgzero=new HashMap<>();
     HashMap<String,Double> significantLORgone=new HashMap<>();
    Iterator<String> it=lor.cont.keySet().iterator(); 
     
    ArrayList<Double> ctSelf=lor.cont.get(goName);
    Double orSelf=((ctSelf.get(0)+0.5)/(ctSelf.get(1)+0.5))/((ctSelf.get(2)+0.5)/(ctSelf.get(3)+0.5));
    Double logOrSelf=Math.log10(orSelf)/Math.log10(2);
     double SESelf=Math.sqrt(1.0/(ctSelf.get(0)+0.5)+1.0/(ctSelf.get(1)+0.5)+1.0/(ctSelf.get(2)+0.5)+1.0/(ctSelf.get(3)+0.5));
     
    for(int goindex=0;goindex<mapFull.GOmap.keySet().size();goindex++){
        String GOname=mapFull.GOmapNumeric.get(goindex);
        
        int iGO=mapFull.GOmap.get(GOname);
        if(GOname.contains("root")){
            System.out.println("skipped: "+GOname);
            continue;
        }
        int aInd=iter;

        ArrayList<Double> ct=lor.cont.get(GOname);
        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));//add continuity correction
        
        int min=Integer.MAX_VALUE;
        
        for(int i=0;i<4;i++)
            if(ct.get(i)<min)
                min=(int)(double)ct.get(i);
       
        double pval=1;
       
        int computeSig=1;
        
        if(min>10 && computeSig!=0){
        int val[][]=new int[2][2];
            val[0][0]=(int)(double)ct.get(0); val[0][1]=(int)(double)ct.get(1);
            val[1][0]=(int)(double)ct.get(2); val[1][1]=(int)(double)ct.get(3);
            ContingencyTable2x2 ctable=new ContingencyTable2x2(val);
            FishersExactTest fst=new FishersExactTest(ctable,jsc.tests.H1.NOT_EQUAL);
        
            pval=fst.getApproxSP();
        }
        else if(computeSig!=0){
                 try {
        myProb = myFisher.exact22((int)(double)ct.get(0), (int)(double)ct.get(1), (int)(double)ct.get(2), (int)(double)ct.get(3));
        myProb = myFisher.right; 
      } catch (IllegalArgumentException e) {
        myProb = Double.NaN;
      }
               pval=myProb;
        }
               significant.put(GOname, pval);
               
            Double logOr=Math.log10(or)/Math.log10(2);
           
         double SE=Math.sqrt(1.0/(ct.get(0)+0.5)+1.0/(ct.get(1)+0.5)+1.0/(ct.get(2)+0.5)+1.0/(ct.get(3)+0.5));
         
            double logdiff=logOr - logOrSelf;
            double jSE=Math.sqrt(SE*SE+SESelf*SESelf);
            double z=logdiff/jSE;
            double pVal=1.0-norm.cumulativeProbability(z);
               significantLORDiff.put(GOname,pVal);
               
              z=logOr/SE;
               pVal=1.0-norm.cumulativeProbability(z);
              significantLORgzero.put(GOname, pVal);    
              
              z=(logOr-1.0)/SE; 
              pVal=1.0-norm.cumulativeProbability(z);
              significantLORgone.put(GOname, pVal);  
    }

    boolean writeOnlyLOROne=false;
    
    if(writeOnlyLOROne==false){
      try{         
                            for(int itt=0;itt<mapFull.GOmap.keySet().size();itt++){
                                if(significant.containsKey(mapFull.GOmapNumeric.get(itt))){
                                    fw.write(significant.get(mapFull.GOmapNumeric.get(itt))+"");
                                fw.write("\n");
                                }
                            }
        }
        catch(Exception e)
                {e.printStackTrace();}    
      
      
       try{
                         
                            
                            for(int itt=0;itt<mapFull.GOmap.keySet().size();itt++){
                                if(significantLORDiff.containsKey(mapFull.GOmapNumeric.get(itt))){
                                    fw1.write(significantLORDiff.get(mapFull.GOmapNumeric.get(itt))+"");
                                fw1.write("\n");
                                }
                            }
        }
        catch(Exception e)
                {e.printStackTrace();}  
       
    
        try{
                         
                            
                            for(int itt=0;itt<mapFull.GOmap.keySet().size();itt++){
                                if(significantLORgzero.containsKey(mapFull.GOmapNumeric.get(itt))){
                                    fw2.write(significantLORgzero.get(mapFull.GOmapNumeric.get(itt))+"");
                                fw2.write("\n");
                                }
                            }
        }
        catch(Exception e)
                {e.printStackTrace();}  
        
        
    }
    
           try{       
                            for(int itt=0;itt<mapFull.GOmap.keySet().size();itt++){
                                if(significantLORgone.containsKey(mapFull.GOmapNumeric.get(itt))){
                                    fw3.write(significantLORgone.get(mapFull.GOmapNumeric.get(itt))+"");
                                fw3.write("\n");
                                }
                            }
        }
        catch(Exception e)
                {e.printStackTrace();}  
       
      
    }
             
            try{
                if(fw!=null)
                     fw.close();}
            catch(IOException e){
                    e.printStackTrace();
            }
            
             try{
                if(fw1!=null)
                     fw1.close();}
            catch(IOException e){
                    e.printStackTrace();
            }
             
              try{
                if(fw2!=null)
                     fw2.close();}
            catch(IOException e){
                    e.printStackTrace();
            }
             
              try{
                if(fw3!=null)
                     fw3.close();}
            catch(IOException e){
                    e.printStackTrace();
            }
                
    }
    
}
