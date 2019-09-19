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
import java.util.Iterator;
import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create heatmap data
 */
public class CreateHeatmapData {
    static public void main(String arg[]){
        
        File prokaryotInput=new File("prokaryotic.txt");
        HashSet<String> pkGOs=new HashSet<>();
        COGGOMap cgmap=new COGGOMap();
           
        
         File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
         File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
         String redGO="GORed5.txt"; String redCog="COGRed150.txt";//CogReduced50,GOReducedDataset5Occ
         
          GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);

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
        
         File semFile=new File("PairwiseGOSim.txt");
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
         
       HashMap<String, Integer> GOoccurence=new HashMap<>();  
       for(String s:cgmap.GOtoIndex.keySet()){
           int count=0;
           for(String cog:cgmap.CogGOmap.keySet()){
               ArrayList<String> GOs=cgmap.CogGOmap.get(cog);
               if(GOs.contains(s))
                   count++;
           }
           GOoccurence.put(s, count);
       } 
       
       HashSet<String> satisfGOs=new HashSet<String>();
       
       for(String s:cgmap.GOtoIndex.keySet())
           if(pkGOs.contains(s) && categories.get(0).contains(s))//>=5
               satisfGOs.add(s);
   
       double logORs[][]=new double[satisfGOs.size()][satisfGOs.size()];
       double semSim[][]=new double[satisfGOs.size()][satisfGOs.size()];
       
       HashMap<String,Integer> tmpIndex=new HashMap<>();
       HashMap<Integer,String> indexTmp=new HashMap<>();
       int ind=0;
       
       for(String s:satisfGOs){
           tmpIndex.put(s, ind);
           indexTmp.put(ind, s);
           ind++;
       }
       
    FisherExactTest myFisher = new FisherExactTest();   
       double myProb;   
       
        HashMap<Integer,ArrayList<Double>> pvalsOR=new HashMap<>();
          
        int computeSig=1;
        
        for(int iter=0;iter<cgmap.GOtoIndex.keySet().size();iter++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
            pvalsOR.put(iter, tmp);
        }
       
        File PvalsOR=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\PvalsAll.txt");
      
       
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
                
    int processed=0;   
  for(String go:satisfGOs){  
     LoadOR lor=new LoadOR();
     File ORInp=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+go+".txt");
     lor.loadContValue(ORInp);
    
    Iterator<String> it=lor.cont.keySet().iterator(); 
    ArrayList<String> goOrdered=new ArrayList<>();
    for(int goindex=0;goindex<lor.cont.keySet().size();goindex++){
        String GOname=mapFull.GOmapNumeric.get(goindex); 
        if(!satisfGOs.contains(GOname))
            continue;
        goOrdered.add(GOname);
        int iGO=mapFull.GOmap.get(GOname);
        if(GOname.contains("root"))//fix this
            continue;
        int aInd=mapFull.GOmap.get(go);
        ArrayList<Double> ct=lor.cont.get(GOname);
        
       
        Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
        Double logOr=Math.log10(or)/Math.log10(2);
        
        double pval=pvalsOR.get(mapFull.GOmap.get(go)).get(goindex);
        
      if(pval>0.05 && computeSig!=0)
          logOr=Double.NaN;

      LoadOR lor1=new LoadOR();
     File ORInp1=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+GOname+".txt");
     lor1.loadContValue(ORInp1);
     
     
      ArrayList<Double> ct1=lor1.cont.get(go);

      
      pval=pvalsOR.get(goindex).get(mapFull.GOmap.get(go));
        
        Double or1=((ct1.get(0)+0.5)/(ct1.get(1)+0.5))/((ct1.get(2)+0.5)/(ct1.get(3)+0.5));
        Double logOr1=Math.log10(or1)/Math.log10(2);
        
      if(pval>0.05 && computeSig!=0)
          logOr1=Double.NaN;
     
        
        logORs[tmpIndex.get(go)][tmpIndex.get(GOname)]=logOr;
        logORs[tmpIndex.get(GOname)][tmpIndex.get(go)]=logOr1;//nije simetricno!, ispraviti!
        semSim[tmpIndex.get(go)][tmpIndex.get(GOname)]=goSem.similarity[iGO][aInd];
        semSim[tmpIndex.get(GOname)][tmpIndex.get(go)]=goSem.similarity[aInd][iGO];
       }
    processed++;
        System.out.println("Processed: "+processed);
     }
       
  String output="HeatLRMat.txt";
  int outputHeader=0;
  
  try
            {
                FileWriter fw = new FileWriter(output);
                
               if(outputHeader==1){ 
                for(int i=0;i<satisfGOs.size();i++)
                    if(i+1<satisfGOs.size())
                        fw.write(indexTmp.get(i)+" ");
                    else fw.write(indexTmp.get(i)+"\n");
               }
                
                for(int i=0;i<satisfGOs.size();i++){
                    if(outputHeader==1)
                    fw.write(indexTmp.get(i)+" ");
                    for(int j=0;j<satisfGOs.size();j++)
                        if((j+1)<satisfGOs.size())
                        fw.write(logORs[i][j]+" ");
                        else
                            fw.write(logORs[i][j]+"\n");
                }
                fw.close();
                
                output="HeatSemSimMat.txt";
                
                fw=new FileWriter(output);
                
                 for(int i=0;i<satisfGOs.size();i++)
                    for(int j=0;j<satisfGOs.size();j++)
                        if((j+1)<satisfGOs.size())
                        fw.write(semSim[i][j]+" ");
                        else
                            fw.write(semSim[i][j]+"\n");
                fw.close();
            }
           catch(Exception e){
                e.printStackTrace();
            }
    }
}
