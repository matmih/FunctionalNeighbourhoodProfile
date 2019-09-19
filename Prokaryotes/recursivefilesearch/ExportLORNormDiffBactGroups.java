/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to export LORs divided by groups of bacteria: free-living, pathogenic in mammals
 */
public class ExportLORNormDiffBactGroups {
    public static void main(String []arg){
           COGGOMap cgmap=new COGGOMap();
        
       File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
       File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
       String redCog="/home/mmihelcic/GOdatasetGenerator/GOanotacije/COGRed150.txt"; String redGO="/home/mmihelcic/GOdatasetGenerator/GOanotacije/GORed5.txt";
         
          GeneOGMapping geneOGMap=new GeneOGMapping();

        geneOGMap.loadGOGMapping(geneOgs);
        cgmap.createCOGGOMapping(OgGOs);
        
        ReducedGOTranslator rgt=new ReducedGOTranslator();
        rgt.ReadAndTranslate(new File(redGO));

    OntologyTools.GeneOntology myGO=null;

    try{
    myGO= new GeneOntology("/home/mmihelcic/GOdatasetGenerator/GOanotacije/go_201401-termdb.obo-xml");
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
                  
         File semFile=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/PairwiseGOSim.txt");
         GOSemanticSimilarity goSem=new GOSemanticSimilarity(mapFull.GOmap.keySet().size());
         goSem.readSimilarity(semFile);
        
          ArrayList<HashSet<String>> categories=new ArrayList<>();
        GODevider gd=new GODevider();
        try{
            categories=gd.ExtractCategories(new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/go_201401-termdb.obo-xml.gz"),cgmap);
        }
        catch(Exception e){
            e.printStackTrace();
        }  
         
        semanticsimilarity.GOMap mapMod=new semanticsimilarity.GOMap();
        mapMod.CreateGOMap(cgmap);
        File fileWithUniprotFrequencyCounts=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/Uniprot-freqs-2070_organisms-Uniprot-GOA-2013-12-10.txt");
        mapMod.loadFrequencies(fileWithUniprotFrequencyCounts);
        File input2=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/PP.arff");
        semanticsimilarity.GOMap mapFull1=new semanticsimilarity.GOMap();
        mapFull1.CreateGOMap(input2);
        mapFull1.printGOMap();
        mapFull1.loadFrequencies(fileWithUniprotFrequencyCounts);

          String path="";

        ArrayList<String> functions=new ArrayList<>();      
       HashMap<Integer,ArrayList<Double>> pvals=new HashMap<>();
        
      int countZeros = 0;

      
     String name = "LOR";
     int randStats = 0;
      FileWriter fw = null;
      FileWriter fw1 = null;
      
      String extension = "";
      
      for(int z=0;z<2;z++){
             countZeros = 0;
          if(z==0)
              extension = "FL";
          else extension = "PM";
          
           try{  
               
                    fw = new FileWriter(name.split("\\.")[0]+"Norm"+extension+".txt");
                    fw1 = new FileWriter(name.split("\\.")[0]+"RandNorm"+extension+".txt");
           }
           catch(IOException e){
               e.printStackTrace();
           }
            LoadOR lor=new LoadOR();
             LoadOR lorRand=new LoadOR();
             
            String contingencyPath = "";
            String contingencyPathRand = "";
        
       if(randStats == 0){
             contingencyPath = "/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContAllSel/";
             contingencyPathRand = "/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContRandAllSel/";
       }
       else if(randStats == 1){
             contingencyPath = "/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContRandAllSel/";
         }
         
         ArrayList<Double> lors = new ArrayList<>();
         ArrayList<Double> lorsRand = new ArrayList<>();

               for(int iter=0;iter<cgmap.GOtoIndex.keySet().size();iter++){ 
                   File ORInp=null;
                   File ORInpRand=null;
                   if(randStats == 0){
                        ORInp=new File(contingencyPath+"ContingencyOut"+extension+cgmap.IndexToGO.get(iter)+".txt");
                        ORInpRand=new File(contingencyPathRand+"RandomContingencyOut"+extension+cgmap.IndexToGO.get(iter)+".txt");
                   }
                   else if(randStats == 1)
                       ORInp=new File(contingencyPath+"RandomContingencyOut"+extension+cgmap.IndexToGO.get(iter)+".txt");
                    lor.loadContValue(ORInp);
                    lorRand.loadContValue(ORInpRand);
                   for(int goindex=0;goindex<cgmap.GOtoIndex.keySet().size();goindex++){
                       ArrayList<Double> ct=lor.cont.get(cgmap.IndexToGO.get(goindex));
                       ArrayList<Double> ct1=lorRand.cont.get(cgmap.IndexToGO.get(goindex));
                         Double or=((ct.get(0)+0.5)/(ct.get(1)+0.5))/((ct.get(2)+0.5)/(ct.get(3)+0.5));
                         Double logOr=Math.log10(or)/Math.log10(2);
                         Double orRand =((ct1.get(0)+0.5)/(ct1.get(1)+0.5))/((ct1.get(2)+0.5)/(ct1.get(3)+0.5));
                         Double logOrRand=Math.log10(orRand)/Math.log10(2);
                            if(ct.get(0)!=0){
                                lors.add(logOr);
                                lorsRand.add(logOrRand);
                            }
                            else countZeros++;
                   }
               }
               
               try{
           
                for(int k=0;k<lors.size();k++){
                    fw.write(lors.get(k)+"\n");
                    fw1.write(lorsRand.get(k)+"\n");
                }
                
                fw.close();
                fw1.close();
          }
          catch(IOException e){
             e.printStackTrace();
          }
        
         System.out.println("NumZeroLORs"+extension+": "+countZeros);
      }
    }
}
