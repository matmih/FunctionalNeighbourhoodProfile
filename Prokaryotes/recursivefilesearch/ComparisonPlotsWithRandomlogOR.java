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
  *@description class create data for comparative plots with LOR obtained on randomised dataset, divided by semantical similarity between GO terms
 */
public class ComparisonPlotsWithRandomlogOR {
    public static void main(String arg[]){

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

        ArrayList<String> functions=new ArrayList<>();      
       HashMap<Integer,ArrayList<Double>> pvals=new HashMap<>();
          
           for(String s:cgmap.GOtoIndex.keySet()){
               if(!categories.get(0).contains(s) || !mapFull1.GOmap.containsKey(s))
                   continue;
               functions.add(s);
           }
         
        Graph g=new Graph(mapFull1.GOmap.keySet().size());
         g.createAdjacency(mapFull1);
        
          System.out.println("Adjacency graph created!");
          
          HashMap<String,HashSet<String>> GOParents=new HashMap<>();
        
          boolean writePVals=false;
          
        for(int i=0;i<functions.size();i++){
            ArrayList<Path> Paths=new ArrayList<Path>();
        semanticsimilarity.Path visited=new semanticsimilarity.Path();
        visited.path.add(mapFull1.GOmap.get(functions.get(i)));
        g.getAllPaths(mapFull1.GOmap.get(functions.get(i)), 0, Paths, visited);
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
            
          HashMap<String,ArrayList<Double>> medClHist=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHist=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHist=new HashMap<>();
          HashMap<String,ArrayList<Double>> medClHistRand=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHistRand=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHistRand=new HashMap<>();
          
        for(int iter=0;iter<functions.size();iter++){ 
            medClHist.put(functions.get(iter),new ArrayList<Double>());
            clParHist.put(functions.get(iter), new ArrayList<Double>());
            distHist.put(functions.get(iter), new ArrayList<Double>());
            medClHistRand.put(functions.get(iter),new ArrayList<Double>());
            clParHistRand.put(functions.get(iter), new ArrayList<Double>());
            distHistRand.put(functions.get(iter), new ArrayList<Double>());
            
            if(!pkGOs.contains(functions.get(iter)))
                continue;
            
     LoadOR lor=new LoadOR();
     File ORInp=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+cgmap.IndexToGO.get(iter)+".txt");
     lor.loadContValue(ORInp);
     
     LoadOR lor1=new LoadOR();
     File ORInp1=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesRandomNew\\"+"RandomContingencyOut"+cgmap.IndexToGO.get(iter)+".txt");
     lor1.loadContValue(ORInp1);
    
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
        
        ArrayList<Double> ctRand=lor1.cont.get(GOname);
        Double Randor=((ctRand.get(0)+0.5)/(ctRand.get(1)+0.5))/((ctRand.get(2)+0.5)/(ctRand.get(3)+0.5));
        Double logOrRand=Math.log10(Randor)/Math.log10(2);

        if(goSem.similarity[iGO][aInd]>2 && !GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
            medClHist.get(functions.get(iter)).add(logOr);
            medClHistRand.get(functions.get(iter)).add(logOrRand);
        }
        
        if(/*goSem.similarity[iGO-1][aInd-1]>2 && !*/GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
            clParHist.get(functions.get(iter)).add(logOr);
            clParHistRand.get(functions.get(iter)).add(logOrRand);
        }
        
        if(goSem.similarity[iGO][aInd]<=2 && !GOParents.get(functions.get(iter)).contains(GOname) /*&& pvals.get(cgmap.GOtoIndex.get(functions.get(iter))).get(goindex)<0.1*/){
            distHist.get(functions.get(iter)).add(logOr);
            distHistRand.get(functions.get(iter)).add(logOrRand);
        }
    }
    System.out.println("Iter complited: "+iter);
     }
       
       String fileName="HistDividedAll";
            try{
                String fileTmp=fileName+"ClPar.txt";
       FileWriter fw = new FileWriter(fileTmp); 
       
       for(int i=0;i<functions.size();i++){
           
         ArrayList<Double> tmp=clParHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1<tmp.size())
                                 fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
       }
       
       fw.close();
       
       fileTmp=fileName+"MedCl.txt";
       fw = new FileWriter(fileTmp);
       for(int i=0;i<functions.size();i++){
         
         
               ArrayList<Double>  tmp=medClHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
       }
       fw.close();
       
       fileTmp=fileName+"Dist.txt";
       fw = new FileWriter(fileTmp);
       
        for(int i=0;i<functions.size();i++){
                         
                ArrayList<Double> tmp=distHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
         
            
        }
         fw.close();
   
       }
            catch(Exception e)
                {e.printStackTrace();}
            
               fileName="HistDividedAllRand";
            try{
            String fileTmp=fileName+"ClPar.txt";
       FileWriter fw = new FileWriter(fileTmp);     
       
       for(int i=0;i<functions.size();i++){
           
         ArrayList<Double> tmp=clParHistRand.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                               if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
       }
       fw.close();
       
        fileTmp=fileName+"MedCl.txt";
        fw = new FileWriter(fileTmp);  
       for(int i=0;i<functions.size();i++){
         
     
               ArrayList<Double>  tmp=medClHistRand.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
       }
       fw.close();
        fileTmp=fileName+"Dist.txt";
        fw = new FileWriter(fileTmp); 
        for(int i=0;i<functions.size();i++){
                          
                ArrayList<Double> tmp=distHistRand.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1<tmp.size())
                                    fw.write(tmp.get(itt)+"\n");
                                else
                                    fw.write(tmp.get(itt)+"\n"); 
                            }
        }
       fw.close();
       }
            catch(Exception e)
                {e.printStackTrace();}  
       }        
    }
