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
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description main class creating data for comparative LOR histogram divided by semantic similarity  
 */
public class HistogramDataForSelectGOs {
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
        cgmap.reduceMap(red);
        geneOGMap.reduceMap(red);
        
        semanticsimilarity.GOMap mapFull=new semanticsimilarity.GOMap();
         mapFull.CreateGOMap(cgmap);
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
         
          String path="";

        ArrayList<String> functions=new ArrayList<>();      
          
        for(int iter=0;iter<cgmap.GOtoIndex.keySet().size();iter++){ 
     
            ArrayList<Double> tmp=new ArrayList<>(Collections.nCopies(cgmap.GOtoIndex.keySet().size(), 0.0));
        }
             
        functions.add("GO0005975");
        functions.add("GO0006260");
        functions.add("GO0006974");
        functions.add("GO0016051");
        functions.add("GO0006457");
        functions.add("GO0046700");
        functions.add("GO0006310");
        functions.add("GO0046903");
        functions.add("GO0008610");
        functions.add("GO0006865");
        functions.add("GO0006974");
        functions.add("GO0006508");
        
          HashMap<String,HashSet<String>> GOParents=new HashMap<>();
        
          boolean writePVals=false;
          
        for(int i=0;i<functions.size();i++){
        HashSet<String> gos=new HashSet<>();
        gos.add(functions.get(i));
        
        int go1=rgt.translateToIntSingle(functions.get(i));
        System.out.println("provjera: "+go1+" "+functions.get(i));
        HashSet<Integer> pars=new HashSet<>();
        Collection<GOTerm> parents = myGO.get(go1).getAllParents();
    
	for ( GOTerm curPar : parents ) {
         if (curPar.getId() != 8150 && curPar.getId() != 3674 && curPar.getId() != 5575)
           pars.add(curPar.getId());
    }
    
    HashSet<String> parsString=rgt.translateToString(pars);
    gos.addAll(parsString);
        
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
          
          HashMap<String,ArrayList<Double>> medClHist=new HashMap<>();
          HashMap<String,ArrayList<Double>> clParHist=new HashMap<>();
          HashMap<String,ArrayList<Double>> distHist=new HashMap<>();
          
        for(int iter=0;iter<functions.size();iter++){ 
            medClHist.put(functions.get(iter),new ArrayList<Double>());
            clParHist.put(functions.get(iter), new ArrayList<Double>());
            distHist.put(functions.get(iter), new ArrayList<Double>());
            
            ArrayList<String> tmp=new ArrayList<>();
          
		  if(!pkGOs.contains(functions.get(iter)))
                continue;
            
     LoadOR lor=new LoadOR();
     File ORInp=new File("C:\\Users\\matej\\Desktop\\Contingency tables\\AllTablesNew\\"+"ContingencyOut"+functions.get(iter)+".txt");
     lor.loadContValue(ORInp);
            
      ArrayList<Double> ctSelf=lor.cont.get(functions.get(iter));
    Double orSelf=((ctSelf.get(0)+0.5)/(ctSelf.get(1)+0.5))/((ctSelf.get(2)+0.5)/(ctSelf.get(3)+0.5));
    Double logOrSelf=Math.log10(orSelf)/Math.log10(2);
     
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

        if(goSem.similarity[iGO][aInd]>2 && !GOParents.get(functions.get(iter)).contains(GOname)){
            medClHist.get(functions.get(iter)).add(logOr);
        }
        
        if(GOParents.get(functions.get(iter)).contains(GOname)){
            clParHist.get(functions.get(iter)).add(logOr);
        }
        
        if(goSem.similarity[iGO][aInd]<=2 && !GOParents.get(functions.get(iter)).contains(GOname)){
            distHist.get(functions.get(iter)).add(logOr);
        }
    }
     }
       
       String fileName="HistDivided";
       
       for(int i=0;i<functions.size();i++){
           String fileNameTmp=fileName+functions.get(i)/*+".txt"*/;
           try{
            String fnTmp=fileNameTmp+"CLPar.txt";
         FileWriter fw = new FileWriter(fnTmp); 
         
         ArrayList<Double> tmp=clParHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1!=tmp.size())
                                    fw.write(tmp.get(itt)+",");
                                if(itt+1==tmp.size()){
                                    fw.write(tmp.get(itt)+"");
                                    fw.write("\n");
                                }
                            }
         
         fw.close();
         fnTmp=fileNameTmp+"MedCl.txt";
         fw = new FileWriter(fnTmp); 
                tmp=medClHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1!=tmp.size())
                                    fw.write(tmp.get(itt)+",");
                                if(itt+1==tmp.size()){
                                    fw.write(tmp.get(itt)+"");
                                    fw.write("\n");
                                }
                            }
                            
         fw.close();
         fnTmp=fileNameTmp+"Dist.txt";
         fw = new FileWriter(fnTmp); 
                tmp=distHist.get(functions.get(i));
         
                            for(int itt=0;itt<tmp.size();itt++){
                                if(itt+1!=tmp.size())
                                    fw.write(tmp.get(itt)+",");
                                if(itt+1==tmp.size()){
                                    fw.write(tmp.get(itt)+"");
                                    fw.write("\n");
                                }
                            }
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}
       }         
       }
        
    }
