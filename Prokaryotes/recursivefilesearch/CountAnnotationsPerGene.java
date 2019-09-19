/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.File;
import OntologyTools.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to count the number of annotations per gene
 */
 
public class CountAnnotationsPerGene {

    static public void main(String args[]){

        //for windows
        //File geneOgs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\pid2ogs.txt");
        //File OgGOs=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");

        //for linux
       File geneOgs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/pid2ogs.txt");
       File OgGOs=new File("/home/mmihelcic/GOdatasetGenerator/GOanotacije/og2gos-Uniprot-GOA-10-12-2013-2070_orgs-OG_has_funcs_that_50_perc_or_more_genes_with_functions_in_that_OG_have.txt");
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
        
        String mappingOut = "ReducedCogOGMappingPhenotypeR.txt";
        
        boolean writecogogmapping = false;
        
        File outMPP = new File(mappingOut);
      
        if(writecogogmapping )
        try
            {
                FileWriter fw = new FileWriter(outMPP);


                   Iterator<String> it = cgmap.CogGOmap.keySet().iterator();
                   while(it.hasNext()){
                       String og = it.next();
                       if(og.contains("NOG"))
                           continue;
                       ArrayList<String> gos = cgmap.CogGOmap.get(og);
                       
                       fw.write(og+":");
                       
                       for(int i=0;i<gos.size();i++)
                           if(i+1<gos.size())
                                fw.write(gos.get(i)+",");
                           else fw.write(gos.get(i)+"\n");
                       
                   }

                fw.close();
                    return;
            }
           catch(Exception e){
                e.printStackTrace();
            }
        
        geneOGMap.reduceMap(red);

        Map<String,Set<Integer>> cogGOmapTmp=new HashMap<>();
        boolean createNewHeader=false;
        boolean createHeaderAllGO=false;
        
        
        if(createHeaderAllGO==true){
            COGGOMap cgmap1=new COGGOMap();
        cgmap1.createCOGGOMapping(OgGOs);
            Iterator<String> it =cgmap.CogGOmap.keySet().iterator();

        while(it.hasNext()){
            String cog=it.next();
            if(!cogGOmapTmp.containsKey(cog))
                cogGOmapTmp.put(cog, new HashSet<Integer>());

            ArrayList<String> gos=cgmap.CogGOmap.get(cog);

            for(String go:gos){
                go=go.replace("GO", "");
                int goInt=Integer.parseInt(go);
                cogGOmapTmp.get(cog).add(goInt);
            }
        }
        
        GoAnnotationsSimple annots = new GoAnnotationsSimple(cogGOmapTmp, myGO);
        String hierarchyHeader = annots.getGoHierarchyAsClusString();
        System.out.println("Header: \n");
        System.out.println(hierarchyHeader);
        System.out.println("\n");

        File out=new File("headerFull.txt");

         try
            {
                FileWriter fw = new FileWriter(out);


                    fw.write(hierarchyHeader);


                fw.close();

            }
           catch(Exception e){
                e.printStackTrace();
            }
         
        return;
           
        }
        
        if(createNewHeader==true){
        Iterator<String> it =cgmap.CogGOmap.keySet().iterator();

        while(it.hasNext()){
            String cog=it.next();
            if(!cogGOmapTmp.containsKey(cog))
                cogGOmapTmp.put(cog, new HashSet<Integer>());

            ArrayList<String> gos=cgmap.CogGOmap.get(cog);

            for(String go:gos){
                go=go.replace("GO", "");
                int goInt=Integer.parseInt(go);
                cogGOmapTmp.get(cog).add(goInt);
            }
        }

        System.out.println("Cog size: "+cogGOmapTmp.size());

        GoAnnotationsSimple annots = new GoAnnotationsSimple(cogGOmapTmp, myGO);
        String hierarchyHeader = annots.getGoHierarchyAsClusString();
        System.out.println("Header: \n");
        System.out.println(hierarchyHeader);
        System.out.println("\n");

        File out=new File("header.txt");

         try
            {
                FileWriter fw = new FileWriter(out);


                    fw.write(hierarchyHeader);


                fw.close();

            }
           catch(Exception e){
                e.printStackTrace();
            }
        }

        String linuxInputFolder="/home/mmihelcic/GOdatasetGenerator/GOanotacije/tmpFiles/proba";
        String[] extensions = {"ptt"};
           String BaselineLinuxOoutputFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/BaselineOrganism1.arff";
           String BaselineLinuxOoutputFileNoOperons="/home/mmihelcic/GOdatasetGenerator/GOanotacije/BaselineOrganismNoOps1.arff";
           String DistanceLinuxOoutputFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/DistanceForest.arff";
           String BaselineLinuxGroupFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/GroupOrganism1.arff";
           String LinuxCombinedFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/CombinedOrganism.arff";
           String taxIDFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/TaxIDs.txt";
           String cogCounts="CogNogCounts.txt";
           String GOCounts="GOCounts.txt";
           String contingencyOut="/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContAll/ContingencyOut";
           String contingencyOutRand="/home/mmihelcic/GOdatasetGenerator/GOanotacije/NewContRandAll/RandomContingencyOut";
           String AssociationFile="/home/mmihelcic/GOdatasetGenerator/GOanotacije/NetworkPropagation.arff";
           String accretionFile = "/home/mmihelcic/GOdatasetGenerator/Recursive file search/PredictionsAtPrecision/IA.csv";
           String AccrtetionFileOutput = "/home/mmihelcic/GOdatasetGenerator/Recursive file search/PredictionsAtPrecision/NetworkPropagation.arff";
           
         boolean recursive = true;
            int k=4, intTrain=1, createDatasets=7,createContingency=0,createGFStatistics=0;

        MultiGenomeSet mgs=new MultiGenomeSet(linuxInputFolder);      
        
         ArrayList<Double> res=mgs.numGenes(extensions, recursive, k, cgmap, geneOGMap, 1, AssociationFile, intTrain);
        
         System.out.println("NumGenes: "+res.get(0));
         System.out.println("NumAnnots: "+res.get(1));
         System.out.println("Annots per gene: "+(res.get(1))/res.get(0));    
    }
}
