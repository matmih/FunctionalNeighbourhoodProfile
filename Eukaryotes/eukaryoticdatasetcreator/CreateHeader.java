/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import OntologyTools.GOTerm;
import OntologyTools.GeneOntology;
import OntologyTools.GoAnnotationsSimple;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to create GO ontology header file for the final dataset (eukaryotic organisms)
 */
public class CreateHeader {
    
    public static void main(String [] args){
        
          Map<String,Set<Integer>> cogGOmapTmp=new HashMap<>();
        boolean createNewHeader=true;
        boolean createHeaderAllGO=false;
        
        GeneOntology myGO=null;
        
        try{
    myGO= new GeneOntology("C:\\Users\\matej\\Downloads\\Eukaryot data\\go_daily-termdb.obo-xml\\go_daily-termdb.obo-xml");
    }
    catch(IOException e){
        e.printStackTrace();
    }
      //  GO:0005625
       //File ogFuncFile=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_10_30NewnonCH.txt");//fungy
       File ogFuncFile=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsProp3_10_30New.txt");//metazoa
        OGGOMapping cgmap=new OGGOMapping();
        
        cgmap.createCOGGOMapping(ogFuncFile);
        
        System.out.println("go number: ");
        System.out.println(cgmap.GOtoIndex.get("GO0005624"));
        
        if(createNewHeader==true){
        Iterator<String> it =cgmap.CogGOmap.keySet().iterator();

        while(it.hasNext()){
            String cog=it.next();
            if(!cogGOmapTmp.containsKey(cog))
                cogGOmapTmp.put(cog, new HashSet<Integer>());

            ArrayList<String> gos=cgmap.CogGOmap.get(cog);

            for(String go:gos){
                if(go.contains("GO0005624"))
                System.out.println("Contained!");
                go=go.replace("GO", "");
                int goInt=Integer.parseInt(go);
                //System.out.println(goInt);
                cogGOmapTmp.get(cog).add(goInt);
            }
        }

        System.out.println("Cog size: "+cogGOmapTmp.size());

        GoAnnotationsSimple annots = new GoAnnotationsSimple(cogGOmapTmp, myGO);
        String hierarchyHeader = annots.getGoHierarchyAsClusString();
        System.out.println("Header: \n");
        System.out.println(hierarchyHeader);
        System.out.println("\n"); 

       // File out=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\headerFungy3_10_30NewnonCH.txt");//fungy
        File out=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\headerFungy3_10_30New.txt");//meatzoa

         try
            {
                FileWriter fw = new FileWriter(out);


                    fw.write(hierarchyHeader);


                fw.close();

            }
           catch(Exception e){
                e.printStackTrace();
            }
         
         for(String s:cogGOmapTmp.keySet()){
             Set<Integer> func=cogGOmapTmp.get(s);
            
             HashSet<Integer> toRemove=new HashSet<>();
             
             for(int x:func){
                 
                   if(!annots.usedGOs.contains(x))
                         toRemove.add(x);
                     }
             
             for(int x:toRemove)
                 func.remove(x);
             
             cogGOmapTmp.put(s, func);
             
         }
         
          //File out1=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\NogFunctionsProp3_10_30NewnonCHF.txt");//fungy
          File out1=new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\NogFunctionsProp3_10_30NewF.txt");//metazoa

         try
            {
                FileWriter fw = new FileWriter(out1);


                for(String s:cogGOmapTmp.keySet()){
                     Set<Integer> gos=cogGOmapTmp.get(s);
                    
                     /*if(gos.size()<5)
                         continue;*/
                     
                     fw.write(s+": ");

                   for(int h:gos){
                       
                       String g=""+h;
                       
                       while(g.length()<7)
                           g="0"+g;
                       
                       g="GO:"+g;
                       
                       fw.write(g+" ");
                   }
                    fw.write("\n");
                }

                fw.close();

            }
           catch(Exception e){
                e.printStackTrace();
            }

        }
        
    }
    
}
