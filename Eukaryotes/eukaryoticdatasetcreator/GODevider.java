/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

 /**
 * @author Vedrana Vidulin, Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to divide a GO ontology on different GO namespaces
 */
 
public class GODevider {
    
      /**
     * Extracts GO category sets.
     * 
     * @param goFile                File with gene ontology - a file ending with ".obo-xml.gz" - download the latest ontology file from the Gene ontology web page.
     *  
     * @throws java.io.IOException
     * @throws java.lang.InterruptedException
     */
    public ArrayList<HashSet<String>> ExtractCategories (File goFile, OGGOMapping map)
                                  throws IOException, InterruptedException
    {

        //Map<Integer, List<String>> res = new TreeMap<>();
       
        //----------------------------------------------------------------------
        //COLLECT GO CATEGORY TYPES FROM THE ONTOLOGY
        //----------------------------------------------------------------------
        /*Map<Integer, Double> typeMolecularFunction = new HashMap<>();
        Map<Integer, Double> typeBiologicalProcess = new HashMap<>();
        Map<Integer, Double> typeCellularComponent = new HashMap<>();*/
        
        HashSet<String> typeMolecularFunction= new HashSet<>();
        HashSet<String> typeBiologicalProcess = new HashSet<>();
        HashSet<String> typeCellularComponent = new HashSet<>();
        ArrayList<HashSet<String>> res=new ArrayList<>();
        
        InputStream fileStream = new FileInputStream(goFile);
        InputStream gzipStream = new GZIPInputStream(fileStream);
        Reader decoder = new InputStreamReader(gzipStream);
        
        try (BufferedReader bufRdr = new BufferedReader(decoder))
        {
            String line;

            Integer go = null;
            String GO="";
            
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                
                if (line.startsWith("<id>GO:"))
                    go = Integer.valueOf(line.substring(line.indexOf(":") + 1, line.lastIndexOf("<")));
                else if (line.startsWith("<namespace>") && go != null)
                {
                    String nspc = line.substring(line.indexOf(">") + 1, line.lastIndexOf("<"));
                    
                    int numDig=0;
                    int tmp=go;
                    
                    while(tmp>0){
                        tmp/=10;
                        numDig++;
                    }
                    
                    GO+="GO";
                    
                    if(numDig<8){
                        for(int i=0;i<(7-numDig);i++)
                            GO+="0";
                        GO+=go;
                    }
                    
                    if (map.GOtoIndex.containsKey(GO)){
                        System.out.println("Contained: "+GO);
                        if (nspc.equals("molecular_function"))
                            typeMolecularFunction.add(GO);
                        else if (nspc.equals("biological_process"))
                            typeBiologicalProcess.add(GO);
                        else if (nspc.equals("cellular_component"))
                            typeCellularComponent.add(GO);
                    }
                    //System.out.println("GO: "+GO);
                    //System.out.println("GO num: "+go);
                    go = null;
                    GO="";
                }
            }
        }
        res.add(typeBiologicalProcess); res.add(typeMolecularFunction); res.add(typeCellularComponent);
        return res;
    }
}
