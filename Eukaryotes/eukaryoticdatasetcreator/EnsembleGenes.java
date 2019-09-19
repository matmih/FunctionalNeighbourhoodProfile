/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import static eukaryoticdatasetcreator.CreateReducedMappingsFile.ENCODING;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description classcontaining .dat file parser
 */
public class EnsembleGenes {
    public HashSet<String> genes=new HashSet<>();
    
    public void loadGenes(File input){
        try{
        FileInputStream inputStream = new FileInputStream(input.getAbsolutePath());
                     BufferedReader reader1 = new BufferedReader(new InputStreamReader(inputStream,ENCODING));
                      String line = null;
                      

                      while ((line = reader1.readLine()) != null) {//parse a .dat file
                        genes.add(line.trim());
                    }
                      reader1.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
    
}
