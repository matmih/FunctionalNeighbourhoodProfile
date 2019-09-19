/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.File;
import org.javatuples.Pair;
 
 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to add index for use with k-NN
 * input: distanceForest dataset
 * replace \t with space
 */
 
public class CreateBBDataTable {
    
    static public void main(String args[]){
        ARFF fRead=new ARFF();
        fRead.loadARFF(new File("DistanceForest.arff"));
        fRead.attributes.add(1,new Pair("index","numeric"));
        
        for(int i=0;i<fRead.data.size();i++)
            fRead.data.get(i).add(0, (double)i);
        
        fRead.writeARFFHierarchy("DistancekNNGFA.arff", true);
        
    }
    
}
