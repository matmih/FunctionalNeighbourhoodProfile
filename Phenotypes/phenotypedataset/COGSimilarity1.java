/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypedataset;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.javatuples.*;
/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to compute spatial and distance similarity
 */
public class COGSimilarity1 {
int[][] DistanceMatrix;
HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>> neighbours=new HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>>();
int numCOGs=0;
public int annotatedGenes = 0, numgenes=0;

void computeAdjacency(COG cog){
   DistanceMatrix=new int[cog.ancogs.size()][cog.ancogs.size()];

    for(int i=0;i<cog.ancogs.size();i++)
        for(int j=0;j<cog.ancogs.size();j++){}
           // DistanceMatrix[i][j]=Math.min((Integer)cog.ancogs.get(j).getValue0()-(Integer)cog.ancogs.get(i).getValue1(),(Integer)cog.cogs.get(cog.cogs.size()-1).getValue+(Integer)cog.ancogs.get(i).getValue0());
}

void countAnnotatedGenes(LoadGeneData1 genes, Mappings map){
    
    for(int kcont=0;kcont<genes.geneContigs.size();kcont++){ 
        ArrayList<Triplet<Integer,Integer,String>> t = genes.geneContigs.get(kcont);
        for(int i=0;i<t.size();i++){
            System.out.println("Cont size: "+t.size());
            Triplet<Integer,Integer,String> ck=t.get(i);
            if(!map.geneOgMapping.containsKey(ck.getValue2())){
                numgenes++;
           continue;
       }
            else{
                numgenes++;
                annotatedGenes++;
            }
        }
    }
    
}


void computeSpatialNeighboors(LoadGeneData1 genes, Mappings map, int k){

    int tmpK=k;
    if(tmpK<=0)
        return;
    
 int shift=0;
 
 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

for(int kcont=0;kcont<genes.geneContigs.size();kcont++){ 
    ArrayList<Triplet<Integer,Integer,String>> t = genes.geneContigs.get(kcont);
    
 for(int i=0;i<t.size();i++){
     System.out.println("Cont size: "+t.size());
     if(t.size()<k)
         continue;
       Triplet<Integer,Integer,String> ck=t.get(i);

       if(!map.geneOgMapping.containsKey(ck.getValue2())){
           continue;
       }

       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         Triplet<Integer,Integer,String> tmp=t.get(j%t.size());

       if(!tmp.getValue2().equals("-") && (map.geneOgMapping.keySet().contains(tmp.getValue2()))){
         neighbours.get(ck).add(t.get(j%t.size()));
         numCC++;
       }
         if(i-(j-i)>=0){
             tmp=t.get(i-(j-i));
              if(!tmp.getValue2().equals("-") && (map.geneOgMapping.keySet().contains(tmp.getValue2()))){
           neighbours.get(ck).add(t.get(i-(j-i)));
           numCC++;
              }
         }
         else{
             tmp=t.get(t.size()+i-(j-i));
              if(!tmp.getValue2().equals("-") && (map.geneOgMapping.keySet().contains(tmp.getValue2()))){
            neighbours.get(ck).add(t.get(t.size()+i-(j-i)));
            numCC++;
              }
         }
     }
       if(numCC>k)
           neighbours.get(ck).remove(neighbours.get(ck).size()-1);
       if(neighbours.get(ck).size()>0)
           numCOGs++;
       else if(neighbours.get(ck).size()==0)
           neighbours.remove(ck);
    }
  }
        System.out.println("neighbours size: "+neighbours.size());
        
        
}
}

