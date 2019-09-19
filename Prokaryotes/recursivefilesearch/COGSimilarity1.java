/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;
import java.lang.Math.*;
import java.util.ArrayList;
import java.util.HashMap;
import org.javatuples.*;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class compute spatial and distance similarity
 */
 
public class COGSimilarity1 {
int[][] DistanceMatrix;
HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>> neighbours=new HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>>();
int numCOGs=0;

void computeAdjacency(COG cog){
   DistanceMatrix=new int[cog.ancogs.size()][cog.ancogs.size()];

    for(int i=0;i<cog.ancogs.size();i++)
        for(int j=0;j<cog.ancogs.size();j++){
}

void computeSpatialNeighboors(COG cog, COGGOMap mapping, GeneOGMapping geneOGMap ,int k){

    int tmpK=Math.min(k,cog.ancogs.size()-1);
    if(tmpK<=0)
        return;
    
 int shift=0;

 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

 for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(!geneOGMap.geneOGsMap.containsKey(ck.getValue2()))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         Triplet<Integer,Integer,String> tmp=cog.cogs.get(j%cog.cogs.size());

       if(!tmp.getValue2().equals("-") && (geneOGMap.geneOGsMap.keySet().contains(tmp.getValue2()))){
         neighbours.get(ck).add(cog.cogs.get(j%cog.cogs.size()));
         numCC++;
       }
         if(i-(j-i)>=0){
             tmp=cog.cogs.get(i-(j-i));
              if(!tmp.getValue2().equals("-") && (geneOGMap.geneOGsMap.keySet().contains(tmp.getValue2()))){
           neighbours.get(ck).add(cog.cogs.get(i-(j-i)));
           numCC++;
              }
         }
         else{
             tmp=cog.cogs.get(cog.cogs.size()+i-(j-i));
              if(!tmp.getValue2().equals("-") && (geneOGMap.geneOGsMap.keySet().contains(tmp.getValue2()))){
            neighbours.get(ck).add(cog.cogs.get(cog.cogs.size()+i-(j-i)));
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


void computeSpatialNeighboorsVisual(COG cog, COGGOMap mapping, GeneOGMapping geneOGMap ,int k){
 int tmpK = k;
    if(tmpK<=0)
        return;
    
 int shift=0;

 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

 for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);       
       if(!geneOGMap.geneOGsMap.containsKey(ck.getValue2()))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         Triplet<Integer,Integer,String> tmp=cog.cogs.get(j%cog.cogs.size());

         neighbours.get(ck).add(cog.cogs.get(j%cog.cogs.size()));
         numCC++;
		 
         if(i-(j-i)>=0){
             tmp=cog.cogs.get(i-(j-i));
           neighbours.get(ck).add(cog.cogs.get(i-(j-i)));
           numCC++;
         }
         else{
             tmp=cog.cogs.get(cog.cogs.size()+i-(j-i));
            neighbours.get(ck).add(cog.cogs.get(cog.cogs.size()+i-(j-i)));
            numCC++;
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

void computeSpatialNeighboorsExact(COG cog, COGGOMap mapping ,int k){
 if(cog.ancogs.size()<k+1)//maknuti ovu liniju i dodati k=velicina cog-ova da se iskoriste sve informacije
     return;

 int shift=0;

 if(k%2==0)
       shift=k/2;
 else shift=k/2+1;

 for(int i=0;i<cog.cogs.size();i++){
       Triplet<Integer,Integer,String> ck=cog.cogs.get(i);
       
       if(ck.getValue2().equals("-"))
           continue;
       if(!(mapping.CogGOmap.keySet().contains(ck.getValue2())))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
        int j=i+shift; 
         Triplet<Integer,Integer,String> tmp=cog.cogs.get(j%cog.cogs.size());

       if(!tmp.getValue2().equals("-") && (mapping.CogGOmap.keySet().contains(tmp.getValue2()))){
         neighbours.get(ck).add(cog.cogs.get(j%cog.cogs.size()));
         numCC++;
       }
         if(i-(j-i)>=0){
             tmp=cog.cogs.get(i-(j-i));
              if(!tmp.getValue2().equals("-") && (mapping.CogGOmap.keySet().contains(tmp.getValue2()))){
           neighbours.get(ck).add(cog.cogs.get(i-(j-i)));
           numCC++;
              }
         }
         else{
             tmp=cog.cogs.get(cog.cogs.size()+i-(j-i));
              if(!tmp.getValue2().equals("-") && (mapping.CogGOmap.keySet().contains(tmp.getValue2()))){
            neighbours.get(ck).add(cog.cogs.get(cog.cogs.size()+i-(j-i)));
            numCC++;
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
}

