/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.util.ArrayList;
import java.util.HashMap;
import org.javatuples.Triplet;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the gene neighbourhood
 */
public class GeneNeighbours {

int[][] DistanceMatrix;
HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>> neighbours=new HashMap<Triplet<Integer,Integer,String>,ArrayList<Triplet<Integer,Integer,String>>>();

 int numCOGs=0;
void computeSpatialNeighboors(Gene genes, int k){
    int tmpK=Math.min(k,genes.anotGen.size()-1);
    if(tmpK<=0)
        return;
    
 int shift=0;
 
 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

 
 for(int i=0;i<genes.genes.size();i++){
       Triplet<Integer,Integer,String> ck=genes.genes.get(i);

       if(!genes.anotGen.contains(ck.getValue2()))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         
      
         Triplet<Integer,Integer,String> tmp=null;
         
         if(j<genes.genes.size()){
                tmp=genes.genes.get(j);

       if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
         neighbours.get(ck).add(genes.genes.get(j));
         numCC++;
       }
     }
         if(i-(j-i)>=0){
             tmp=genes.genes.get(i-(j-i));
              if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
           neighbours.get(ck).add(genes.genes.get(i-(j-i)));
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


void computeSpatialNeighboorsDeep(Gene genes, int k){
    int tmpK=Math.min(k,genes.anotGen.size()-1);
    if(tmpK<=0)
        return;
    
 int shift=0;
 
 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

 
 for(int i=0;i<genes.genes.size();i++){
       Triplet<Integer,Integer,String> ck=new Triplet<>(genes.genes.get(i).getValue0(),genes.genes.get(i).getValue1(),genes.genes.get(i).getValue2());

       if(!genes.anotGen.contains(ck.getValue2()))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(new Triplet<>(ck.getValue0(),ck.getValue1(),ck.getValue2()), elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         
      
         Triplet<Integer,Integer,String> tmp=null;
         
         if(j<genes.genes.size()){
                tmp = new Triplet<>(genes.genes.get(j).getValue0(),genes.genes.get(j).getValue1(),genes.genes.get(j).getValue2());

       if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
         neighbours.get(ck).add(/*genes.genes.get(j)*/tmp);
         numCC++;
       }
     }
         if(i-(j-i)>=0){
             tmp = new Triplet<>(genes.genes.get(i-(j-i)).getValue0(),genes.genes.get(i-(j-i)).getValue1(),genes.genes.get(i-(j-i)).getValue2());

              if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
           neighbours.get(ck).add(/*genes.genes.get(i-(j-i))*/tmp);
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


void computeSpatialNeighboorsnonCH(Gene genes, int k){
    int tmpK=Math.min(k,genes.anotGen.size()-1);
    if(tmpK<=0)
        return;
    
 int shift=0;
 
 if(tmpK%2==0)
     shift=tmpK/2;
 else shift=tmpK/2+1;

 
 for(int g=0;g<genes.genesContig.size();g++){
     genes.genes = genes.genesContig.get(g);
     
 for(int i=0;i<genes.genes.size();i++){
       Triplet<Integer,Integer,String> ck=genes.genes.get(i);

       if(!genes.anotGen.contains(ck.getValue2()))
           continue;
       
       if(!neighbours.containsKey(ck)){
          ArrayList<Triplet<Integer,Integer,String>> elem=new  ArrayList<Triplet<Integer,Integer,String>>();
          neighbours.put(ck, elem);
       }
int numCC=0;
     for(int j=i+1;j<i+shift+1;j++){
         
      
         Triplet<Integer,Integer,String> tmp=null;
         
         if(j<genes.genes.size()){
                tmp=genes.genes.get(j);

       if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
         neighbours.get(ck).add(genes.genes.get(j));
         numCC++;
       }
     }
         if(i-(j-i)>=0){
             tmp=genes.genes.get(i-(j-i));
              if(!tmp.getValue2().equals("-") && (genes.anotGen.contains(tmp.getValue2()))){
           neighbours.get(ck).add(genes.genes.get(i-(j-i)));
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
}

}
