/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package semanticsimilarity;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the Wang semantic similarity between a pair of GO terms
  *https://www.bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#wang-method
 */
public class Wang {
    
    public double sA(int GOA, int reqGO,double we, Graph g){
        double res=0.0;
        
        if(reqGO==GOA)
            return 1.0;
        
        ArrayList<Integer> children=g.childList.get(reqGO);
        HashSet<Integer> inters=new HashSet<Integer>();
        Iterator<Integer> it;
        
        for(int i=0;i<children.size();i++)
            if(g.dagList.get(GOA).contains(children.get(i)))
                inters.add(children.get(i));
        
        it=inters.iterator();
        
                double temp=0.0;
        while(it.hasNext()){
			int el = it.next();
			if(el == reqGO)
                continue;
           temp=we*sA(GOA,it.next(),we,g);
           if(temp>res){
               res=temp;
           }
        }
		
        return res;
    }
    
    public double sVA(int GOA, double we, Graph g){
        double res=0.0;
        
        Iterator<Integer> it=g.dagList.get(GOA).iterator();
        
        while(it.hasNext()){
            res+=sA(GOA,it.next(),we,g);
        }
        
        return res;
    }
    
    public double sGO(int GOA, int GOB, double we, Graph g){
        double res=0.0;
        HashSet<Integer> inters=new HashSet<Integer>();
        Iterator<Integer> it=g.dagList.get(GOA).iterator();
        HashSet<Integer> lB=g.dagList.get(GOB);
        
        while(it.hasNext()){
            int elem=it.next();

            if(lB.contains(elem)){
                inters.add(elem);
            }
        }
        
        double denominator=0.0,numerator=0.0;
        
        it=inters.iterator();
        
        while(it.hasNext()){
            int elem=it.next();
            denominator+=(sA(GOA,elem,we,g)+sA(GOB,elem,we,g));
        }
        
        numerator=sVA(GOA,we,g)+sVA(GOB,we,g);
        res=denominator/numerator;
        
        return res;
    }   
}
