/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import javafx.util.Pair;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to compute a gene sequence from the gene data and the sequence data
 */
public class ComputeGeneSequence {
    HashMap<String,String> geneSequences;
    
    public ComputeGeneSequence(){
        geneSequences = new HashMap<>();
    }
    
    public void computeGeneSeqInfo(LoadGeneData g, LoadSequenceData l){
        
        Iterator<String> it = g.genes.keySet().iterator();
        
        while(it.hasNext()){
            String gid = it.next();
            
            if(geneSequences.containsKey(gid))
                continue;
            
            Pair<ArrayList<Integer>,String> gInfo = g.genes.get(gid);
            int minLoc = gInfo.getKey().get(0);
            int maxLoc = gInfo.getKey().get(1);
            int strand = gInfo.getKey().get(2);
            String cont = gInfo.getValue();
            if(l.contigSequence.get(cont).length()<=(maxLoc-1)){
                System.out.println("******************************");
                  System.out.println("seq size: "+l.contigSequence.get(cont).length());
                  System.out.println("gid: "+gid);
                  System.out.println("max: "+maxLoc);
                  System.out.println("min: "+minLoc);
                 System.out.println("******************************");
                 continue;
            }
            String geneSeq = l.contigSequence.get(cont).substring(minLoc-1,maxLoc);
            if(strand==0){
                geneSeq = new StringBuilder(geneSeq).reverse().toString();
                StringBuilder complement = new StringBuilder(geneSeq);
                for(int k=0;k<complement.length();k++){
                    if(complement.charAt(k)=='A')
                        complement.setCharAt(k, 'T');
                   else if(complement.charAt(k)=='T')
                        complement.setCharAt(k, 'A');
                   else if(complement.charAt(k)=='G')
                        complement.setCharAt(k, 'C');
                   else if(complement.charAt(k)=='C')
                        complement.setCharAt(k, 'G');
                    
                }
                
                geneSeq = complement.toString();
                
            }
            geneSequences.put(gid, geneSeq);           
        }
    }
    
    public void shortOutput(){
        
        Iterator<String> it = geneSequences.keySet().iterator();
        int count = 0;
        while(it.hasNext()){
            String gid = it.next();
            
            String seq = geneSequences.get(gid);
            
            if(seq.length()>20)
                System.out.println(gid+" : "+seq.substring(0,20));
            else System.out.println(seq);
            
            if(count>10)
                return;
            count++;
            
        }
        
    }
    
    
}
