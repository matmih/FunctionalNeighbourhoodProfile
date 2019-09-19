/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package semanticsimilarity;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about GOSimilarity
 */
public class GOSimilarity {
    double similarity[][];
    double similarityLin[][];
    double similarityRel[][];
    double similarityWang[][];
    int numElem;
    
    public GOSimilarity(int numGOs){
        similarity=new double[numGOs][numGOs];
        similarityLin=new double[numGOs][numGOs];
        similarityRel=new double[numGOs][numGOs];
        similarityWang=new double[numGOs][numGOs];
        numElem=numGOs;
    }
    
    public void addSimilarity(int GO1, int GO2,double[] sim){
      
        similarity[GO1-1][GO2-1]=sim[0];
        similarityLin[GO1-1][GO2-1]=sim[1];
        similarityRel[GO1-1][GO2-1]=sim[2];
          if(GO1!=GO2){
        similarity[GO2-1][GO1-1]=sim[0];
        similarityLin[GO2-1][GO1-1]=sim[1];
        similarityRel[GO2-1][GO1-1]=sim[2];
          }
    }
    
    public void addWangSimilarity(int GO1, int GO2, double sim){
         similarityWang[GO1-1][GO2-1]=sim;
         if(GO1!=GO2)
              similarityWang[GO2-1][GO1-1]=sim;
    }
    
    ArrayList<Double> getSimilarities(int GO){
        ArrayList<Double> l=new ArrayList<Double>();
        
        for(int i=0;i<numElem;i++)
            l.add(similarity[GO][i]);
        
        return l;
    }
    
    void writeSimilarity(File[] output){
          try
            {
                FileWriter fw = new FileWriter(output[0]);
                
                for(int i=0;i<numElem;i++)
                    for(int j=0;j<numElem;j++)
                        if((j+1)<numElem)
                        fw.write(similarity[i][j]+" ");
                        else
                            fw.write(similarity[i][j]+"\n");
                fw.close();
                
                
                fw=new FileWriter(output[1]);
                
                 for(int i=0;i<numElem;i++)
                    for(int j=0;j<numElem;j++)
                        if((j+1)<numElem)
                        fw.write(similarityLin[i][j]+" ");
                        else
                            fw.write(similarityLin[i][j]+"\n");
                fw.close();
                
                
                fw=new FileWriter(output[2]);
                
                 for(int i=0;i<numElem;i++)
                    for(int j=0;j<numElem;j++)
                        if((j+1)<numElem)
                        fw.write(similarityRel[i][j]+" ");
                        else
                            fw.write(similarityRel[i][j]+"\n");
                fw.close();
                
                 fw=new FileWriter(output[3]);
                
                 for(int i=0;i<numElem;i++)
                    for(int j=0;j<numElem;j++)
                        if((j+1)<numElem)
                        fw.write(similarityWang[i][j]+" ");
                        else
                            fw.write(similarityWang[i][j]+"\n");
                fw.close();
                
            }
           catch(Exception e){
                e.printStackTrace();
            }
    }
    
}
