/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to create FASTA file from the gene sequence
 */
public class CreateFastaFile {
    File output;
    
    public CreateFastaFile(File out){
        output = out;
    }
    
    
    public void createFile(ComputeGeneSequence cg){
        
        try{
 
            FileWriter fw = new FileWriter(output); 
        
            Iterator<String> it=cg.geneSequences.keySet().iterator();
            
            while(it.hasNext()){
                   String gid = it.next();
                   String seq = cg.geneSequences.get(gid).replaceAll("[Y,N,M,R,S,K,W,B,D,H,V]", "");//delete Ns from output
                   fw.write(">"+gid+"\n");
                   for(int i=0;i<((seq.length()/60)+1);i++){
                       if(seq.length()>(60*(i+1)))
                            fw.write(seq.substring(i*60,(i+1)*60)+"\n");
                       else fw.write(seq.substring(i*60,seq.length())+"\n");
                   }
          }

        fw.close();

   }catch (Exception e) {
            e.printStackTrace();
        }       
    }
    
    public void createFileParts(ComputeGeneSequence cg, int numPerFile){
        
        try{
 
            int numSeq = cg.geneSequences.keySet().size();
            int numIter = (int)Math.ceil((double)numSeq/numPerFile);
            
            Object gids[] = cg.geneSequences.keySet().toArray();
            
          for(int iter=0;iter<numIter;iter++){ 
              String outPart = output.getAbsolutePath().split("\\.")[0]+iter+".fa";
            FileWriter fw = new FileWriter(new File(outPart)); 
        

           for(int j=0;j<numPerFile;j++){
               if(iter*numPerFile+j>=gids.length)
                   break;
                   
                   
                   String gid = (String)gids[iter*numPerFile+j];
                   String seq = cg.geneSequences.get(gid).replaceAll("[Y,N,M,R,S,K,W,B,D,H,V]", "");//delete Ns from output
                   fw.write(">"+gid+"\n");
                   for(int i=0;i<((seq.length()/60)+1);i++){
                       if(seq.length()>(60*(i+1)))
                            fw.write(seq.substring(i*60,(i+1)*60)+"\n");
                       else fw.write(seq.substring(i*60,seq.length())+"\n");
                   }
          }

        fw.close();
          }
   }catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    public void createFileParts(ArrayList<String> order, ComputeGeneSequence cg, int numPerFile){
        
        try{
 
            int numSeq = cg.geneSequences.keySet().size();
            int numIter = (int)Math.ceil((double)numSeq/numPerFile);
                        
          for(int iter=0;iter<numIter;iter++){ 
              String outPart = output.getAbsolutePath().split("\\.")[0]+iter+".fa";
            FileWriter fw = new FileWriter(new File(outPart)); 
        

           for(int j=0;j<numPerFile;j++){
               if(iter*numPerFile+j>=order.size())
                   break;
                   
                   
                  String gid = order.get(iter*numPerFile+j);
                   String seq = cg.geneSequences.get(gid).replaceAll("[Y,N,M,R,S,K,W,B,D,H,V]", "");//delete Ns from output
                   fw.write(">"+gid+"\n");
                   for(int i=0;i<((seq.length()/60)+1);i++){
                       if(seq.length()>(60*(i+1)))
                            fw.write(seq.substring(i*60,(i+1)*60)+"\n");
                       else fw.write(seq.substring(i*60,seq.length())+"\n");
                   }
          }

        fw.close();
          }
   }catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    public void createFileParts(ArrayList<String> order, ComputeGeneSequence cg, int numPerFile, int BatchNum){
        
        try{
 
            int numSeq = cg.geneSequences.keySet().size();
            int numIter = 1;//(int)Math.ceil((double)numSeq/numPerFile);
            
          for(int iter=0;iter<numIter;iter++){ 
              String outPart = output.getAbsolutePath().split("\\.")[0]+BatchNum+".fa";
            FileWriter fw = new FileWriter(new File(outPart)); 
        

           for(int j=0;j<numPerFile;j++){
               if(iter*numPerFile+j>=order.size())
                   break;
                   
                  String gid = order.get(iter*numPerFile+j);
                   String seq = cg.geneSequences.get(gid).replaceAll("[Y,N,M,R,S,K,W,B,D,H,V]", "");//delete Ns from output
                   fw.write(">"+gid+"\n");
                   for(int i=0;i<((seq.length()/60)+1);i++){
                       if(seq.length()>(60*(i+1)))
                            fw.write(seq.substring(i*60,(i+1)*60)+"\n");
                       else if(i*60<seq.length()) fw.write(seq.substring(i*60,seq.length())+"\n");
                   }
          }

        fw.close();
          }
   }catch (Exception e) {
            e.printStackTrace();
        }
    }
}
