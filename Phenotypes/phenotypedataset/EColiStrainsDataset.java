/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypedataset;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  * work was performed while MM was a visiting Phd student at IRB Barcelona
  *@mail matmih1@gmail.com
  *@description class to compute bacterial gene sequence
 */
public class EColiStrainsDataset {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {        
        String genomesPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains genomes\\Bacteria";
        String sequencesPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains sequences\\Sequences";
        String outputPath = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains sequences\\bactSeq.fa";
        
       
        File inputGenomes = new File(genomesPath);
        File inputSequences = new File(sequencesPath);
        File output = new File(outputPath);
        String extension[] = {"gff"};
        
         ComputeGeneSequence gs = new ComputeGeneSequence();
         CreateFastaFile cff = new CreateFastaFile(output);
         String pathT = "C:\\Users\\mmihelcic\\Downloads\\Bacterial strains sequences\\Sequences\\";
         
         int count=0, it=Integer.MAX_VALUE;
         
         try{
         Collection files = FileUtils.listFiles(inputGenomes, extension, true);

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
        LoadGeneData g = new LoadGeneData();
        LoadSequenceData s = new LoadSequenceData();
     
        
        g.loadData(input);
        g.shortOutput();
        
        System.out.println("Gene output finished.....");
        
        String name = input.getName();
        String t[] = name.split("\\.");
        name = t[0]+".fasta";
        String path = pathT+name;
        
        s.loadData(new File(path));
        s.shortOutput();
        
        System.out.println("Sequence output finished.....");
        
        gs.computeGeneSeqInfo(g, s);
        gs.shortOutput();
        count++;
        if(count>it)
            break;
            }
        
                } catch (Exception e) {
            e.printStackTrace();
        }
         
         System.out.println("Gene sequence number: "+gs.geneSequences.keySet().size());
         
         
         int batchNum = 0, numPerFile = 250000, added=0;
         ArrayList<String> ordered = new ArrayList<>(numPerFile/*gs.geneSequences.keySet().size()*/);
        
         Iterator<String> itOrd = gs.geneSequences.keySet().iterator();
         
         while(itOrd.hasNext()){
             if(ordered.size()==0)
                 ordered.add(itOrd.next());
             else{
                 String key = itOrd.next();
                 for(int i=0;i<ordered.size();i++){
                     if(gs.geneSequences.get(ordered.get(i)).length()<gs.geneSequences.get(key).length()){
                         ordered.add(i, key);
                         added=1;
                         if(ordered.size()%100 ==0)
                             System.out.println("ordered size: "+ordered.size());
                         break;
                     }
                 }
                 if(added==0){
                     ordered.add(ordered.size(), key);
                 }
                 else added=0;
             }
             
             if(ordered.size() == numPerFile){
                 cff.createFileParts(ordered, gs, numPerFile, batchNum);
                 batchNum++;
                 ordered.clear();
             }
             
         }
        
         cff.createFileParts(ordered, gs, numPerFile, batchNum);
         
        System.out.println("Gene sequences computed.....");
    }
    
}
