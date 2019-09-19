/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import org.apache.commons.io.FileUtils;
import org.javatuples.Pair;
import org.javatuples.Triplet;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing functions to create contingency tables and datasets for different models (eukaryotic organisms)
 */
public class ConstructionAlgorithms {
     final static Charset ENCODING = StandardCharsets.UTF_8;
     
     File root;

    ConstructionAlgorithms(String rootPath){
         root=new File(rootPath);
    }
    
     void createContingency(String taxIDFilePath,String[] extensions, boolean recursive, int k,OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, Mappings map ,String output,int isTrain, boolean randomize , int OrgType,HashMap<Integer,Integer> taxTranslation){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
         HashSet<Integer> okTaxes = new HashSet<>();
          if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
         Random rand = new Random();
         
          int taxid = -1;
         
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);

         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);//nebitno
         cot.initializeContingency(cgmap);
         int numGenomes=0;
         int batchSize=30, batchCount=0;
         final ArrayList<GeneNeighbours> simList=new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                
                //add all cases and files for fungy, metazoa + parameters
                
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                 int found=0;
                   String fileName=input.getName().substring(0,input.getName().length()-4);
                   taxid=taxIDMap.get(fileName);
                 
                 if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){
                  Gene cg=new Gene();
                  if(OrgType==0)
                     cg.findCogsNew(input,geneOGMap,map);//fungy
                  else if(OrgType ==1)
                       cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);
                 
                  
                  if(randomize == true)
                      cg.randomize(rand);
                  
                 GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                 

                if(sim.numCOGs==0){//unsesn bias u OR-ove
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                 
                
                simList.add(sim);
                 batchCount++;
                 }
                 else{
                     continue;
                 }  

                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final OGGOMapping c=cgmap;
                    final HashMap<String,HashSet<Pair<String,String>>> cGO=geneOGMap;
                    final HashMap<Integer,Integer> taxTranslationF  = taxTranslation;
                    final int taxidF = taxid;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_ind = index;
                        @Override
                        public void run() {
                         Contingency cotTmp=new Contingency();
                         cotTmp.createGOCog(c);
                         cotTmp.initializeContingency(c);
                         cotTmp.addContingency(simList.get(proc_ind),c,cGO,taxidF,taxTranslationF);
                         
                                  try{    
                                    lock.lock();
                             try{
                                 contList.set(proc_ind,cotTmp);
                             }finally{
                                 lock.unlock();
                             }
                        }
                        catch(Exception e){
                            e.printStackTrace();
                         }                             
                        }
                          });
                    }
                    }catch(Exception e){
                        e.printStackTrace();
                    }
                finally {
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                    
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);

            }
            
            if(batchCount!=0){
                batchSize = batchCount;
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final OGGOMapping c=cgmap;
                    final HashMap<String,HashSet<Pair<String,String>>> cGO=geneOGMap;
                     final HashMap<Integer,Integer> taxTranslationF  = taxTranslation;
                    final int taxidF = taxid;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                        @Override
                        public void run() {
                            Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(c);
                        cotTmp.initializeContingency(c);
                        cotTmp.addContingency(simList.get(index),c,cGO,taxidF,taxTranslationF);
                        try{    
                        lock.lock();
                             try{
                                 contList.set(index,cotTmp);
                             }finally{
                                 lock.unlock();
                             }
                        }
                        catch(Exception e){
                            e.printStackTrace();
                         }
                        }
                          });
                    }
                    }catch(Exception e){
                        e.printStackTrace();
                    }
                finally {
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
            }
            
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
     
     
    void createContingencyNoSinc(String taxIDFilePath,String[] extensions, boolean recursive, int k,OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap, Mappings map ,String output,int isTrain, boolean randomize , int OrgType,HashMap<Integer,Integer> taxTranslation){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
         HashSet<Integer> okTaxes = new HashSet<>();
          if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
         Random rand = new Random();
         
          int taxid = -1;
         
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);

         Contingency cot=new Contingency();
         cot.createGOCog(cgmap);//nebitno
         cot.initializeContingency(cgmap);
         int numGenomes=0, numLocs = 0;
         int batchSize=30, batchCount=0;
         final ArrayList<GeneNeighbours> simList=new ArrayList<>();
         final ArrayList<Integer> taxesUsed = new ArrayList<>();
          final ArrayList<Contingency> contList=new ArrayList<>();
          ArrayList<Gene> genes = new ArrayList<>();
          int termination=50;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                               
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                 int found=0;
                   String fileName=input.getName().substring(0,input.getName().length()-4);
                   taxid=taxIDMap.get(fileName);
                 
                 if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){
                  Gene cg=new Gene();
                  if(OrgType==0)
                     cg.findCogsNew(input,geneOGMap,map);//fungy
                  else if(OrgType ==1)
                       cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);
                 
                  
                  if(randomize == true)
                      cg.randomize(rand);
                  
                 GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                 

                if(sim.numCOGs==0){//unsesn bias u OR-ove
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                 
                
                taxesUsed.add(taxid);
                simList.add(sim);
                genes.add(cg);
                numLocs+=sim.neighbours.keySet().size();
                 batchCount++;
                 }
                 else{
                     continue;
                 }  

                if(batchCount==(batchSize)){
                    
                    for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(30);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final OGGOMapping c=cgmap;
                    final HashMap<String,HashSet<Pair<String,String>>> cGO=geneOGMap;
                    final HashMap<Integer,Integer> taxTranslationF  = taxTranslation;
                    final int taxidF = taxid;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_num = index;
                        @Override
                        public void run() {
                            System.out.println("Proc num: "+proc_num);
                            System.out.println("Sim size: "+simList.get(proc_num).neighbours.keySet().size());
                            System.out.println("TaxID: "+taxesUsed.get(proc_num));
                            try{
                               contList.get(proc_num).addContingency(simList.get(proc_num),c,cGO,taxesUsed.get(proc_num),taxTranslationF);
                            }
                            catch(Exception e){
                               System.out.println("Something went wrong with proc: "+proc_num+" "+taxesUsed.get(proc_num)); 
                               e.printStackTrace();
                            }
                         
                        }
                          });
                    }
                    }catch(Exception e){
                        e.printStackTrace();
                    }
                finally {
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                taxesUsed.clear();
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
            }
            
            if(batchCount!=0){
                batchSize = batchCount;
                for(int i=0;i<batchSize;i++){
                        Contingency cotTmp=new Contingency();
                       cotTmp.createGOCog(cgmap);
                        cotTmp.initializeContingency(cgmap);
                        contList.add(cotTmp);
                    }
                    
                ExecutorService exec = Executors.newFixedThreadPool(batchCount);
                
                final Lock lock = new ReentrantLock();
                
                try {
                    final OGGOMapping c=cgmap;
                    final HashMap<String,HashSet<Pair<String,String>>> cGO=geneOGMap;
                     final HashMap<Integer,Integer> taxTranslationF  = taxTranslation;
                    final int taxidF = taxid;
                    for (int i=0;i<batchSize;i++) {
                        final int index=i;
                        exec.submit(new Runnable() {
                            int proc_ind = index;
                        @Override
                        public void run() {
                             System.out.println("Proc num: "+proc_ind);
                            System.out.println("Sim size: "+simList.get(proc_ind).neighbours.keySet().size());
                            System.out.println("TaxID: "+taxesUsed.get(proc_ind));
                            try{
                                 contList.get(proc_ind).addContingency(simList.get(proc_ind),c,cGO,taxesUsed.get(proc_ind),taxTranslationF);
                            }
                            catch(Exception e){
                                System.out.println("Something went wrong: "+proc_ind+" "+taxesUsed.get(proc_ind));
                                e.printStackTrace();
                            }
                        }
                          });
                    }
                    }catch(Exception e){
                        e.printStackTrace();
                    }
                finally {
                        exec.shutdown();
                        
                    }
                
                try {
                    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    
                for(int i=0;i<batchSize;i++){
                    for(Integer ind:cot.contTable.keySet())
                        for(Integer ind1:cot.contTable.keySet())
                            for(k=0;k<4;k++){
                        int val=cot.contTable.get(ind).get(ind1).get(k);
                        val=val+contList.get(i).contTable.get(ind).get(ind1).get(k);
                        cot.contTable.get(ind).get(ind1).set(k, val);
                                }
                }
                batchCount=0;
                contList.clear();
                simList.clear();
                taxesUsed.clear();
            }
            System.out.println("Number of locations: "+numLocs);
            cot.writeContingency(cgmap, output);
        } catch (Exception e) {
            e.printStackTrace();
        }
       
   }
        
     void createBaselineOrganisms(String taxIDFilePath,String[] extensions, boolean recursive, int k,OGGOMapping cgmap,HashMap<String,HashSet<String>> geneOGMap, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal."))
                     continue;

                  Gene cg=new Gene();
                  cg.findCogs(input,geneOGMap,map); //change to 0 to get COGs from file
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }

                  GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganism(sim, cgmap,geneOGMap,map);
                }
                else{
                cbf.createFeatures(sim, cgmap,geneOGMap,map);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                cg.genes.clear();
                cg.anotGen.clear();
                sim.neighbours.clear();
            }
            cbf.normalize();
            
         
            System.out.println("Total number of organizms: "+allTax.size());
            
            HashSet<Integer> taxSP=new HashSet();
            
            for(int s:allTax){
                taxSP.add(taxTranslation.get(s));
                System.out.println("taxID: "+taxTranslation.get(s));
            }
            System.out.println("Number of species: "+taxSP.size());
            
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createTrainHeader(output, cgmap, headerInput, isTrain);
            tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
     
     public static double calculateMean(double[] m){
    if (m.length == 0) return 0; // don't divide by zero!

    double mean=0;
    
   for(int i=0;i<m.length;i++)
       mean=mean+m[i];
   
   return mean/m.length;
}

public static double calculateStandardDeviation(double[] sd){

    if (sd.length == 0) return 0; // don't divide by zero!

    double sum = 0;

    double mean = calculateMean(sd);

    for (int i = 0; i<sd.length; i++){
        sum = sum + (sd[i] - mean) * (sd[i] - mean);
    }
    double squaredDiffMean = (sum) / (sd.length); 
    double standardDev = (Math.sqrt(squaredDiffMean)); 

    return standardDev; 
}
     
     //guava
     void countOGOccurence(String taxIDFilePath, File ensembleTIDsFile,String[] extensions, boolean recursive, CreateReducedMappingsFile rmf, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, String output, HashMap<Integer,Integer> taxTranslation,  Mappings map){
         
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
         
         
           try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;         
         
          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
             
             
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
                  
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal."))
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;
                
            if(taxTranslation.containsKey(taxid)){
                if(eggnogTaxs.contains(taxid)){found=1;}
                else if(eggnogTaxs.contains(taxTranslation.get(taxid)) && taxTranslation.get(taxid)==taxid && ensembleTaxs.contains(taxid)) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && eggnogTaxs.contains(taxTranslation.get(taxid))) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && !eggnogTaxs.contains(taxTranslation.get(taxid))){

                    for(int td:eggnogTaxs){
                        if(!taxTranslation.containsKey(td))
                            continue;
                        if(taxTranslation.get(td).equals(taxTranslation.get(taxid))){
                            found=1;
                            break;
                        }
                  }
                    
                }
            }
                System.out.println("found: "+found);
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
                
                  Gene cg=new Gene();
                  cg.findCogsNew(input,geneOGMap,map); //change to 0 to get COGs from file
                  usedTIDs.add(taxid);
                 
                  for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                              continue;
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
                  
                usedSpecies.add(taxTranslation.get(taxid));
            }

         FileWriter fw = new FileWriter(output); //the true will append the new data

        it = ogOccCount.keySet().iterator();
        
        while(it.hasNext()){
            String og= it.next();
            ArrayList<Integer> ocOg = new ArrayList();
            
            HashMap<Integer,Integer> taxOcc = ogOccCount.get(og);
            
            Iterator<Integer> taxit = taxOcc.keySet().iterator();
            
            while(taxit.hasNext()){
                int ti = taxit.next();
                if(taxOcc.get(ti)>0)
                ocOg.add(taxOcc.get(ti));
            }
            
            
            double tmp[] = new double[ocOg.size()];
            
            for(int kt=0;kt<ocOg.size();kt++)
                tmp[kt]=ocOg.get(kt);
            
            double oc = calculateMean(tmp);
            
            int norg = ogOrgCount.get(og).size();
            
            for(int kt=0;kt<ocOg.size();kt++)
                fw.write(tmp[kt]+" ");
            fw.write("\n");
            
            fw.write(og+" "+(oc)+" "+calculateStandardDeviation(tmp)+"\n");
            
        }
        fw.close();
         
        } catch (Exception e) {
            e.printStackTrace();
        }      
     }
     
 void createBaselineOrganismsNewLeaveOut(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String outputTrain, String outputTest,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, HashSet<String> leaveoutOGs, int useNC, int OrgType){
           HashSet<Integer> okTaxes = new HashSet<>();
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
      
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
         
         if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();

         int numGenomes=0, numFiles=0;
         
          HashMap<Integer, HashSet<String>> organismOGs = new HashMap<>();

          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                

                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                System.out.println("found: "+found);
                
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                  
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){   
                  Gene cg=new Gene();
                  if(OrgType==0)
                     cg.findCogsNew(input,geneOGMap,map);//fungy
                  else if(OrgType ==1)
                       cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }

                  GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                
                System.out.println("NN: "+sim.neighbours.keySet().size());
                
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCountLeaveOut(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,leaveoutOGs,taxid);
                }
                else{
                cbf.createFeaturesNewCountLeaveOut(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,leaveoutOGs,taxid);//modify so that genes containing COGs from leaveout cogs are not counted in feature construction
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                 
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                              continue;
                              
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
                 
                cg.genes.clear();
                cg.anotGen.clear();
                sim.neighbours.clear();
                 }
                 else{
                  
                   Gene cg=new Gene();

                   int maxs=-1;
                   if(OrgType==0)
                        maxs=cg.findCogsNewnonCH(input,geneOGMap,map);//fungy
                   if(OrgType==1)
                       maxs=cg.findCogsNewnonCHMetazoa(input,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova .nonchromosomal");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
                  
                  if(maxs<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova u contigu");
                      continue;
                  }
                  
                  if(cg.genesContig.size()>1000 && OrgType==0){
                      System.out.println("To many contigs!");
                      continue;
                  }
                  
                  ArrayList<GeneNeighbours> nSims = new ArrayList<>();
                  
                  for(int ccs = 0; ccs < cg.genesContig.size(); ccs++){
                      cg.genes = cg.genesContig.get(ccs);
                     GeneNeighbours sim1=new GeneNeighbours();
                     sim1.computeSpatialNeighboors(cg, k);
                     nSims.add(sim1);
                  }

                  int exist  = 0;
                  for(int ccs = 0; ccs < nSims.size();ccs++)
                      if(nSims.get(ccs).numCOGs>0)
                          exist=1;
                  
                  if(exist==0)
                      continue;
                  
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCountnonCHLeaveOut(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,leaveoutOGs,taxid);
                    System.out.println("Success .nonchromosomal");
                }
                else{
                cbf.createFeaturesNewCountNonCHLeaveOut(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,leaveoutOGs,taxid);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                 System.out.println("Success .nonchromosomal");
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                
               for(int gg=0;gg<cg.genesContig.size();gg++){ 
                   cg.genes = cg.genesContig.get(gg);
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                               HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                              ogOccCount.put(ogsOCC, taxOcc);
                              
                          }
                              
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
               }
                 
                 
                cg.genes.clear();
                cg.anotGen.clear();
               // sim.neighbours.clear();
                nSims.clear();
                 }
                usedSpecies.add(taxTranslation.get(taxid));
            }
            cbf.normalize();
            
            FileWriter fw = new FileWriter("NOGCountnonCH.txt"); //the true will append the new data

        it = ogOccCount.keySet().iterator();
        
        while(it.hasNext()){
            String og= it.next();
            ArrayList<Integer> ocOg = new ArrayList();
            
            HashMap<Integer,Integer> taxOcc = ogOccCount.get(og);
            
            Iterator<Integer> taxit = taxOcc.keySet().iterator();
            
            while(taxit.hasNext()){
                int ti = taxit.next();
                if(taxOcc.get(ti)>0)
                ocOg.add(taxOcc.get(ti));
            }
            
            
            double tmp[] = new double[ocOg.size()];
            
            for(int kt=0;kt<ocOg.size();kt++)
                tmp[kt]=ocOg.get(kt);
            
            double oc = calculateMean(tmp);
            
            int norg = ogOrgCount.get(og).size();
            
            for(int kt=0;kt<ocOg.size();kt++)
                fw.write(tmp[kt]+" ");
            fw.write("\n");
            
            fw.write(og+" "+(oc)+" "+calculateStandardDeviation(tmp)+"\n");
            
        }
        fw.close();
        
         try{
                fw = new FileWriter("usedTaxFilenonCH.txt"); 
                
               /* Iterator<String>*/ it=taxFNMapping.keySet().iterator();
                 
                 while(it.hasNext()){
                     String fn=it.next();
                     
                     int f=taxFNMapping.get(fn);
                     
                     fw.write(fn+" "+f);
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }
                 
                 
                     try{
                fw = new FileWriter("usedTaxesnonCH.txt"); 
                
                Iterator<Integer> it1=usedTIDs.iterator();
                 
                 while(it1.hasNext()){
                     Integer tid=it1.next();
                     
                     fw.write(tid+"");
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }             
        
         
            System.out.println("Total number of organizms: "+allTax.size());
            
            HashSet<Integer> taxSP=new HashSet();
            
            for(int s:allTax){
                taxSP.add(taxTranslation.get(s));
                System.out.println("taxID: "+taxTranslation.get(s));
            }
            System.out.println("Number of species: "+taxSP.size());
            
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createTrainHeader(outputTrain, cgmap, headerInput, isTrain);
            tt.appendBaselineRowsToSetTrain(outputTrain, cgmap, cbf, leaveoutOGs ,isTrain);//create test data in the same way
             tt.createTrainHeader(outputTest, cgmap, headerInput, isTrain);
            tt.appendBaselineRowsToSetTest(outputTest, cgmap, cbf, leaveoutOGs ,isTrain);//create test data in the same way
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }     
     
 
  void createDistanceForestNewLeaveOut(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k, int useLogDist,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String outputTrain, String outputTest,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, HashSet<String> leaveoutOGs, int useNC, int OrgType){
           HashSet<Integer> okTaxes = new HashSet<>();
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
      
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
         
         if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();

         int numGenomes=0, numFiles=0;
         
          HashMap<Integer, HashSet<String>> organismOGs = new HashMap<>();

          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
              
               LocationSimilarity sim=new LocationSimilarity();
                      sim.initializeSimilarity(cgmap);
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;           
                  
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                System.out.println("found: "+found);
                if(found==0 /*|| (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid))*/)
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  if(OrgType == 0)
                    cg.findCogsNew(input, geneOGMap, map); 
                  if(OrgType==1)
                    cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
              
                  if(useLogDist==0)//distance on circular genome
                    sim.computeSimilarity(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                      sim.computeLogSimilarity(cg, cgmap,geneOGMap, taxTranslation ,taxid);
                                
            }          
            else{
                   Gene cg=new Gene();
                   if(OrgType==0)
                     cg.findCogsNewnonCH(input, geneOGMap, map); 
                   if(OrgType==1)
                      cg.findCogsNewnonCHMetazoa(input, geneOGMap, map,taxid);
                   
                    if(cg.genesContig.size()>1000 && OrgType ==0){
                        System.out.println("To many contigs!");
                        continue;
                    }
                    
                    if(useLogDist==0)//distance on circular genome
                            sim.computeSimilaritynonCH(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                            sim.computeLogSimilaritynonCH(cg, cgmap,geneOGMap, taxTranslation ,taxid);                   
                    }
            usedSpecies.add(taxTranslation.get(taxid));
                    }
            
            Iterator<String> it1 = sim.locationN.keySet().iterator();
            
            System.out.println("Before normalize");
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         sim.normalize(useLogDist);
         
         System.out.println("After normalize");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         
         sim.removeEmpty();
         
         System.out.println("After remove");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createLocationTrainHeaderHoldout(outputTrain, cgmap, sim, leaveoutOGs, headerInput, isTrain);
            tt.appendLocationRowsToSetHoldout(outputTrain, cgmap, sim, leaveoutOGs, isTrain, 0);
            tt.createLocationTrainHeaderHoldout(outputTest, cgmap, sim, leaveoutOGs, headerInput, isTrain);
            tt.appendLocationRowsToSetHoldout(outputTest, cgmap, sim, leaveoutOGs, isTrain, 1);
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        } 
   }     
 
     
      void createBaselineOrganismsNew(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC, int OrgType){
           HashSet<Integer> okTaxes = new HashSet<>();
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
      
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
         
         if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         
          HashMap<Integer, HashSet<String>> organismOGs = new HashMap<>();
     
          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                

                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                System.out.println("found: "+found);
                
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                  
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){   
                  Gene cg=new Gene();
                  if(OrgType==0)
                     cg.findCogsNew(input,geneOGMap,map);//fungy
                  else if(OrgType ==1)
                       cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }

                  GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                
                System.out.println("NN: "+sim.neighbours.keySet().size());
                
                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCount(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                }
                else{
                cbf.createFeaturesNewCount(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                 
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                              continue;                             
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
                 
                cg.genes.clear();
                cg.anotGen.clear();
                sim.neighbours.clear();
                 }
                 else{
                  
                   Gene cg=new Gene();
                   int maxs=-1;
                   if(OrgType==0)
                        maxs=cg.findCogsNewnonCH(input,geneOGMap,map);//fungy
                   if(OrgType==1)
                       maxs=cg.findCogsNewnonCHMetazoa(input,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova .nonchromosomal");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
                  
                  if(maxs<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova u contigu");
                      continue;
                  }
                  
                  if(cg.genesContig.size()>1000 && OrgType==0){
                      System.out.println("To many contigs!");
                      continue;
                  }
                  
                  ArrayList<GeneNeighbours> nSims = new ArrayList<>();
                  
                  for(int ccs = 0; ccs < cg.genesContig.size(); ccs++){
                      cg.genes = cg.genesContig.get(ccs);
                     GeneNeighbours sim1=new GeneNeighbours();
                     sim1.computeSpatialNeighboors(cg, k);
                     nSims.add(sim1);
                  }
				  
                  int exist  = 0;
                  for(int ccs = 0; ccs < nSims.size();ccs++)
                      if(nSims.get(ccs).numCOGs>0)
                          exist=1;
                  
                  if(exist==0)
                      continue;
                  
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCountnonCH(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                    System.out.println("Success .nonchromosomal");
                }
                else{
                cbf.createFeaturesNewCountNonCH(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                 System.out.println("Success .nonchromosomal");
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
               
               for(int gg=0;gg<cg.genesContig.size();gg++){ 
                   cg.genes = cg.genesContig.get(gg);
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                               HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                              ogOccCount.put(ogsOCC, taxOcc);
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
               }
                 
                 
                cg.genes.clear();
                cg.anotGen.clear();
                nSims.clear();
                 }
                usedSpecies.add(taxTranslation.get(taxid));
            }
            cbf.normalize();
            
            FileWriter fw = new FileWriter("NOGCountnonCH.txt"); //the true will append the new data

        it = ogOccCount.keySet().iterator();
        
        while(it.hasNext()){
            String og= it.next();
            ArrayList<Integer> ocOg = new ArrayList();
            
            HashMap<Integer,Integer> taxOcc = ogOccCount.get(og);
            
            Iterator<Integer> taxit = taxOcc.keySet().iterator();
            
            while(taxit.hasNext()){
                int ti = taxit.next();
                if(taxOcc.get(ti)>0)
                ocOg.add(taxOcc.get(ti));
            }
            
            
            double tmp[] = new double[ocOg.size()];
            
            for(int kt=0;kt<ocOg.size();kt++)
                tmp[kt]=ocOg.get(kt);
            
            double oc = calculateMean(tmp);
            
            int norg = ogOrgCount.get(og).size();
            
            for(int kt=0;kt<ocOg.size();kt++)
                fw.write(tmp[kt]+" ");
            fw.write("\n");
            
            fw.write(og+" "+(oc)+" "+calculateStandardDeviation(tmp)+"\n");
            
        }
        fw.close();
        
         try{
                fw = new FileWriter("usedTaxFilenonCH.txt"); 
                
               /* Iterator<String>*/ it=taxFNMapping.keySet().iterator();
                 
                 while(it.hasNext()){
                     String fn=it.next();
                     
                     int f=taxFNMapping.get(fn);
                     
                     fw.write(fn+" "+f);
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }
                 
                 
                     try{
                fw = new FileWriter("usedTaxesnonCH.txt"); 
                
                Iterator<Integer> it1=usedTIDs.iterator();
                 
                 while(it1.hasNext()){
                     Integer tid=it1.next();
                     
                     fw.write(tid+"");
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }             
        
         
            System.out.println("Total number of organizms: "+allTax.size());
            
            HashSet<Integer> taxSP=new HashSet();
            
            for(int s:allTax){
                taxSP.add(taxTranslation.get(s));
                System.out.println("taxID: "+taxTranslation.get(s));
            }
            System.out.println("Number of species: "+taxSP.size());
            
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createTrainHeader(output, cgmap, headerInput, isTrain);
            tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
      
       void createBaselineOrganismsNewNoSintheny(String taxIDFilePath, File ensembleTIDsFile, SynthenicBlocks blocks, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC, int OrgType){
    //OrgType - 0-fungy, 1-metazoa
           HashSet<Integer> okTaxes = new HashSet<>();
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
       //maps file name with the organism taxID
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         //maps taxID with the first organism it occurs
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null){
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
      
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
         
         if(OrgType ==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
         
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         
          HashMap<Integer, HashSet<String>> organismOGs = new HashMap<>();

          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
              
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }

      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                
                
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                System.out.println("found: "+found);
                
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                  
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){   
                  Gene cg=new Gene();
                  if(OrgType==0)
                     cg.findCogsNewNoSyntheni(input,blocks,taxTranslation,taxid,geneOGMap,map);//fungy
                  else if(OrgType ==1)
                       cg.findCogsNewMetazoaNoSyntheni(input,blocks,taxTranslation,taxid,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }

                  GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                
                System.out.println("NN: "+sim.neighbours.keySet().size());

                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCount(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                }
                else{
                cbf.createFeaturesNewCount(sim, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                 
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                              continue;                             
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
                 
                cg.genes.clear();
                cg.anotGen.clear();
                sim.neighbours.clear();
                 }
                 else{
                  
                   Gene cg=new Gene();
                   int maxs=-1;
                   if(OrgType==0)
                        maxs=cg.findCogsNewnonCHNoSyntheni(input,blocks,taxTranslation,taxid,geneOGMap,map);//fungy
                   if(OrgType==1)
                       maxs=cg.findCogsNewnonCHMetazoaNoSyntheni(input,blocks,taxTranslation,taxid,geneOGMap,map,taxid);
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova .nonchromosomal");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
                  
                  if(maxs<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova u contigu");
                      continue;
                  }
                  
                  if(cg.genesContig.size()>1000 && OrgType==0){
                      System.out.println("To many contigs!");
                      continue;
                  }
                  
                  ArrayList<GeneNeighbours> nSims = new ArrayList<>();
                  
                  for(int ccs = 0; ccs < cg.genesContig.size(); ccs++){
                      cg.genes = cg.genesContig.get(ccs);
                     GeneNeighbours sim1=new GeneNeighbours();
                     sim1.computeSpatialNeighboors(cg, k);
                     nSims.add(sim1);
                  }

                  int exist  = 0;
                  for(int ccs = 0; ccs < nSims.size();ccs++)
                      if(nSims.get(ccs).numCOGs>0)
                          exist=1;
                  
                  if(exist==0)
                      continue;
                  
                System.out.println("Neighbours computed!");
                if(taxIDMap.containsKey(input.getName().substring(0,input.getName().length()-4)))
                 allTax.add(taxIDMap.get(input.getName().substring(0,input.getName().length()-4)));
                    
                    OKTIDs.add(taxid);
                if(usedTIDs.contains(taxid)){
                    cbf.createFeaturesOrganismNewCountnonCH(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                    System.out.println("Success .nonchromosomal");
                }
                else{
                cbf.createFeaturesNewCountNonCH(nSims, cgmap,geneOGMap,map,taxTranslation,organismOGs,taxid);
                 usedTIDs.add(taxid);
                 System.out.println("Broj obradenih organizama: "+usedTIDs.size());
                 System.out.println("Success .nonchromosomal");
                }
                
                System.out.println("Features computed!");
                
                numGenomes++;
                System.out.println("Obradeno dokumenata: "+numGenomes);
                System.out.println("Broj OK organizama: "+OKTIDs.size());
                
               for(int gg=0;gg<cg.genesContig.size();gg++){ 
                   cg.genes = cg.genesContig.get(gg);
                 for(int gi=0;gi<cg.genes.size();gi++){
                      
                      if(!geneOGMap.containsKey(cg.genes.get(gi).getValue2()))
                          continue;
                      
                      HashSet<Pair<String,String>> ogs = geneOGMap.get(cg.genes.get(gi).getValue2());
                      
                      HashSet<String> ogsg = new HashSet<>();
                      
                      for(Pair<String,String> p:ogs)
                          if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                                 ogsg.add(p.getValue0());
                      
                     for(String ogsOCC:ogsg){ 
                      if(ogOrgCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                          ogOrgCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/).add(taxTranslation.get(taxid));
                          
                          if(!ogOccCount.containsKey(ogsOCC/*cg.genes.get(gi).getValue2()*/)){
                               HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                              ogOccCount.put(ogsOCC, taxOcc);
                          }
                              
                        
                          HashMap<Integer,Integer> taxOcc = ogOccCount.get(ogsOCC/*cg.genes.get(gi).getValue2()*/);
                          
                          if(!taxOcc.containsKey(taxTranslation.get(taxid)))
                              taxOcc.put(taxTranslation.get(taxid), 0);
                          
                          int occ = taxOcc.get(taxTranslation.get(taxid));
                          occ++;
                          System.out.println("occ: "+occ);
                          taxOcc.put(taxTranslation.get(taxid), occ);
                          ogOccCount.put(ogsOCC/*cg.genes.get(gi).getValue2()*/, taxOcc);
                      }
                    }
                  }
               }
                 
                 
                cg.genes.clear();
                cg.anotGen.clear();
                nSims.clear();
                 }
                usedSpecies.add(taxTranslation.get(taxid));
            }
            cbf.normalize();
            
            FileWriter fw = new FileWriter("NOGCountnonCH.txt"); //the true will append the new data

        it = ogOccCount.keySet().iterator();
        
        while(it.hasNext()){
            String og= it.next();
            ArrayList<Integer> ocOg = new ArrayList();
            
            HashMap<Integer,Integer> taxOcc = ogOccCount.get(og);
            
            Iterator<Integer> taxit = taxOcc.keySet().iterator();
            
            while(taxit.hasNext()){
                int ti = taxit.next();
                if(taxOcc.get(ti)>0)
                ocOg.add(taxOcc.get(ti));
            }
            
            
            double tmp[] = new double[ocOg.size()];
            
            for(int kt=0;kt<ocOg.size();kt++)
                tmp[kt]=ocOg.get(kt);
            
            double oc = calculateMean(tmp);
            
            int norg = ogOrgCount.get(og).size();
            
            for(int kt=0;kt<ocOg.size();kt++)
                fw.write(tmp[kt]+" ");
            fw.write("\n");
            
            fw.write(og+" "+(oc)+" "+calculateStandardDeviation(tmp)+"\n");
            
        }
        fw.close();
        
         try{
                fw = new FileWriter("usedTaxFilenonCH.txt"); 
                
               /* Iterator<String>*/ it=taxFNMapping.keySet().iterator();
                 
                 while(it.hasNext()){
                     String fn=it.next();
                     
                     int f=taxFNMapping.get(fn);
                     
                     fw.write(fn+" "+f);
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }
                 
                 
                     try{
                fw = new FileWriter("usedTaxesnonCH.txt"); 
                
                Iterator<Integer> it1=usedTIDs.iterator();
                 
                 while(it1.hasNext()){
                     Integer tid=it1.next();
                     
                     fw.write(tid+"");
                   
                     fw.write("\n");
                 }
                
                fw.close();
              }
               catch(IOException e){
              e.printStackTrace();
            }             
        
         
            System.out.println("Total number of organizms: "+allTax.size());
            
            HashSet<Integer> taxSP=new HashSet();
            
            for(int s:allTax){
                taxSP.add(taxTranslation.get(s));
                System.out.println("taxID: "+taxTranslation.get(s));
            }
            System.out.println("Number of species: "+taxSP.size());
            
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createTrainHeader(output, cgmap, headerInput, isTrain);
            tt.appendBaselineRowsToSet(output, cgmap, cbf, isTrain);
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
      
      
      public void countGOFrequency(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf, String[] extensions, boolean recursive, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap,Mappings map, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC){
          
          HashMap<String,Integer> goFrequency=new HashMap<>();
          int NumGenesWithOG = 0;
          
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         HashSet<Integer> ensembleTaxs = new HashSet();
         HashSet<Integer> eggnogTaxs = new HashSet();
         BufferedReader reader1;
         try {
             Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
             reader1 = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
             while ((line = reader1.readLine()) != null) {
                 
                 ensembleTaxs.add(Integer.parseInt(line.trim()));
             }
             reader1.close();
         }
         catch(Exception e){e.printStackTrace();}
         Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
         while(it.hasNext()){
             
             HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
             
             if(ns.size()==0)
                 continue;
             
             for(Pair<String,String> nss:ns){
                 eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
             }
         }
         HashSet<Integer> usedSpecies = new HashSet();
         System.out.println("eggnogTaxs: "+eggnogTaxs.size());
         System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         for (Iterator iterator = files.iterator(); iterator.hasNext();){//read the files
             numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
               if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                

            if(taxTranslation.containsKey(taxid)){
                if(eggnogTaxs.contains(taxid)){found=1;}
                else if(eggnogTaxs.contains(taxTranslation.get(taxid)) && taxTranslation.get(taxid)==taxid && ensembleTaxs.contains(taxid)) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && eggnogTaxs.contains(taxTranslation.get(taxid))) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && !eggnogTaxs.contains(taxTranslation.get(taxid))){

                    for(int td:eggnogTaxs){
                        if(!taxTranslation.containsKey(td))
                            continue;
                        if(taxTranslation.get(td).equals(taxTranslation.get(taxid))){
                            found=1;
                            break;
                        }
                  }
                    
                }
            }
            else continue;
            
             if(taxTranslation.get(taxid) == 4932 && taxid!=764097)
                          found=0;
                      else if(taxid==764097) found=1;
            
                System.out.println("found: "+found);
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
            
                 taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  cg.findCogsNew(input,geneOGMap,map); //change to 0 to get COGs from file
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                                   //the code goes here
                 //geneogmap
                 //cgmap
                for(int i=0;i<cg.genes.size();i++){
                    if(!geneOGMap.containsKey(cg.genes.get(i).getValue2()))
                        continue;
                    HashSet<Pair<String,String>> tmp = geneOGMap.get(cg.genes.get(i).getValue2());
                     HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:tmp){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
            HashSet<String> gos = new HashSet<>();
            
            for(String og:OGs){
                ArrayList<String> t = cgmap.CogGOmap.get(og);
                if(t==null)
                    continue;
                gos.addAll(t);
            }
            
            for(String go:gos){
                if(!goFrequency.containsKey(go)){
                    goFrequency.put(go, 1);
                }
                else{
                    int f = goFrequency.get(go);
                    f=f+1;
                    goFrequency.put(go, f);
                }
            }
            NumGenesWithOG++;
                }
                  
         }       
       }
         
         Iterator<String> itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             System.out.print(go+" "+freq+"\n");
         }  
         
         try{
           FileWriter fw = new FileWriter(new File("GOFrequency.txt")); 
          
           itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             fw.write(go+" "+freq+"\n");
         }  
         fw.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }
         
      }
      
        public void countGOFrequencyNew(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf, String[] extensions, boolean recursive, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap,Mappings map, boolean randomize, HashMap<Integer,Integer> taxTranslation, int k, int useNC){
          
          HashSet<Integer> okTaxes = new HashSet<>();
                  
                  
         okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
            
            
            
          HashMap<String,Double> goFrequency=new HashMap<>();
          double NumGenesWithOG = 0.0;
          
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         HashSet<Integer> ensembleTaxs = new HashSet();
         HashSet<Integer> eggnogTaxs = new HashSet();
         BufferedReader reader1;
         try {
             Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
             reader1 = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
             while ((line = reader1.readLine()) != null) {
                 
                 ensembleTaxs.add(Integer.parseInt(line.trim()));
             }
             reader1.close();
         }
         catch(Exception e){e.printStackTrace();}
         Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
         while(it.hasNext()){
             
             HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
             
             if(ns.size()==0)
                 continue;
             
             for(Pair<String,String> nss:ns){
                 eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
             }
         }
         HashSet<Integer> usedSpecies = new HashSet();
         System.out.println("eggnogTaxs: "+eggnogTaxs.size());
         System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         int numLoc = 0;
         for (Iterator iterator = files.iterator(); iterator.hasNext();){//read the files
             numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
               if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                
                
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                if(found==0)
                    continue;

                 taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){
                      System.out.println("found: "+found);
                  Gene cg=new Gene();
                  cg.findCogsNew(input,geneOGMap,map); //change to 0 to get COGs from file
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                   GeneNeighbours sim=new GeneNeighbours();
                sim.computeSpatialNeighboors(cg,k);
                
                System.out.println("NN: "+sim.neighbours.keySet().size());

                if(sim.numCOGs==0){
                    System.out.println("Nema COG-ova s funkcijom u susjedstvu!");
                    continue;
                }
                                  
                 Iterator<Triplet<Integer,Integer,String>> keySetIterator = sim.neighbours.keySet().iterator();
               numLoc+=sim.neighbours.keySet().size();
               System.out.println("Sim size: "+sim.neighbours.keySet().size());
                 while(keySetIterator.hasNext()){
            Triplet<Integer,Integer,String> key = keySetIterator.next();
                
                    if(!geneOGMap.containsKey(key.getValue2()))
                        continue;
                    HashSet<Pair<String,String>> tmp = geneOGMap.get(key.getValue2());
                     HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:tmp){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
            HashSet<String> gos = new HashSet<>();
            
            for(String og:OGs){
                ArrayList<String> t = cgmap.CogGOmap.get(og);
                if(t==null)
                    continue;
                gos.addAll(t);
            }
            
            for(String go:gos){
                if(!goFrequency.containsKey(go)){
                    goFrequency.put(go, 1.0);
                }
                else{
                    double f = goFrequency.get(go);
                    f=f+1.0;
                    goFrequency.put(go, f);
                }
            }
            NumGenesWithOG++;
                }
                  
         }       
       }
         
         Iterator<String> itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             System.out.print(go+" "+freq+"\n");
         }  
         
         System.out.println("NumGenesWithOG: "+NumGenesWithOG);
         System.out.println("Num locations1"+numLoc);
         
         try{
           FileWriter fw = new FileWriter(new File("GOFrequencyConsistent.txt")); 
          
           itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             fw.write(go+" "+freq+" "+goFrequency.get(go)+"\n");
         }  
         fw.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }       
      }    
      
       public void countGOFrequencyMetazoa(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf, String[] extensions, boolean recursive, OGGOMapping cgmap, HashMap<String,HashSet<Pair<String,String>>> geneOGMap,Mappings map, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC){
          
          HashMap<String,Integer> goFrequency=new HashMap<>();
          int NumGenesWithOG = 0;
          
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;
         HashSet<Integer> ensembleTaxs = new HashSet();
         HashSet<Integer> eggnogTaxs = new HashSet();
         BufferedReader reader1;
         try {
             Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
             reader1 = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
             while ((line = reader1.readLine()) != null) {
                 
                 ensembleTaxs.add(Integer.parseInt(line.trim()));
             }
             reader1.close();
         }
         catch(Exception e){e.printStackTrace();}
         Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
         while(it.hasNext()){
             
             HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
             
             if(ns.size()==0)
                 continue;
             
             for(Pair<String,String> nss:ns){
                 eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
             }
         }
         HashSet<Integer> usedSpecies = new HashSet();
         System.out.println("eggnogTaxs: "+eggnogTaxs.size());
         System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         for (Iterator iterator = files.iterator(); iterator.hasNext();){//read the files
             numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
               if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;                

            if(taxTranslation.containsKey(taxid)){
                if(eggnogTaxs.contains(taxid)){found=1;}
                else if(eggnogTaxs.contains(taxTranslation.get(taxid)) && taxTranslation.get(taxid)==taxid && ensembleTaxs.contains(taxid)) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && eggnogTaxs.contains(taxTranslation.get(taxid))) found=1;
                else if(!ensembleTaxs.contains(taxTranslation.get(taxid)) && !eggnogTaxs.contains(taxTranslation.get(taxid))){

                    for(int td:eggnogTaxs){
                        if(!taxTranslation.containsKey(td))
                            continue;
                        if(taxTranslation.get(td).equals(taxTranslation.get(taxid))){
                            found=1;
                            break;
                        }
                  }
                    
                }
            }
            else continue;
            
             if(taxTranslation.get(taxid) == 4932 && taxid!=764097)
                          found=0;
                      else if(taxid==764097) found=1;
            
                System.out.println("found: "+found);
                if(found==0 || (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid)))
                    continue;
            
                 taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  cg.findCogsNewMetazoa(input,geneOGMap,map,taxid); //change to 0 to get COGs from file
                  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                for(int i=0;i<cg.genes.size();i++){
                    if(!geneOGMap.containsKey(cg.genes.get(i).getValue2()))
                        continue;
                    HashSet<Pair<String,String>> tmp = geneOGMap.get(cg.genes.get(i).getValue2());
                     HashSet<String> OGs = new HashSet<>();
            
            for(Pair<String,String> p:tmp){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    OGs.add(p.getValue0());
            }
            
            HashSet<String> gos = new HashSet<>();
            
            for(String og:OGs){
                ArrayList<String> t = cgmap.CogGOmap.get(og);
                if(t==null)
                    continue;
                gos.addAll(t);
            }
            
            for(String go:gos){
                if(!goFrequency.containsKey(go)){
                    goFrequency.put(go, 1);
                }
                else{
                    int f = goFrequency.get(go);
                    f=f+1;
                    goFrequency.put(go, f);
                }
            }
            NumGenesWithOG++;
                }
                  
         }       
       }
         
         Iterator<String> itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             System.out.print(go+" "+freq+"\n");
         }  
         
         try{
           FileWriter fw = new FileWriter(new File("GOFrequency.txt")); 
          
           itS = goFrequency.keySet().iterator();
         
         while(itS.hasNext()){
             String go = itS.next();
             double freq = goFrequency.get(go);
             freq/=NumGenesWithOG;
             fw.write(go+" "+freq+"\n");
         }  
         fw.close();
         }
         catch(IOException e){
             e.printStackTrace();
         }
         
      }
       
       
       
        void computeAccrition(String taxIDFilePath,String accretionFilePath,String predictionFilePath,String[] extensions, boolean recursive, int k, int organismType, Mappings map, OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation , String output,int isTrain, boolean randomize){
             
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
         
         HashMap<String,Double> accretions = new HashMap<>();
          HashMap<String,HashSet<String>> predictions = new HashMap<>();
         
            try {
                File accretionFile = new File(accretionFilePath);
             Path path =Paths.get(accretionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(",");
          String og = tok[0].trim();
          int num = Integer.parseInt(og);
            
         int  nDigGo=0;
         
               while(num>0){
                  nDigGo++;
                  num/=10;
             }
         
         String GOName = "GO";
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=og;
          
          Double accr = Double.parseDouble(tok[1].trim());
          
          accretions.put(GOName, accr);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
            
            HashSet<Integer> okTaxes = new HashSet<>();
            
             if(organismType==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(organismType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
            
            
       try {
                File predictionFile = new File(predictionFilePath);
             Path path =Paths.get(predictionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null){
          String[] tok=line.split("\t");
          String og = tok[0].trim();

          for(int j=1;j<tok.length;j++){
              if(!predictions.containsKey(og))
                  predictions.put(og, new HashSet<String>());
              
              predictions.get(og).add(tok[j].trim());
          }
      }
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}     
       
       System.out.println("predictions size: "+predictions.keySet().size());
       System.out.println("accretion size: "+accretions.keySet().size());
         
       double onlyKnownAnnotationAccretion=0.0, newlyPredictedAnnotationAccretion=0.0, knownAndPredictedAccretion = 0.0;
       double numGenes = 0.0; 
       
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         Random rand = new Random();
         
         int numGenomes=0, numFiles=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                  String fileName=input.getName().substring(0,input.getName().length()-4);
                   int taxid=taxIDMap.get(fileName);
                   
                   if(!okTaxes.contains(taxid))
                       continue;
                 
                  Gene cg=new Gene();
                  if(organismType == 0)
                      cg.findCogsNewNoDuplicates(input, geneOGMap, map);
                  else if (organismType == 1)
                      cg.findCogsNewMetazoaNoDuplicates(input, geneOGMap, map, taxid);
                  
                  if(randomize==true)
                      cg.randomize(rand);
                                    
                  for(String gen:cg.anotGen){
                       HashSet<Pair<String,String>> tmp = geneOGMap.get(gen);
                     HashSet<String> ogs = new HashSet<>();
            
            for(Pair<String,String> p:tmp){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    ogs.add(p.getValue0());
            }
                    
                      HashSet<String> allgos = new HashSet<>();
                      HashSet<String> predictedgos = new HashSet<>();
                      
                      for(String og:ogs){
                         ArrayList<String> go = cgmap.CogGOmap.get(og);
                         if(!cgmap.CogGOmap.containsKey(og))
                             continue;
                     
                         allgos.addAll(go);
                         if(predictions.containsKey(og))
                             predictedgos.addAll(predictions.get(og));
                      }
                      
                       numGenes++;//add outside loop, proteins should be merged with corresponding genes (not to be counted multiple times)
                      
                      for(String go:allgos){
                          if(!accretions.containsKey(go)){
                              System.out.println("Error allgos: "+go);
                              continue;
                          }
                          if(!predictedgos.contains(go)){
                              onlyKnownAnnotationAccretion+=accretions.get(go);
                          }
                          else if(predictedgos.contains(go)){
                              knownAndPredictedAccretion+=accretions.get(go);
                          }
                      }
                      
                      for(String go:predictedgos){
                           if(!accretions.containsKey(go)){
                              System.out.println("Error predictions: "+go);
                              continue;
                          }
                          if(!allgos.contains(go)){
                              newlyPredictedAnnotationAccretion+=accretions.get(go);
                          }
                      }
                      
                  }
            }
            
            
            FileWriter fw = new FileWriter(output);
            
            fw.write("Only known: "+(onlyKnownAnnotationAccretion/numGenes)+"\n");
            fw.write("Known and predicted: "+(knownAndPredictedAccretion/numGenes)+"\n");
            fw.write("Only predicted: "+(newlyPredictedAnnotationAccretion/numGenes)+"\n");
            fw.close();
            
            System.out.println("Only known: "+(onlyKnownAnnotationAccretion/numGenes));
            System.out.println("Known and predicted: "+(knownAndPredictedAccretion/numGenes));
            System.out.println("Only predicted: "+(newlyPredictedAnnotationAccretion/numGenes));

        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
    
     void computeAccritionDiffOnt(String taxIDFilePath,String accretionFilePath,String predictionFilePath,String[] extensions, boolean recursive, ArrayList<HashSet<String>> categories ,int k, int organismType, Mappings map, OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<Integer,Integer> taxTranslation , String output,int isTrain, boolean randomize){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
          
         
         HashMap<String,Double> accretions = new HashMap<>();
          HashMap<String,HashSet<String>> predictions = new HashMap<>();
         
            try {
                File accretionFile = new File(accretionFilePath);
             Path path =Paths.get(accretionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(",");
          String og = tok[0].trim();
          //fix og
          int num = Integer.parseInt(og);
            
         int  nDigGo=0;
         
               while(num>0){
                  nDigGo++;
                  num/=10;
             }
         
         String GOName = "GO";
              
              while(nDigGo<7){
                  GOName+="0";
                  nDigGo++;
              }
              
              GOName+=og;
          
          Double accr = Double.parseDouble(tok[1].trim());
          
          accretions.put(GOName, accr);
      }
   
      reader.close();
      
         }
         catch(Exception e){e.printStackTrace();}
            
            HashSet<Integer> okTaxes = new HashSet<>();
            
             if(organismType==0){
          okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);//fungy
         }
         if(organismType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);
         }
            
            
       try {
                File predictionFile = new File(predictionFilePath);
             Path path =Paths.get(predictionFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
             
      while ((line = reader.readLine()) != null){
          String[] tok=line.split("\t");
          String og = tok[0].trim();

          for(int j=1;j<tok.length;j++){
              if(!predictions.containsKey(og))
                  predictions.put(og, new HashSet<String>());
              
              predictions.get(og).add(tok[j].trim());
          }
      }
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}     
       
       System.out.println("predictions size: "+predictions.keySet().size());
       System.out.println("accretion size: "+accretions.keySet().size());
         
       double onlyKnownAnnotationAccretionBP=0.0, onlyKnownAnnotationAccretionMF=0.0, onlyKnownAnnotationAccretionCC=0.0, 
              newlyPredictedAnnotationAccretionBP=0.0,newlyPredictedAnnotationAccretionMF=0.0, newlyPredictedAnnotationAccretionCC=0.0,
              knownAndPredictedAccretionBP = 0.0, knownAndPredictedAccretionMF = 0.0, knownAndPredictedAccretionCC = 0.0;
       double numGenes = 0.0; 
       
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         Random rand = new Random();
         
         int numGenomes=0, numFiles=0;
            for (Iterator iterator = files.iterator(); iterator.hasNext();){
                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());

                 
                  String fileName=input.getName().substring(0,input.getName().length()-4);
                   int taxid=taxIDMap.get(fileName);
                   
                   if(!okTaxes.contains(taxid))
                       continue;
                 
                  Gene cg=new Gene();
                  if(organismType == 0)
                      cg.findCogsNewNoDuplicates(input, geneOGMap, map);
                  else if (organismType == 1)
                      cg.findCogsNewMetazoaNoDuplicates(input, geneOGMap, map, taxid);
                  
                  if(randomize==true)
                      cg.randomize(rand);
                                    
                  for(String gen:cg.anotGen){
                       HashSet<Pair<String,String>> tmp = geneOGMap.get(gen);
                     HashSet<String> ogs = new HashSet<>();
            
            for(Pair<String,String> p:tmp){
                if(Integer.parseInt(p.getValue1())==taxid || taxTranslation.get(Integer.parseInt(p.getValue1())).equals(taxTranslation.get(taxid)))
                    ogs.add(p.getValue0());
            }
              
                      ArrayList<HashSet<String>> allgos = new ArrayList<>();
                      ArrayList<HashSet<String>> predictedgos = new ArrayList<>();
                      
                      for(int z=0;z<categories.size();z++){
                          allgos.add(new HashSet<String>());
                          predictedgos.add(new HashSet<String>());
                      }
                      
                      for(String og:ogs){
                         ArrayList<String> go = cgmap.CogGOmap.get(og);
                         if(!cgmap.CogGOmap.containsKey(og))
                             continue;
                     
                         for(String g:go){
                             if(categories.get(0).contains(g))
                                     allgos.get(0).add(g);
                             else if(categories.get(1).contains(g))
                                     allgos.get(1).add(g);
                             else if(categories.get(2).contains(g))
                                     allgos.get(2).add(g);
                         }
                         
                         if(predictions.containsKey(og)){
                              HashSet<String> pr = predictions.get(og);
                              for(String g:pr){
                              if(categories.get(0).contains(g))
                                     predictedgos.get(0).add(g);
                             else if(categories.get(1).contains(g))
                                     predictedgos.get(1).add(g);
                             else if(categories.get(2).contains(g))
                                     predictedgos.get(2).add(g);
                              }
                             
                         }
                         
                      }
                      
                       numGenes++;//add outside loop, proteins should be merged with corresponding genes (not to be counted multiple times)
                    
                    for(int z=0;z<3;z++){   
                      for(String go:allgos.get(z)){
                          if(!accretions.containsKey(go)){
                              System.out.println("Error allgos: "+go);
                              continue;
                          }
                          if(!predictedgos.get(z).contains(go)){
                              if(z==0)
                              onlyKnownAnnotationAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   onlyKnownAnnotationAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   onlyKnownAnnotationAccretionCC+=accretions.get(go);
                          }
                          else if(predictedgos.get(z).contains(go)){
                               if(z==0)
                              knownAndPredictedAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   knownAndPredictedAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   knownAndPredictedAccretionCC+=accretions.get(go);
                          }
                      }
                    
                      
                      for(String go:predictedgos.get(z)){
                           if(!accretions.containsKey(go)){
                              System.out.println("Error predictions: "+go);
                              continue;
                          }
                          if(!allgos.get(z).contains(go)){
                               if(z==0)
                              newlyPredictedAnnotationAccretionBP+=accretions.get(go);
                              else if(z==1)
                                   newlyPredictedAnnotationAccretionMF+=accretions.get(go);
                              else if(z==2)
                                   newlyPredictedAnnotationAccretionCC+=accretions.get(go);
                          }
                      }
                    }
                  }
            }
          
            String outputTmp = output;
            
          for(int i=0;i<3;i++){  
              String o[] = outputTmp.split("\\.");
              String t = "";
              
              if(i==0)
                  t="BP";
              else if(i==1)
                  t="MF";
              else if(i==2)
                  t="CC";
              
              output = o[0]+t+"."+o[1];
            FileWriter fw = new FileWriter(output);
            
            double knownAA = 0.0, knownPredAA = 0.0, newlyPredAA = 0.0;
            
            if(i==0){
                  knownAA = onlyKnownAnnotationAccretionBP;
                  knownPredAA = knownAndPredictedAccretionBP;
                  newlyPredAA = newlyPredictedAnnotationAccretionBP;
            }
              else if(i==1){
                   knownAA = onlyKnownAnnotationAccretionMF;
                   knownPredAA = knownAndPredictedAccretionMF;
                   newlyPredAA = newlyPredictedAnnotationAccretionMF;
              }
              else if(i==2){
                   knownAA = onlyKnownAnnotationAccretionCC;
                   knownPredAA = knownAndPredictedAccretionCC;
                   newlyPredAA = newlyPredictedAnnotationAccretionCC;
              }
            
            fw.write("Only known: "+(knownAA/numGenes)+"\n");
            fw.write("Known and predicted: "+(knownPredAA/numGenes)+"\n");
            fw.write("Only predicted: "+(newlyPredAA/numGenes)+"\n");
            fw.close();
            
            System.out.println("Only known: "+(knownAA/numGenes));
            System.out.println("Known and predicted: "+(knownPredAA/numGenes));
            System.out.println("Only predicted: "+(newlyPredAA/numGenes));
          }

        } catch (Exception e) {
            e.printStackTrace();
        }   
   }
          
        void createAssociationGraph( String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k, int useLogDist,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC, int OrgType){
      
             BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         
         HashSet<Integer> okTaxes = new HashSet<>();

         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);}
         //////////////////////////////////////////////////////////////////////////////////////
         
         ///////////////////////////////////////////fungy
         if(OrgType==0){
         okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);}
         
          try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;         
         
          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
             
             
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
                
               LocationSimilarity sim=new LocationSimilarity();
                      sim.initializeSimilarity(cgmap);
              
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;           
   
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
            
                System.out.println("found: "+found);
                if(found==0 /*|| (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid))*/)
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  if(OrgType == 0)
                    cg.findCogsNew(input, geneOGMap, map); 
                  if(OrgType==1)
                    cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
              
                  if(useLogDist==0)//distance on circular genome
                    sim.computeSimilarity(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                      sim.computeLogSimilarity(cg, cgmap,geneOGMap, taxTranslation ,taxid);
                                
            }          
            else{
                   Gene cg=new Gene();
                   if(OrgType==0)
                     cg.findCogsNewnonCH(input, geneOGMap, map); 
                   if(OrgType==1)
                      cg.findCogsNewnonCHMetazoa(input, geneOGMap, map,taxid);
                   
                    if(cg.genesContig.size()>1000 && OrgType ==0){
                        System.out.println("To many contigs!");
                        continue;
                    }
                    
                    if(useLogDist==0)//distance on circular genome
                            sim.computeSimilaritynonCH(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                            sim.computeLogSimilaritynonCH(cg, cgmap,geneOGMap, taxTranslation ,taxid);                   
                    }
            usedSpecies.add(taxTranslation.get(taxid));
                    }

            Iterator<String> it1 = sim.locationN.keySet().iterator();
            
            System.out.println("Before normalize");
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
             sim.normalize1(useLogDist);
            sim.transformToWeights();
         
         System.out.println("After normalize");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
                  
         System.out.println("After remove");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());            
            tt.createAssociationGraphHeader(output, cgmap, sim, isTrain);
            tt.appendWeightsToGraph(output, cgmap, sim, isTrain);

            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
 }
              
         ArrayList<Double> numGenes( String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k, int useLogDist,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC, int OrgType){
      
             ArrayList<Double> res=new ArrayList<>();
             res.add(0.0); res.add(0.0);
             
             BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         
         HashSet<Integer> okTaxes = new HashSet<>();

         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);}
         //////////////////////////////////////////////////////////////////////////////////////
         
         ///////////////////////////////////////////fungy
         if(OrgType==0){
         okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);}
         
          try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;         
         
          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
             
             
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
       
               LocationSimilarity sim=new LocationSimilarity();
                      sim.initializeSimilarity(cgmap);
              
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;           
   
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
                       
                System.out.println("found: "+found);
                if(found==0 /*|| (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid))*/)
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  if(OrgType == 0)
                    cg.findCogsNew(input, geneOGMap, map); 
                  if(OrgType==1)
                    cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
                  
                  
                  res.set(0, res.get(0)+cg.genes.size());
                   HashSet<String> gos = new HashSet<>();
                  for(Triplet<Integer,Integer,String> t: cg.genes){
                      if(geneOGMap.containsKey(t.getValue2())){
                          HashSet<Pair<String,String>> ogs = geneOGMap.get(t.getValue2());
                          for(Pair<String,String> p:ogs){
                              if(Integer.parseInt(p.getValue1())==taxid){
                                  if(!cgmap.CogGOmap.containsKey(p.getValue0()))
                                      continue;
                                    gos.addAll(cgmap.CogGOmap.get(p.getValue0()));
                              }
                          }
                            res.set(1, res.get(1)+gos.size());
                          gos.clear();
                      }
                  }
                  
                                
            }          
            else{
                   Gene cg=new Gene();
                   if(OrgType==0)
                     cg.findCogsNewnonCH(input, geneOGMap, map); 
                   if(OrgType==1)
                      cg.findCogsNewnonCHMetazoa(input, geneOGMap, map,taxid);
                   
                    if(cg.genesContig.size()>1000 && OrgType ==0){
                        System.out.println("To many contigs!");
                        continue;
                    }
                    
                    HashSet<String> gos = new HashSet<>();
                    
                    for(int z=0;z<cg.genesContig.size();z++){
                         res.set(0, res.get(0)+cg.genesContig.get(z).size());

                     for(Triplet<Integer,Integer,String> t:cg.genesContig.get(z)){
                      if(geneOGMap.containsKey(t.getValue2())){
                          HashSet<Pair<String,String>> ogs = geneOGMap.get(t.getValue2());
                          for(Pair<String,String> p:ogs){
                              if(Integer.parseInt(p.getValue1())==taxid){
                                   if(!cgmap.CogGOmap.containsKey(p.getValue0()))
                                      continue;
                                   gos.addAll(cgmap.CogGOmap.get(p.getValue0()));
                              }
                          }
                          res.set(1, res.get(1)+gos.size());
                          gos.clear();
                      }
                  }
               }

            }
          }
        } catch (Exception e) {
            e.printStackTrace();
        }
       return res;
    }
        
        
            
        void createDistanceForest(String taxIDFilePath, File ensembleTIDsFile, CreateReducedMappingsFile rmf,String[] extensions, boolean recursive, int k, int useLogDist,OGGOMapping cgmap,HashMap<String,HashSet<Pair<String,String>>> geneOGMap, HashMap<String,HashSet<Integer>> ogOrgCount, HashMap<String,HashMap<Integer,Integer>> ogOccCount, Mappings map, File headerInput, String output,int isTrain, boolean randomize, HashMap<Integer,Integer> taxTranslation, int useNC, int OrgType){
    
          BufferedReader reader;
       File taxIDFile=new File(taxIDFilePath);
         HashMap<String,Integer> taxIDMap=new HashMap<>();
         HashSet<Integer> allTax=new HashSet();
         HashMap<String,Integer> taxFNMapping = new HashMap<>();
         HashMap<Integer,Integer> tidOccurence=new HashMap<>();
         
         HashSet<Integer> okTaxes = new HashSet<>();

         if(OrgType==1){
         okTaxes.add(30608); okTaxes.add(30611);  okTaxes.add(8128); okTaxes.add(6334); //metazoa only
         okTaxes.add(7897); okTaxes.add(13616);  okTaxes.add(46245); okTaxes.add(9601);
         okTaxes.add(7070); okTaxes.add(7668);  okTaxes.add(9361); okTaxes.add(9483);
         okTaxes.add(7757); okTaxes.add(9785);  okTaxes.add(6183); okTaxes.add(13735);
         okTaxes.add(9606); okTaxes.add(59729);  okTaxes.add(9371); okTaxes.add(10020);
         okTaxes.add(9913); okTaxes.add(59463);  okTaxes.add(9615); okTaxes.add(7260);
         okTaxes.add(9258); okTaxes.add(6945);  okTaxes.add(13037); okTaxes.add(9739);
         okTaxes.add(10141); okTaxes.add(7245);  okTaxes.add(7244); okTaxes.add(10116);
         okTaxes.add(9031); okTaxes.add(8049);  okTaxes.add(51511); okTaxes.add(9646);
         okTaxes.add(7091); okTaxes.add(9986);  okTaxes.add(121224); okTaxes.add(7227);
         okTaxes.add(9685); okTaxes.add(10228);  okTaxes.add(34740); okTaxes.add(281687);
         okTaxes.add(7217); okTaxes.add(69293);  okTaxes.add(7222); okTaxes.add(10090);
         okTaxes.add(9669); okTaxes.add(31234);  okTaxes.add(8083); okTaxes.add(9544);
         okTaxes.add(9305); okTaxes.add(37347);  okTaxes.add(45351); okTaxes.add(99883);
         okTaxes.add(8090); okTaxes.add(7955);  okTaxes.add(31033); okTaxes.add(9598);
         okTaxes.add(9595); okTaxes.add(8364);  okTaxes.add(61853); okTaxes.add(6238);
         okTaxes.add(7425); okTaxes.add(6239);  okTaxes.add(9796); okTaxes.add(28377);
         okTaxes.add(7176); okTaxes.add(132908);  okTaxes.add(9103); okTaxes.add(54126);
         okTaxes.add(7165); okTaxes.add(9823);  okTaxes.add(7029); okTaxes.add(7719);
         okTaxes.add(43179); okTaxes.add(135651);  okTaxes.add(7159); okTaxes.add(9813);}
         //////////////////////////////////////////////////////////////////////////////////////
         
         ///////////////////////////////////////////fungy
         if(OrgType==0){
         okTaxes.add(578458); okTaxes.add(222929); okTaxes.add(426418); okTaxes.add(4950);
         okTaxes.add(663331); okTaxes.add(644223); okTaxes.add(322104); okTaxes.add(554155);
         okTaxes.add(431241); okTaxes.add(665079); okTaxes.add(426428); okTaxes.add(341663);
         okTaxes.add(578455); okTaxes.add(4956); okTaxes.add(367110); okTaxes.add(214684);
         okTaxes.add(294746); okTaxes.add(294747); okTaxes.add(1064592); okTaxes.add(559305);
         okTaxes.add(402676); okTaxes.add(336963); okTaxes.add(246410); okTaxes.add(573826);
         okTaxes.add(227321); okTaxes.add(379508); okTaxes.add(5061); okTaxes.add(1071381);
         okTaxes.add(242507); okTaxes.add(559295); okTaxes.add(1071378); okTaxes.add(306902);
         okTaxes.add(306901); okTaxes.add(764097); okTaxes.add(876142); okTaxes.add(931890);
         okTaxes.add(413071); okTaxes.add(535722); okTaxes.add(907965); okTaxes.add(526221);
         okTaxes.add(367775); okTaxes.add(284811); okTaxes.add(284812); okTaxes.add(573729);
         okTaxes.add(510516); okTaxes.add(660122); okTaxes.add(985895); okTaxes.add(334819);
         okTaxes.add(240176);}
   
         try {
             Path path =Paths.get(taxIDFile.getAbsolutePath());
             reader = Files.newBufferedReader(path,ENCODING);
             String line = null;
      while ((line = reader.readLine()) != null) {
          String[] tok=line.split(" ");
          String name=tok[0].substring(0,tok[0].length()-4);
          int tid=Integer.parseInt(tok[1]);
          taxIDMap.put(name, tid);
      }
   
      reader.close();
         }
         catch(Exception e){e.printStackTrace();}
          
       try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         DatasetWriter tt=new DatasetWriter();
         OGFeatures cbf=new OGFeatures();
         HashSet<Integer> usedTIDs=new HashSet<>();
         Random rand = new Random();
         HashSet<Integer> OKTIDs=new HashSet<>();
         int numGenomes=0, numFiles=0;         
         
          HashSet<Integer> ensembleTaxs = new HashSet();
              HashSet<Integer> eggnogTaxs = new HashSet();
              
             BufferedReader reader1;
             
             
              try {
                     Path path =Paths.get(ensembleTIDsFile.getAbsolutePath());
                     reader1 = Files.newBufferedReader(path,ENCODING);
                      String line = null;
                     
                      while ((line = reader1.readLine()) != null) {
                            
                           ensembleTaxs.add(Integer.parseInt(line.trim()));
                    }
      reader1.close();
       }
         catch(Exception e){e.printStackTrace();}
              
              Iterator<String> it= rmf.ProteinNOG.keySet().iterator();
              
              while(it.hasNext()){
                  
                  HashSet<Pair<String,String>> ns = rmf.ProteinNOG.get(it.next());
                  
                  if(ns.size()==0)
                     continue;

                 for(Pair<String,String> nss:ns){
                     eggnogTaxs.add(Integer.parseInt(nss.getValue1()));
                 }
              }
              
              HashSet<Integer> usedSpecies = new HashSet();
         
              System.out.println("eggnogTaxs: "+eggnogTaxs.size());
              System.out.println("ensembleTaxs: "+ensembleTaxs.size());
         
              
               LocationSimilarity sim=new LocationSimilarity();
                      sim.initializeSimilarity(cgmap);
              
            for (Iterator iterator = files.iterator(); iterator.hasNext();){

                 numFiles++;
                 System.out.println("Processing file: "+numFiles);
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                 
                 if(input.getAbsolutePath().contains(".nonchromosomal.") && useNC==0)
                     continue;

                 String fileName=input.getName().substring(0,input.getName().length()-4);
                int taxid=taxIDMap.get(fileName);
                int found=0;           
                
                if(okTaxes.contains(taxid)) found=1; //metazoa  + fungy
                    else continue;
                       
                System.out.println("found: "+found);
                if(found==0 /*|| (usedSpecies.contains(taxTranslation.get(taxid)) && !usedTIDs.contains(taxid))*/)
                    continue;
                         
                    taxFNMapping.put(fileName, taxid);
                 
                 if(!input.getAbsolutePath().contains(".nonchromosomal.")){      
                  Gene cg=new Gene();
                  if(OrgType == 0)
                    cg.findCogsNew(input, geneOGMap, map); 
                  if(OrgType==1)
                    cg.findCogsNewMetazoa(input,geneOGMap,map,taxid);  
                  System.out.println("gene and anotated gene size!");
                  System.out.println(cg.genes.size());
                  System.out.println(cg.anotGen.size());
                  
                  if(randomize==true)
                      cg.randomize(rand);
                  
                  if(cg.anotGen.size()<k+1){
                      System.out.println("Nema dovoljno anotiranih COG-ova");
                      System.out.println(cg.anotGen.size());
                      continue;
                  }
              
                  if(useLogDist==0)//distance on circular genome
                    sim.computeSimilarity(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                      sim.computeLogSimilarity(cg, cgmap,geneOGMap, taxTranslation ,taxid);
                                
            }          
            else{
                   Gene cg=new Gene();
                   if(OrgType==0)
                     cg.findCogsNewnonCH(input, geneOGMap, map); 
                   if(OrgType==1)
                      cg.findCogsNewnonCHMetazoa(input, geneOGMap, map,taxid);
                   
                    if(cg.genesContig.size()>1000 && OrgType ==0){
                        System.out.println("To many contigs!");
                        continue;
                    }
                    
                    if(useLogDist==0)//distance on circular genome
                            sim.computeSimilaritynonCH(cg, cgmap, geneOGMap, taxTranslation,taxid); 
                  else if(useLogDist==1)//logarithm of distances
                            sim.computeLogSimilaritynonCH(cg, cgmap,geneOGMap, taxTranslation ,taxid);                   
                    }
            usedSpecies.add(taxTranslation.get(taxid));
                    }
            
            Iterator<String> it1 = sim.locationN.keySet().iterator();
            
            System.out.println("Before normalize");
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         sim.normalize(useLogDist);
         
         System.out.println("After normalize");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         
         sim.removeEmpty();
         
         System.out.println("After remove");
         it1 = sim.locationN.keySet().iterator();
            while(it1.hasNext()){
                String og = it1.next();
                Pair<ArrayList<Double>, ArrayList<Double>> p = sim.locationN.get(og);
                int countNum=0, countDist=0;
                for(int i=0;i<p.getValue0().size();i++){
                    if(p.getValue0().get(i)!=0)
                        countNum++;
                    if(p.getValue1().get(i)!=0)
                        countDist++;
                }
                
                System.out.println("cog: "+og+" CC: "+countNum+" CD: "+countDist);
            }
         
            System.out.println("cbf size: "+cbf.COGbaselinefeaturesmap.keySet().size());
            tt.createLocationTrainHeader(output, cgmap, sim ,headerInput, isTrain);
            tt.appendLocationRowsToSet(output, cgmap, sim, isTrain);
            cbf.COGbaselinefeaturesmap.clear();
        } catch (Exception e) {
            e.printStackTrace();
        }   
   }  
}
