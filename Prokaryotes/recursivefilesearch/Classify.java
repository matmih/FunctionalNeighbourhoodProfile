/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Random;
import weka.classifiers.Evaluation;
import weka.core.Instances;
import hr.irb.fastRandomForest.*;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import weka.classifiers.evaluation.output.prediction.PlainText;
import weka.classifiers.meta.FilteredClassifier;
import weka.core.Range;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.AddID;
import weka.filters.unsupervised.attribute.Remove;

 /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute single-class predictions using Fast Random Forest algorithm
  *https://github.com/sdvillal/fast-random-forest
 */
public class Classify {
     FileReader fr=null;
     
     public void classify(String input, String output, Random randGen, int run, int iter, int type, int seed,ArrayList<HashMap<Integer,ArrayList<Double>>> results){
         try{
         fr=new FileReader(input);
         }
         catch(Exception e){
             e.printStackTrace();
         }
          BufferedReader reader = new BufferedReader(fr);
          Instances data = null;
          try{
          data=new Instances(reader);
          reader.close();
          fr.close();
          }
          catch(Exception e){
              e.printStackTrace();
          }
         
 // setting class attribute
          data.setClassIndex(data.numAttributes() - 1);
         
          FastRandomForest frf=new FastRandomForest();
          
          frf.setSeed(seed/*randGen.nextInt((100 - 1) + 1) + 1*/);
          frf.setDebug(true);
          frf.setNumThreads(30);
          frf.setNumTrees(600);
          frf.setNumFeatures((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
          String[] options = new String[4];
          options[0]="-I";
          options[1]="5";
          options[2]="-K";
          options[3]=""+((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
          try{
          frf.setOptions(options);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          Evaluation eval=null;
          try{
            eval=new Evaluation(data);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          int classIndex=data.numAttributes()-1;
          
          try{
          eval.crossValidateModel(frf,data,5,randGen);
          }
          catch(Exception e){
              e.printStackTrace();
          }

          HashMap<Integer,ArrayList<Double>> values=results.get(iter);
          if(!values.containsKey(type))
              values.put(type, new ArrayList<Double>());
          
         double auprc=eval.areaUnderPRC(0);
         double auc=eval.areaUnderROC(0);
         values.get(type).add(auprc);
         double auprc1=eval.areaUnderPRC(1);
         double auc1=eval.areaUnderROC(1);
        
         try{
         FileWriter fw = new FileWriter(output);
         
          fw.write("AUC0: "+auc+"\n");
          fw.write("AUPRC0: "+auprc+"\n");
          fw.write("AUC1: "+auc1+"\n");
          fw.write("AUPRC1: "+auprc1+"\n");
          fw.write("Confusion matrix: \n\n");
          fw.write(eval.confusionMatrix()[0][0]+"\t"+eval.confusionMatrix()[0][1]+"\n");
          fw.write(eval.confusionMatrix()[1][0]+"\t"+eval.confusionMatrix()[1][1]+"\n");
         
         fw.close();
        }
        catch(Exception e)
                {e.printStackTrace();}  
     }
     
      public void classifyPred(String input, String output, Random randGen, int run, int iter, int type, int seed,ArrayList<HashMap<Integer,StringBuffer>> results){
         try{
         fr=new FileReader(input);
         }
         catch(Exception e){
             e.printStackTrace();
         }
          BufferedReader reader = new BufferedReader(fr);
          Instances data = null;
          try{
          data=new Instances(reader);
          reader.close();
          fr.close();
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          
           Remove remove = new Remove();
           remove.setAttributeIndices("1");
           
 // setting class attribute
          data.setClassIndex(data.numAttributes() - 1);
          
            AddID addId = new AddID();
            try{
            addId.setInputFormat(data); 
            data = Filter.useFilter(data, addId);
            }
            catch(Exception e){
                e.printStackTrace();
            }
         
          FastRandomForest frf=new FastRandomForest();
          
          frf.setSeed(seed/*randGen.nextInt((100 - 1) + 1) + 1*/);
          frf.setDebug(true);
          frf.setNumThreads(30);
          frf.setNumFeatures((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
          if(data.numAttributes()<2)
              frf.setNumFeatures(1);
          String[] options = new String[4];
          options[0]="-I";
          options[1]="5";
          options[2]="-K";
          options[3]=""+((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
       
          try{
          frf.setOptions(options);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          Evaluation eval=null;
          try{
            eval=new Evaluation(data);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          int classIndex=data.numAttributes()-1;
          PlainText pt=new PlainText();
          StringBuffer buff=new StringBuffer();
          pt.setBuffer(buff);
          weka.core.Range rwng = new weka.core.Range();
          rwng.setRanges("");
          java.lang.Boolean bl = false;
          ArrayList<Object> inputArr = new ArrayList<>();
          inputArr.add(pt); inputArr.add(rwng); inputArr.add(bl);
           FilteredClassifier fc = new FilteredClassifier();
          fc.setFilter(remove);
          fc.setClassifier(frf);
          Range attsToOutput = new Range("first-last");
          Boolean outputDist = new Boolean(true);

          try{
          eval.crossValidateModel(fc,data,5/*100*/,randGen,pt,attsToOutput,outputDist);
          }
          catch(Exception e){
              e.printStackTrace();
          }

          HashMap<Integer,StringBuffer> values=results.get(iter);
          if(!values.containsKey(type))
              values.put(type, pt.getBuffer());
          
          System.out.println("Predictions: ");
          System.out.println(values.get(type).toString());
          
         double auprc=eval.areaUnderPRC(0);
         double auc=eval.areaUnderROC(0);
         double auprc1=eval.areaUnderPRC(1);
         double auc1=eval.areaUnderROC(1); 
     }
     
      public void classifyCMSum(String GO, GOMap mapFull,String input, String output, Random randGen, int run, int iter, int type, int seed,HashMap<Integer,ArrayList<Double>> results){
         try{
         fr=new FileReader(input);
         }
         catch(Exception e){
             e.printStackTrace();
         }
          BufferedReader reader = new BufferedReader(fr);
          Instances data = null;
          try{
          data=new Instances(reader);
          reader.close();
          fr.close();
          }
          catch(Exception e){
              e.printStackTrace();
          }
         
 // setting class attribute
          data.setClassIndex(data.numAttributes() - 1);
         
          FastRandomForest frf=new FastRandomForest();
          
          frf.setSeed(seed/*randGen.nextInt((100 - 1) + 1) + 1*/);
          frf.setDebug(true);
          frf.setNumThreads(25);
          frf.setNumTrees(600);
          frf.setNumFeatures((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
          String[] options = new String[4];
          options[0]="-I";
          options[1]="5";
          options[2]="-K";
          options[3]=""+((int)Math.floor(Math.sqrt(data.numAttributes()))+1);
        
          try{
          frf.setOptions(options);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          Evaluation eval=null;
          try{
            eval=new Evaluation(data);
          }
          catch(Exception e){
              e.printStackTrace();
          }
          
          int classIndex=data.numAttributes()-1;
          
          try{
          eval.crossValidateModel(frf,data,5,randGen);
          }
          catch(Exception e){
              e.printStackTrace();
          }
         
         if(!results.containsKey(mapFull.GOmap.get(GO)))
         {
             ArrayList<Double> tmpAr=new ArrayList<>();
            tmpAr.add(eval.confusionMatrix()[0][0]);
            tmpAr.add(eval.confusionMatrix()[0][1]);
            tmpAr.add(eval.confusionMatrix()[1][0]);
            tmpAr.add(eval.confusionMatrix()[1][1]);
             results.put(mapFull.GOmap.get(GO), tmpAr);
         }
         else{
             ArrayList<Double> tmpAr=results.get(mapFull.GOmap.get(GO));
             tmpAr.add(eval.confusionMatrix()[0][0]);
             tmpAr.add(eval.confusionMatrix()[0][1]);
             tmpAr.add(eval.confusionMatrix()[1][0]);
             tmpAr.add(eval.confusionMatrix()[1][1]);
             results.put(mapFull.GOmap.get(GO), tmpAr); 
         }
     }    
}
