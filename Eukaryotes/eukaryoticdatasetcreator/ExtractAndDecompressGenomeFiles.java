/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eukaryoticdatasetcreator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import org.apache.commons.io.FileUtils;
import org.rauschig.jarchivelib.Archiver;
import org.rauschig.jarchivelib.ArchiverFactory;
/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to extract and decompress eukaryotic genome files
 */
public class ExtractAndDecompressGenomeFiles {
    
    public static void gunzipIt(String INPUT_GZIP_FILE, String OUTPUT_FILE){

     byte[] buffer = new byte[1024];
     
     
     try{

    	 GZIPInputStream gzis =
    		new GZIPInputStream(new FileInputStream(INPUT_GZIP_FILE));

    	 FileOutputStream out =
            new FileOutputStream(OUTPUT_FILE);

        int len;
        while ((len = gzis.read(buffer)) > 0) {
        	out.write(buffer, 0, len);
        }

        gzis.close();
    	out.close();

    	System.out.println("Done");

    }catch(IOException ex){
       ex.printStackTrace();
    }
   }
    
    public static void main(String [] args){
        
        //String root="C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiGenomes";
        String root="C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\AnimalGenomesHigh";
        File folder=new File(root);
        String extensions[]={"gz"};
        boolean recursive=true;

        //File destination = new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\FungiExtracted");
        File destination = new File("C:\\Users\\matej\\Downloads\\Eukaryot data\\Metazoa\\AnimalHighExtracted");
         
        try{
         Collection files = FileUtils.listFiles(folder, extensions, recursive);
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File input = (File) iterator.next();
                System.out.println("File = " + input.getAbsolutePath());
                 System.out.println("File = " + input.getName());
                String outName=input.getName();
                outName=outName.substring(0, outName.length()-3);
                 System.out.println("OutName: "+outName);
               gunzipIt(input.getAbsolutePath(),destination.getAbsolutePath()+"\\"+outName);
                
                /*Archiver archiver = ArchiverFactory.createArchiver(input);
                archiver.extract(input, destination);*/
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
}
