/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.Collection;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;

  /**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description this class is used to copy all the .ptt files from the gene directory to a working local directory
 * (can be used to copy only .ptt files with a given number of annotated COGs)
 */
 
public class Copyptt {

     void copy(File root, String dest, String[] extensions, boolean recursive){

         try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);

            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File file = (File) iterator.next();
                System.out.println("File = " + file.getAbsolutePath());
                 System.out.println("File = " + file.getName());
                File sourceLocation=file,targetLocation=new File(dest+file.getName());
                InputStream in = new FileInputStream(sourceLocation);
            OutputStream out = new FileOutputStream(targetLocation);

            // Copy the bits from instream to outstream
            byte[] buf = new byte[1024];
            int len;
            while ((len = in.read(buf)) > 0) {
                out.write(buf, 0, len);
            }
            in.close();
            out.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
     }
     
     void TaxIDFile(File root, String dest, String[] extensions, boolean recursive){
         try{
         Collection files = FileUtils.listFiles(root, extensions, recursive);
         File targetLocation=new File(dest+"TaxIDs.txt");
         BufferedWriter out = new BufferedWriter(new OutputStreamWriter(
          new FileOutputStream(dest+"TaxIDs.txt"), "utf-8"));
         
            for (Iterator iterator = files.iterator(); iterator.hasNext();) {
                File file = (File) iterator.next();
                System.out.println("File = " + file.getAbsolutePath());
                 System.out.println("File = " + file.getName());
                File sourceLocation=file;
                InputStream in = new FileInputStream(sourceLocation);
            TaxIDExtractor ex=new TaxIDExtractor();
            ex.ExtractID(sourceLocation);
            out.write(file.getName()+" "+ex.TaxID+"\n");
            in.close();
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
     }
}
