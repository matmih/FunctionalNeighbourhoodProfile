/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package recursivefilesearch;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing function to count the number of annotated COGs in a .ptt file (use for prokaryotic genomes)
 */
public class ParsePttFile {
final static Charset ENCODING = StandardCharsets.UTF_8;

int countCOGannot(File input){
        int count=0,lineNum=0;

        BufferedReader reader;
         try {
      Path path =Paths.get(input.getAbsolutePath());
      System.out.println("Path: "+input.getAbsolutePath());
      reader = Files.newBufferedReader(path,ENCODING);
      String line = null;
      while ((line = reader.readLine()) != null) {
          lineNum++;
          if(lineNum>3){
          String[] st=line.split("\t");
          if(lineNum==4){
               for(int i=0;i<st.length;i++)
                   System.out.println(st[i]);
          }
          if(!st[7].equals("-"))
              count++;
          }
      }
         }catch(IOException ioe)
            {
              System.err.println("IOException: " + ioe.getMessage());
            }
       return count;
    }
}
