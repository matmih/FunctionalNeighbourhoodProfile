/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

import java.io.File;
import java.nio.file.Files;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class used to delete Files from a system
 */
public class FileDeleter {
    
    public void deleteFile(File toDelete){
        
        try{
            Files.delete(toDelete.toPath());
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
}
