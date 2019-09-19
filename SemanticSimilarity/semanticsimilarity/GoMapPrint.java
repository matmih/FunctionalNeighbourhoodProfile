/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package semanticsimilarity;

import java.io.File;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to print GO map
 */
public class GoMapPrint {
        public static void main(String[] args) {
        // TODO code application logic here
        File input=new File("C:\\Users\\matej\\Documents\\NetBeansProjects\\Plot error measures\\Data\\PhyleticProfiles-5-Binary-Hierarchical-reduced-10.arff");
        GOMap map=new GOMap();
        map.CreateGOMap(input);
        map.printGOMapNumeric();
        map.saveGOMap(new File("C:\\Users\\matej\\Desktop\\GOMap.txt"));
            
        }
}
