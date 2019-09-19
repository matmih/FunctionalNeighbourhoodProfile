/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package semanticsimilarity;

import java.util.ArrayList;

/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to store information about GO ancestors
 */
public class Path {
    public ArrayList<Integer> path;
    
   public Path(){
        path=new ArrayList<>();
    }
    
  public  Path(Path p){
        path=new ArrayList<>();
        for(int i=0;i<p.path.size();i++)
            path.add(p.path.get(i));
    }
    
  public int findIntersection(Path p){
       int index=-1;
       
       for(int i=path.size()-1;i>=0;i--)
           for(int j=p.path.size()-1;j>=0;j--)
               if(p.path.get(j).compareTo(path.get(i))==0)
                   return path.get(i);         
       return index;
    }
    
    public int findIntersectionStart(Path p){
       int index=-1;
       
       for(int i=0;i<path.size();i++)
           for(int j=0;j<p.path.size();j++)
               if(p.path.get(j).compareTo(path.get(i))==0)
                   return path.get(i);         
       return index;
    }
    
}
