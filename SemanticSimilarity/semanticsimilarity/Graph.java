package semanticsimilarity;


import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class containing structure and functions to deal with Graphs
 */
public class Graph {
    
    int adjacency[][];
    ArrayList<ArrayList<Integer>> adjList=new ArrayList<ArrayList<Integer>>();
    ArrayList<ArrayList<Integer>> childList=new ArrayList<ArrayList<Integer>>();
    ArrayList<HashSet<Integer>> dagList=new ArrayList<HashSet<Integer>>();
    int numGO;
    
    public Graph(int numGOs){
        adjacency=new int[numGOs][numGOs];
        numGO=numGOs;
    }
 
    public void createAdjacency(GOMap map){
        for(int i=0;i<numGO;i++){
            adjList.add(i,new ArrayList<Integer>());
            childList.add(i,new ArrayList<Integer>());
            dagList.add(i,new HashSet<Integer>());
        }
        
            String f[]=map.hierarchy.split(",");

            int index[]=new int[2];
            for(int i=0;i<f.length;i++){
                String fun[]=f[i].split("/");
                for(int j=0;j<fun.length;j++){
                    if(fun[j].length()<7 && !fun[j].contains("root")){
                    while(fun[j].length()<7)
                        fun[j]="0"+fun[j];
                    if(!fun[j].contains("GO"))
                    fun[j]="GO"+fun[j];
                    index[j]=map.GOmap.get(fun[j]);
                    }
                    else if(!fun[j].contains("GO") && !fun[j].contains("root")){
                       fun[j]="GO"+fun[j];
                    index[j]=map.GOmap.get(fun[j]); 
                    }
                    else if(fun[j].contains("GO") && !fun[j].contains("root") && fun[j].length()>=9){
                         index[j]=map.GOmap.get(fun[j]); 
                    }
                    if(j==1){
                     adjacency[index[1]][index[0]]=1;//reverse the graph
                    adjList.get(index[1]).add(index[0]);//reverse the graph
                    childList.get(index[0]).add(index[1]);
                    dagList.get(index[0]).add(index[0]);
                     dagList.get(index[1]).add(index[1]);
                     dagList.get(index[1]).add(index[0]);
                    }
                }
            }
            
        }
    
     public void createAdjacency1(GOMap map){
        
        for(int i=0;i<numGO;i++){
            adjList.add(i,new ArrayList<Integer>());
            childList.add(i,new ArrayList<Integer>());
            dagList.add(i,new HashSet<Integer>());
        }
        
            String f[]=map.hierarchy.split(",");

            int index[]=new int[2];
            for(int i=0;i<f.length;i++){
                String fun[]=f[i].split("/");
                for(int j=0;j<fun.length;j++){
                    if(!fun[j].contains("root")){
                        if(!map.GOmap.containsKey(fun[j]))
                            continue;
                    index[j]=map.GOmap.get(fun[j]);
                    }
                    else if(!fun[j].contains("root")){
                       fun[j]=fun[j];
                    index[j]=map.GOmap.get(fun[j]); 
                    }
                    if(j==1){
                     adjacency[index[1]][index[0]]=1;//reverse the graph

                    adjList.get(index[1]).add(index[0]);//reverse the graph
                    childList.get(index[0]).add(index[1]);
                    dagList.get(index[0]).add(index[0]);
                     dagList.get(index[1]).add(index[1]);
                     dagList.get(index[1]).add(index[0]);
                    }
                }
            }       
        }
    
   public void saveAdjacency(File output, int numGO){
        
        try
            {
                FileWriter fw = new FileWriter(output);
                
                fw.write("D\n");
                fw.write(numGO+"\n");
                
                 for(int i=0;i<numGO;i++){
            for(int j=0;j<numGO;j++)
                fw.write(adjacency[i][j]+" ");
            fw.write("\n");
                 }          
                fw.close();
            }
           catch(Exception e){
                e.printStackTrace();
            } 
    }
    
    ArrayList<Integer> getAdjNodes(int node){
        ArrayList<Integer> res=new ArrayList<>();
        
        for(int i=0;i<numGO;i++)
            if(adjacency[node][i]==1 && node!=i)
                res.add(i);
        
        return res;
    }
    
   public boolean isConnected(int node1, int node2){
        
        return (adjacency[node1][node2]==1);
    }
    
   public boolean Contains(Path p, int node){
        for(int i=0;i<p.path.size();i++)
            if(p.path.get(i).compareTo(node)==0)
                return true;
        return false;
    }
    
   public void getAllPaths(int startNodeId, int endNodeId, ArrayList<Path> paths, Path visited){
        
        int back = visited.path.get(visited.path.size()-1);

    ArrayList<Integer> adjNode = adjList.get(back);

    // Examine adjacent nodes
    for ( int i=0;i<adjNode.size();i++)
    {   
        int node = adjNode.get(i);
        dagList.get(startNodeId).add(node);
        if ( Contains( visited, node ) ) continue;
 
        if ( node == endNodeId )
        {
            visited.path.add(node);
            // Get hop count for this path
            int size = visited.path.size();
            int hops = size - 1; 
          
                Path path=new Path( visited );
                paths.add(path );          
            visited.path.remove(hops);
        }           
    }
 
 
    // in breadth-first, recursion needs to come after visiting adjacent nodes
    for ( int i=0;i<adjNode.size();i++)
    { 
        int node = adjNode.get(i);
       if ( Contains( visited, node )  || node == endNodeId )     
            continue;
 
        visited.path.add(node);
 
        getAllPaths( startNodeId,endNodeId,paths, visited);        
 
        int n = (int) visited.path.size() - 1;     
        visited.path.remove(n);
    }       
         
    }
    
    
    public double[] computeSimilarity(int GO1, int GO2, GOMap map){
        double sim=0.0;
        
        
        if(GO1==GO2){
            double resE[]=new double[3];
           sim=map.computeScore(GO1);
            resE[0]=sim;
            resE[1]=2*sim/(map.computeScore(GO1)+map.computeScore(GO2));
            resE[2]=2*sim*(1-map.frequency.get(GO1))/(map.computeScore(GO1)+map.computeScore(GO2));
            return resE;
        }
        
        ArrayList<Path> paths1=new ArrayList<>();
         ArrayList<Path> paths2=new ArrayList<>();
        Path visited=new Path();

       visited.path.add(GO1);
        getAllPaths(GO1, 0, paths1, visited);
       visited.path.clear();
       visited.path.add(GO2);
       getAllPaths(GO2,0,paths2,visited);
       visited.path.clear();
       
       double tsim=0.0;
       double res[]=new double[3];
       int Gind=0;
       for(int i=0;i<paths1.size();i++){
           for(int j=0;j<paths2.size();j++){  

                int Cindex=paths2.get(j).findIntersectionStart(paths1.get(i));

                tsim=map.computeScore(Cindex);
				
                if(tsim>sim){
                    sim=tsim;
                 Gind=Cindex;
                }
              }
       }   
       res[0]=sim;
       res[1]=2*sim/(map.computeScore(GO1)+map.computeScore(GO2));
       res[2]=2*sim*(1-map.frequency.get(Gind))/(map.computeScore(GO1)+map.computeScore(GO2));
        return res;
    }
    
   public double computeSimilarityLin(int GO1, int GO2, GOMap map){
        double sim=0.0;
        
        
        if(GO1==GO2)
            return map.computeScore(GO1);
        
        ArrayList<Path> paths1=new ArrayList<>();
         ArrayList<Path> paths2=new ArrayList<>();
        Path visited=new Path();
       visited.path.add(GO1);
        getAllPaths(GO1, 0, paths1, visited);
       visited.path.clear();
       visited.path.add(GO2);
       getAllPaths(GO2,0,paths2,visited);
       visited.path.clear();
       
       double tsim=0.0;
       for(int i=0;i<paths1.size();i++){
           for(int j=0;j<paths2.size();j++){  
                int Cindex=paths2.get(j).findIntersectionStart(paths1.get(i));
                tsim=map.computeScore(Cindex);
                if(tsim>sim)
                    sim=tsim;
              }
       }          
        return 2*sim/(map.computeScore(GO1)+map.computeScore(GO2));
    }
    
    public double computeSimilarityRel(int GO1, int GO2, GOMap map){
        double sim=0.0;
        
        
        if(GO1==GO2)
            return map.computeScore(GO1);
        
        ArrayList<Path> paths1=new ArrayList<>();
         ArrayList<Path> paths2=new ArrayList<>();
        Path visited=new Path();
       visited.path.add(GO1);
        getAllPaths(GO1, 0, paths1, visited);
       visited.path.clear();
       visited.path.add(GO2);
       getAllPaths(GO2,0,paths2,visited);
       visited.path.clear();
       
       double tsim=0.0;
       int Gind=0;
       for(int i=0;i<paths1.size();i++){
           for(int j=0;j<paths2.size();j++){  
                int Cindex=paths2.get(j).findIntersectionStart(paths1.get(i));
                tsim=map.computeScore(Cindex);
                if(tsim>sim){
                    sim=tsim;
                    Gind=Cindex;
                }
              }
       }          
        return 2*sim*(1-map.frequency.get(Gind))/(map.computeScore(GO1)+map.computeScore(GO2));
    }
    
}
