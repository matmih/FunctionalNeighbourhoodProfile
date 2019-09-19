/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package recursivefilesearch;

/**
 * @author Matej Mihelcic
  *@institution Rudjer Boskovic Institute, Zagreb, Croatia
  *@mail matmih1@gmail.com
  *@description class to compute the circular distance in prokaryotic genomes
 */
public class CircularDistance {
    
            double result=Double.MAX_VALUE;
            
         double computeDistance(double coordXL1,double coordXR1,double coordXL2,double coordXR2,double genomeSize){ 
                 
                    result=0.0;
                double max=genomeSize;
                if(coordXL1>=coordXL2 && coordXR1>=coordXR2){
                    if(coordXL1>coordXR2)
                        result+=Math.min(Math.abs(coordXL1-coordXR2),Math.abs(coordXR1-max-coordXL2))+Math.pow(10.0, -10.0);//added max
                    else
                        result+=Math.pow(10.0, -10.0);
                }
                else{
                    if(coordXL2>=coordXL1 && coordXR2>=coordXR1){
                        if(coordXL2>coordXR1)
                     result+=Math.min(Math.abs(coordXL2-coordXR1),Math.abs(coordXR2-max-coordXL1))+Math.pow(10.0, -10.0);
                        else
                            result+=Math.pow(10.0, -10.0);
                    }
                    else
                            result+=Math.pow(10.0, -10.0);
                }
                return result;
            }
         }
