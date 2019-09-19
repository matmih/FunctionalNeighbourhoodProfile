package recursivefilesearch;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Locale;

/**
 * Performs the Fisher Exact Test to determine if there is a relationship 
 * between two categorical variables.
 * 
 * Adapted from JavaScript code by Oyvind Langsrud, see
 * <a href="http://www.langsrud.com/fisher.htm">
 * http://www.langsrud.com/fisher.htm</a>
 * 
 * @author Fran Supek (fran.supek[AT]irb.hr)
 */
public class FisherExactTest {


  /**
   * Reference: "Lanczos, C. 'A precision approximation
   * of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
   * Translation of  Alan Miller's FORTRAN-implementation
   * See http://lib.stat.cmu.edu/apstat/245
   */
  private double lnfact(int n) {

    if (n<=1) return(0.0);
    
    int z = n+1;
    //return(lngamm(n+1));

    double x = 0;
    x += 0.1659470187408462e-06/(z+7);
    x += 0.9934937113930748e-05/(z+6);
    x -= 0.1385710331296526    /(z+5);
    x += 12.50734324009056     /(z+4);
    x -= 176.6150291498386     /(z+3);
    x += 771.3234287757674     /(z+2);
    x -= 1259.139216722289     /(z+1);
    x += 676.5203681218835     /(z);
    x += 0.9999999999995183;
    return(Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));

  }

  //private double lnbico(double n, double k) {
  private double lnbico(int n, int k) {
    return( lnfact(n)-lnfact(k)-lnfact(n-k) );
  }


  
  private int sn11,sn1_,sn_1,sn;
  //private double sn11,sn1_,sn_1,sn;
  private double sprob;
  
  

  private double hyper0(int n11i, int n1_i, int n_1i, int ni) {
    if ( ! ( n1_i != 0 || n_1i != 0 || ni != 0 ) ) {
      if ( ! (n11i % 10 == 0) ) {
        if ( n11i==sn11+1 ) {
          sprob *= ( (double) (sn1_-sn11)/(n11i) )*( (double) (sn_1-sn11)/(n11i+sn-sn1_-sn_1) );
          sn11 = n11i;
          return sprob;
        }
        if (n11i==sn11-1 ) {
          sprob *= ( (double) (sn11)/(sn1_-n11i) )*( (double) (sn11+sn-sn1_-sn_1)/(sn_1-n11i) );
          sn11 = n11i;
          return sprob;
        }
      }
      sn11 = n11i;
    } else {
      sn11 = n11i;
      sn1_ = n1_i;
      sn_1 = n_1i;
      sn   = ni;
    }

    sprob = Math.exp(lnbico(sn1_,sn11)+lnbico(sn-sn1_,sn_1-sn11)-lnbico(sn,sn_1));

    return sprob;
  }

 

  public double left, right, twotail;
  
  final public double exact22(int n11_, int n12_, int n21_, int n22_) {

    if ( n11_<0 || n12_<0 || n21_<0 || n22_<0 )
      throw new IllegalArgumentException("Counts passed to exact22() must all" +
              " be nonnegative.");

    int n1_ = n11_ + n12_;
    int n_1 = n11_ + n21_;
    int n   = n11_ + n12_ + n21_ + n22_;

    //double prob = exact(n11_, n1_, n_1, n);
    double sleft, sright, sless, slarg;
    double prob;

    int i, j;

    int max = n1_;
    if (n_1<max) max=n_1;

    int min = n1_ + n_1 - n;
    if (min < 0) min=0;
    if (min == max) {

      sless = 1; sright= 1;
      sleft = 1; slarg = 1;
      prob = 1;

    } else {

      prob = hyper0(n11_,n1_,n_1,n);
      double probReduced = 0.99999999*prob;

      sleft = 0;
      double p = hyper0(min,0,0,0);
      for (i = min+1; p < probReduced; i++) {
        sleft += p;
        p = hyper0(i,0,0,0);
      }
      i--;
      if (p < 1.00000001*prob) sleft += p;
        else i--;

      sright = 0;
      p = hyper0(max,0,0,0);
      for (j = max-1; p < probReduced; j--) {
        sright += p;
        p = hyper0(j,0,0,0);
      }
      j++;
      if (p<1.00000001*prob) sright += p;
        else j++;

      if ( Math.abs(i-n11_) < Math.abs(j-n11_) ) {
        sless = sleft;
        slarg = 1 - sleft + prob;
      } else {
        sless = 1 - sright + prob;
        slarg = sright;
      }

    }

    this.left    = sless;
    this.right   = slarg;
    this.twotail = sleft + sright;
    if (twotail>1) twotail=1;
    return prob;
  }
  
  
  /**
   * Runs the Fisher's exact test for many 2x2 contingency tables given in a 
   * text file and reports the p-values (two-tailed), and also the odds ratio
   * with its 95% CI. Prints the original lines, with the additional columns appended.
   * 
   * @param args
   * @throws Exception 
   */
  
  public static void main(String[] args) throws Exception {
   
    /*
     * lines are "Specie_number     A    B    C    D"
     * A= Genes in recombined areas for some particular functional category
     * B= Genes outside recombined areas in a particular category
     * C= Genes in recombined areas for the rest of functional categories
     * D= Genes outside recombined areas for the rest of functional categories.
     */
    String inFile = "C:/Users/fsupek/Desktop/Annotation.Matrix.PedroOK.txt";
    
    BufferedReader bufRdr = new BufferedReader(new FileReader(inFile));
    String line = null;
    while ((line = bufRdr.readLine()) != null) {      
      line=line.trim();
      if ( line.isEmpty() ) {
        continue;
      }
      if ( line.startsWith("FILE")  ) {   // header
        System.out.printf( Locale.US, "%s\t%s\t%s\t%s\t%s\t%s\n", line, "odds_ratio", "OR_95%_conf_int", "Fisher_p_left_tail", "Fisher_p_right_tail", "Fisher_p_2tailed" );
        continue;
      }
      String[] cols = line.split("\t");
      // A, B, C and D are 4, 5, 6 and 7
      int a = Integer.parseInt(cols[4]), b = Integer.parseInt(cols[5]), c = Integer.parseInt(cols[6]), d = Integer.parseInt(cols[7]);
      FisherExactTest myFisher = new FisherExactTest();
      double myProb;
      try {
        myProb = myFisher.exact22(a, b, c, d);
      } catch (IllegalArgumentException e) {
        myProb = Double.NaN;
      }
      double OR = ((double)a/b)/((double)c/d);

      double SE_of_logOR = Math.sqrt( 1.0/a + 1.0/b + 1.0/c + 1.0/d ); 
      //double zSc = Math.log(OR) / SE_of_logOR; // could be used to find the p-value, if reqd
      
      double confInterval = 1.96 * SE_of_logOR ;
      if ( Double.isNaN(myProb) ) {
        System.out.printf("%s\t%.3f\t[%.3f-%.3f]\t%s\t%s\t%s\n", line, OR, Math.exp( Math.log(OR) - confInterval ), Math.exp( Math.log(OR) + confInterval ),
                "ERROR", "ERROR", "ERROR" );
      } else {
        System.out.printf("%s\t%.3f\t[%.3f-%.3f]\t%e\t%e\t%e\n", line, OR, Math.exp( Math.log(OR) - confInterval ), Math.exp( Math.log(OR) + confInterval ),
                myFisher.left, myFisher.right, myFisher.twotail );
      }
      
    }

    // some code for testing the class
    /*
    FisherExactTest myFisher;
    double myProb;
            
    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(5000, 3000, 2000, 1000);
    System.out.println("\nmyFisher.exact22(5000, 3000, 2000, 1000)");

    System.out.printf("Left p-value: %6.4f  Right p-value: %6.4f  " +
            "2-Tail p-value: %6.4f  Exact prob: %6.4f",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   

    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(500, 300, 200, 100);
    System.out.println("\nmyFisher.exact22(500, 300, 200, 100)");

    System.out.printf("Left p-value: %6.4f  Right p-value: %6.4f  " +
            "2-Tail p-value: %6.4f  Exact prob: %6.4f",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   
    
    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(50, 30, 20, 10);
    System.out.println("\nmyFisher.exact22(50, 30, 20, 10)");

    System.out.printf("Left p-value: %6.4f  Right p-value: %6.4f  " +
            "2-Tail p-value: %6.4f  Exact prob: %6.4f",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   

    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(5, 3, 2, 1);
    System.out.println("\nmyFisher.exact22(5, 3, 2, 1)");

    System.out.printf("Left p-value: %6.4f  Right p-value: %6.4f  " +
            "2-Tail p-value: %6.4f  Exact prob: %6.4f",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   
    
    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(87398, 402663, 4022, 10835);
    System.out.println("\nmyFisher.exact22(87398, 402663, 4022, 10835)");

    System.out.printf("Left p-value: %e  Right p-value: %e  " +
            "2-Tail p-value: %e  Exact prob: %e",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   
    
    myFisher = new FisherExactTest();
    myProb = myFisher.exact22(402663, 87398, 10835, 4022);
    System.out.println("\nmyFisher.exact22(402663, 87398, 10835, 4022)");

    System.out.printf("Left p-value: %e  Right p-value: %e  " +
            "2-Tail p-value: %e  Exact prob: %e\n",
            myFisher.left, myFisher.right, myFisher.twotail, myProb);   

    System.out.println("Starting speed test.");
    long millis = System.currentTimeMillis();
    double foolTheOptimizer = 0.0;

    for (int i = 0; i < 5e3; i++) {
       myFisher.exact22(40+i, 8+i, 1+i, 4+i);
       foolTheOptimizer += myFisher.right;
    }
    System.out.println("foolTheOptimizer: " + foolTheOptimizer + " Time(millis): " + (System.currentTimeMillis() - millis) );

    System.out.println("Starting speed test, run 2.");
    millis = System.currentTimeMillis();
    foolTheOptimizer = 0.0;

    for (int i = 0; i < 5e3; i++) {
       myFisher.exact22(40+i, 8+i, 1+i, 4+i);
       foolTheOptimizer += myFisher.right;
    }
    System.out.println("foolTheOptimizer: " + foolTheOptimizer + " Time(millis): " + (System.currentTimeMillis() - millis) );

    System.out.println("Starting speed test, run 2.");
    millis = System.currentTimeMillis();
    foolTheOptimizer = 0.0;

    for (int i = 0; i < 5e3; i++) {
       myFisher.exact22(40+i, 8+i, 1+i, 4+i);
       foolTheOptimizer += myFisher.right;
    }
    System.out.println("foolTheOptimizer: " + foolTheOptimizer + " Time(millis): " + (System.currentTimeMillis() - millis) );
    */
  }
}
