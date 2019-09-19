package OntologyTools;

//import hr.irb.geneMonkey.geneGroups.GenePlex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 * This class containes removed functions from the Revigo's TermSimilarityMatrix;
 * in particular, this concerns the code that can (a) find the similarity matrix
 * using the FunSimMat server via XML-RPC, and (b) find the similarity matrix from
 * gene ovelap in the GenePlex class.
 * 
 * The DistanceMatrix class is a nested HashMap. It is constructed by 
 * sending a provided Collection of GOTerms to the FunSimMat server. 
 * 
 * Note: the server may not support some of the GOTerms; all terms that the
 * server does support (and, consequentially, that the TermSimilarityMatrix
 * contains) are stored in the supportedTerms list.
 * 
 * @author Fran Supek
 */
public class TermSimilarityMatrix_Backup extends HashMap < Integer , HashMap < Integer, Double > > {

  /** All terms supported by the FunSimMat server. */
  public GoTermList supportedTerms;
  
  
  /**
   * Constructs a distance matrix of all the GO terms present in termList.
   * This is done by making an XML-RPC call to the FunSimMat database and
   * asking for the specified simScore. <p>
   * 
   * NOT SUPPORTED ANYMORE, to avoid the dependance on the Redstone library.
   * 
   * @return A DistanceMatrix object.
   */
  /*
  public TermSimilarityMatrix (
          Collection<GOTerm> termList,
          SemanticSimilarityScore simScore,
          GeneOntology myGO,
          String urlOfServer) throws Exception {
      
    throw new IllegalArgumentException("This constructor not supported.");
      
    /*
    if ( termList.size() == 0 ) {
      //throw new IllegalStateException( "A termList of size 0 was passed to " +
      //        "the TermSimilarityMatrix constructor." );
      this.supportedTerms = new GoTermList(); // empty list
      return;
    }

    // we need one parameter for the RPC call - a String containg a
    // comma-separated list of GO cats
    StringBuilder sb = new StringBuilder();
    for ( GOTerm aGoTerm : termList ) {
      sb.append( String.format( "GO:%07d,", aGoTerm.getId() )  );
    }
    Vector rpcParams = new Vector();
    rpcParams.add( sb.toString().substring(0, sb.toString().length()-1) );

    // now contact the server
    Object rpcResult;
    try {
      XmlRpcClient client = new XmlRpcClient( urlOfServer, true );
      rpcResult = client.invoke("Semantic.getSemSims", rpcParams);       
    } catch ( XmlRpcException ex ) {
      System.err.println("Redstone XML-RPC API threw an exception:" + ex);

      return;
    } catch ( XmlRpcFault ex ) {
      System.err.println("Redstone XML-RPC API threw an exception:" + ex);
      return;
    }

    XmlRpcArray rpcArray;
    
    // rpcResult contains an array of arrays of strings
    try {
      rpcArray = (XmlRpcArray) rpcResult;
      //rpcArray = (XmlRpcArray) myParser.parsedValue;
      XmlRpcArray tempArray = rpcArray.getArray(0);
      
      //for ( int i = 0; i<tempArray.size(); i++ )
      //  System.out.print( tempArray.getString(i) + "\t" );
      //System.out.print("\n");

      if ( tempArray.size() != 6 ) throw new XmlRpcException("First row" +
              "does not contain six members, as it should.");
    } catch ( Exception ex ) {
      System.err.print( "Result of unusual format returned from FunSimMat" +
              "server. Exception: " + ex ); return;
    }
  
    // go through table line by line and fill the distMatrix
    // also create a set of all terms returned by server (some might not be
    // supported if FunSimMat server uses an older version of the Gene Ontology,
    // where those terms might not exist
    Set<GOTerm> supportedTermSet = new HashSet<GOTerm>();
    
    for (int i = 1; i < rpcArray.size(); i++) {  // first item contains header
      XmlRpcArray curRow = rpcArray.getArray(i);

      Integer GO1 = Integer.parseInt( curRow.getString(0).replace("GO:", "") );
      Integer GO2 = Integer.parseInt( curRow.getString(1).replace("GO:", "") );
      double similarity = Double.parseDouble( curRow.getString(
              simScore.getColumnIndex() ) );

      // distances are symmetrical
      if ( ! this.containsKey( GO1 ) )
        this.put(GO1, new HashMap<Integer, Double>() );
      this.get( GO1 ).put(GO2, similarity);
      if ( ! this.containsKey( GO2 ) )
        this.put(GO2, new HashMap<Integer, Double>() );
      this.get( GO2 ).put(GO1, similarity);
      
      if ( !myGO.containsKey(GO1) ||
              !myGO.containsKey(GO2) )
        throw new Exception( String.format("The FunSimMat server has returned a" +
                " pair of GO Term IDs (%d,%d) among which at least one is not " +
                "present among the known GO terms", GO1, GO2));
      
      supportedTermSet.add( myGO.get(GO1) );
      supportedTermSet.add( myGO.get(GO2) );
      
    }

    // Report if some terms were not properly handled by FunSimMat server.
    for ( GOTerm aGoTerm : termList )
      if ( !supportedTermSet.contains( aGoTerm ) )
        System.err.printf("FunSimMat server does not recognize GO:" +
                "%s. Term will be omitted in the output CSV file.\n",
                aGoTerm.toString() );
    
    // Continue with only the supported terms 
    this.supportedTerms = new GoTermList(supportedTermSet);
        
  }
  */

  
  
  

  /**
   * Constructs a distance matrix of all the GO terms present in termList.
   * This is done by using the goCategorySimilarity() function of the
   * supplied GenePlex, which must be constructed beforehand.
   *
   * @return A DistanceMatrix object.
   */
   /*  
  public TermSimilarityMatrix (
          Collection<GOTerm> termList,
          //GenePlex genePlex,
          Object genePlex,
          GeneOntology myGO) throws Exception {

   throw new IllegalArgumentException("This constructor not supported.");
 
    if ( termList.size() == 0 )
      throw new IllegalStateException( "A termList of size 0 was passed to " +
              "the TermSimilarityMatrix constructor." );

    // go through table line by line and fill the distMatrix

    Set<GOTerm> supportedTermSet = new HashSet<GOTerm>();

    GOTerm[] termArr = termList.toArray(new GOTerm[0]);

    System.out.println( Arrays.toString(termArr) );

    for (int i = 0; i < termArr.length; i++) {

      for (int j = i; j < termArr.length; j++  ) {

        double simil;
        Integer go1 = termArr[i].getId();
        Integer go2 = termArr[j].getId();

        //System.out.printf( "term1: %d, term2: %d, ", go1, go2 );

        if ( i == j ) {
          simil = 1.0;
        } else {
          simil = genePlex.goCategorySimilarity( go1, go2 );
        }

        // distances are symmetrical
        if ( ! this.containsKey( go1 ) )
          this.put(go1, new HashMap<Integer, Double>() );
        this.get( go1 ).put(go2, simil);
        if ( ! this.containsKey( go2 ) )
          this.put(go2, new HashMap<Integer, Double>() );
        this.get( go2 ).put(go1, simil);

        supportedTermSet.add( myGO.get(go1) );

        //System.out.printf( "similarity: %f\n", simil );

      }

    }

    // All terms will normally be supported by the GenePlex; if some of them
    // aren't, an exception will be thrown before this point is reached.
    this.supportedTerms = new GoTermList(supportedTermSet);

    //System.out.println("Done with creation of TermSimilarityMatrix.");
  }
  */




  /**
   * Constructs a distance matrix of all the GO terms present in termList.
   * This is done by using the routines in SemanticSimilarityScore enum, but
   * a GoTermSizes object must also be provided in this case.
   *
   * @return A DistanceMatrix object.
   */
  /*
  public TermSimilarityMatrix (
          Collection<GOTerm> termList,
          SemanticSimilarityScore simScore,
          GeneOntology myGO,
          GoTermSizes sizes) throws Exception {

    if ( termList.size() == 0 )
      throw new IllegalStateException( "A termList of size 0 was passed to " +
              "the TermSimilarityMatrix constructor." );

    // go through table line by line and fill the distMatrix 

    Set<GOTerm> supportedTermSet = new HashSet<GOTerm>();

    GOTerm[] termArr = termList.toArray(new GOTerm[0]);

    //System.out.println( Arrays.toString(termArr) );

    for (int i = 0; i < termArr.length; i++) {

      for (int j = i; j < termArr.length; j++  ) {

        double simil;
        Integer go1 = termArr[i].getId();
        Integer go2 = termArr[j].getId();

        //System.out.printf( "term1: %d, term2: %d, ", go1, go2 );
        simil = simScore.calculateDistance(termArr[i], termArr[j], sizes, myGO);

        // distances are symmetrical
        if ( ! this.containsKey( go1 ) )
          this.put(go1, new HashMap<Integer, Double>() );
        this.get( go1 ).put(go2, simil);
        if ( ! this.containsKey( go2 ) )
          this.put(go2, new HashMap<Integer, Double>() );
        this.get( go2 ).put(go1, simil);

        supportedTermSet.add( myGO.get(go1) );

        //System.out.printf( "similarity: %f\n", simil );

      }

    }

    // All terms will normally be supported by the GoTermSizes; if some of them
    // aren't, an exception will be thrown before this point is reached. 
    this.supportedTerms = new GoTermList(supportedTermSet);

    //System.out.println("Done with creation of TermSimilarityMatrix.");

  }
  */

  /**
   * Uses this TermSimilaryMatrix to calculate "uniqueness" of each GO term,
   * i.e. its average semantic distance from all other GO terms in the matrix.
   * <p>
   * The "commonness" will be stored as a property in the provided
   * GoTermProperties object (this will overwrite old values of "uniqueness",
   * if present).
   *
   * @param termProps
   */
  public void calculateUniqueness( GoTermProperties termProps ) {

    for ( Integer curGoId : this.keySet() ) {

      double sum = 0.0;
      int count = 0;

      Map < Integer, Double > innerMap = this.get(curGoId);

      for ( Integer curGoIdInner : innerMap.keySet() ) {
        if ( curGoId == curGoIdInner )
          continue;
        if ( innerMap.get(curGoIdInner) != null && !Double.isNaN(innerMap.get(curGoIdInner)) ) {
          sum += innerMap.get(curGoIdInner);
          count++;
        }
      }

      if ( innerMap.size() <= 1 ) {
        termProps.put( curGoId, "uniqueness", 1.0 );
      } else {
        // this assumes the semantic similarity measure in the matrix ranges
        // from 0 to 1, so that 1-x can be used to convert 'similarity' to
        // 'distance'
        termProps.put( curGoId, "uniqueness",
                Math.pow( 1.0 - (sum / (double) (count-1)) , 2.0 ) );
      }

    }

  }





  /**
   * Determines the null distribution for a GeneOntology distance metric.
   * 
   * Returns an array of i members (i=100), where each index represents the 
   * value of the distance metric at the i-th percentile of the null 
   * distribution. E.g. the result[95] returns the 95th percentile.
   * 
   * @param topLevel All randomly selected terms which will be used in 
   * generating the null distribution must be children of this top level node.
   * @param sampleSize How many terms to use in sampling. Over 300 is not 
   * recommended when using FunSimMat server, and will throw an
   * IllegalArgumentException.
   * @param prokaryoticOnly Consider only terms valid in prokaryotes.
   * 
   */
  /*
  public static double[] nullDistribution(GoNamespace namespace, int sampleSize,
          SemanticSimilarityScore similScore, boolean prokaryoticOnly)
          throws Exception {
    
    if ( sampleSize > 300 ) 
      throw new IllegalArgumentException("Sample sizes greater than 300 are" +
              "not supported.");
    
    GeneOntology myGO = new GeneOntology();
    
    Set<GOTerm> terms = new HashSet<GOTerm>();
    while ( terms.size() < sampleSize ) {
      GOTerm aTerm = myGO.getRandomNode();
      if ( ! namespace.isParentOf(aTerm) )
        continue;
      if ( prokaryoticOnly && ! aTerm.isInProkaryoticSubset() )
        continue;
      terms.add( aTerm );
    }
    

    TermSimilarityMatrix myMatrix = new TermSimilarityMatrix( terms,
            similScore, myGO, "http://funsimmat.bioinf.mpi-inf.mpg.de/xmlrpc.php" );
    Collection<GOTerm> termList =  myMatrix.supportedTerms;    
    
    List<Double> simList = new ArrayList<Double>( sampleSize*sampleSize );
    for (GOTerm curGoTerm : termList) {
      int outerGoId = curGoTerm.getId();
      for (GOTerm innerGoTerm : termList) {
        int innerGoId = innerGoTerm.getId();
        if ( innerGoId == outerGoId ) continue;
        simList.add( myMatrix.get( outerGoId ).get( innerGoId ) );        
      }
    }

    Double[] simArray = simList.toArray(new Double[simList.size()]);
    Arrays.sort(simArray);
    
    double[] simPercentiles = new double[101];
    for ( int i = 0; i <= 100; i++ ) {
      int idxInSimArray = (int) ((double) i / 100 * (simArray.length - 1) );
      simPercentiles[i] = simArray[ idxInSimArray ]; 
    }
    
    return simPercentiles;
    
    
  }
  */
  


  /**
   * Filters the rows and columns in the matrix by the specified property.
   * Warning, this alters the original TermSimilarityMatrix!
   */
  protected void filterByProperty( GoTermProperties termProperties,
          String propertyForFiltering,
          double propertyLowBound, double propertyHighBound ){

    Iterator< Integer > itOuter =
            this.keySet().iterator();

    Map<Integer, GOTerm> termCache = new HashMap<Integer, GOTerm>();
    for ( GOTerm term : supportedTerms )
      termCache.put(term.getId(), term);

    while ( itOuter.hasNext() ) {

      Integer curKey = itOuter.next();

      HashMap< Integer, Double > curValue = this.get( curKey );

      Double propValueForTerm = termProperties.get( termCache.get(curKey), propertyForFiltering );
      if ( propValueForTerm != null
              &&
              ( propValueForTerm < propertyLowBound || propValueForTerm > propertyHighBound )
              ) { //remove it!
        itOuter.remove();
        continue;
      }

      Iterator< Integer > itInner = curValue.keySet().iterator();

      while ( itInner.hasNext() ) {

        Integer curKeyInner = itInner.next();

        Double propValueForInnerTerm = termProperties.get( termCache.get(curKeyInner), propertyForFiltering );
        
        if ( propValueForInnerTerm != null
                && 
                ( propValueForInnerTerm < propertyLowBound || propValueForInnerTerm > propertyHighBound )
                ) { //remove it!
          itInner.remove();
        }

      } // itInner loop

    } // itOuter loop

    // also remove all dispensable terms from "supportedTerms"
    Iterator<GOTerm> termsIter = supportedTerms.iterator();
    while (termsIter.hasNext()) {
      Double curTermPropertyVal = termProperties.get(termsIter.next(), propertyForFiltering);
      if ( curTermPropertyVal < propertyLowBound || curTermPropertyVal > propertyHighBound )
        termsIter.remove();
    }

  }

  
  
  /**
   * Converts the similarity matrix into a Weka Instances object.
   * 
   * @param datasetDesc A description to be put into header of arff file.
   * @param exponent Raise each distance in the newly created
   * weka Instances to this power (does not mess with original data in similarity
   * matrix).
   * @return Weka Instances with the distance matrix.
   */
  protected Instances toWekaInstances(String datasetDesc, double exponent ) {
    
    Instances myInst;
    double[] vals;
    FastVector atts = new FastVector();
        
    for (int i = 0; i < supportedTerms.size(); i++) {
      Attribute curAtt = new Attribute( 
              String.format( "dist_%07d", supportedTerms.get(i).getId() ) );
      atts.addElement(curAtt);
    }

    myInst = new Instances( datasetDesc, atts, supportedTerms.size() );

    for (int i = 0; i < supportedTerms.size(); i++) {  // instance by instance
    
      vals = new double[ supportedTerms.size() ]; // as many attributes as instances
      for (int j = 0; j < supportedTerms.size(); j++) {  // attribute by attribute 
        vals[j] = this.get(
                supportedTerms.get(i).getId() ).get( supportedTerms.get(j).getId() );
        vals[j] = Math.pow(vals[j], exponent);
      }
      // add instance to the dataset
      Instance curInst = new DenseInstance(1.0, vals);
      myInst.add( curInst );
    }

    myInst.setClassIndex(-1);  // no class
    return myInst;
    
  }  
  

  /**
   * Converts the similarity matrix into a Weka Instances object.
   * 
   * @return Weka Instances with the distance matrix.
   */
  protected Instances toWekaInstances() {
    return this.toWekaInstances("", 1.0);
  }
  
  
  
  
  /**
   * A test routine to see if the Redstone XML-RPC API functions correctly.
   * Uses the WordTracker test server. <p>
   * 
   * NOT SUPPORTED ANYMORE to remove the dependency on the Redstone API.
   */
  @Deprecated
  private void testRedstoneAPI() throws Exception  {
    /*XmlRpcClient client = new XmlRpcClient(
            "http://test.xmlrpc.wordtracker.com/", true);
    // keyphrases parameter
    Vector keyphrases = new Vector();
    keyphrases.add("mp3"); keyphrases.add("britney spears");

    // start building the parameter list
    Vector params = new Vector();
    params.add("guest");
    params.add(keyphrases);
    params.add("case_distinct");
    params.add(Boolean.FALSE);  params.add(Boolean.TRUE);
    params.add("exclude_adult");
    params.add(new Integer(100)); params.add(new Integer(10));

    Object rpcResult = client.invoke("get_exact_phrase_popularity", params);
    System.out.println(rpcResult.getClass());
    System.out.println(rpcResult);*/
  }  
  
  
  
  
  public static void main(String[] args) throws Exception {
    
    
    /* Finds null distributions for the SimRel semantic similarity score. */
    /*
    double[] percentiles;
    
    percentiles = nullDistribution(GoNamespace.BIOLOGICAL_PROCESS, 200,
            SemanticSimilarityScore.SIMREL, true);
    System.out.println( "Biological process, 200 samples, SimRel metric:\n"
            + Arrays.toString(percentiles) );
    percentiles = nullDistribution(GoNamespace.CELLULAR_COMPONENT, 200,
            SemanticSimilarityScore.SIMREL, true);
    System.out.println( "Cellular component, 200 samples, SimRel metric:\n"
            + Arrays.toString(percentiles) );
    percentiles = nullDistribution(GoNamespace.MOLECULAR_FUNCTION, 200,
            SemanticSimilarityScore.SIMREL, true);
    System.out.println( "Molecular function, 200 samples, SimRel metric:\n"
            + Arrays.toString(percentiles) );
    */
    
    
  }


  @Override
  public TermSimilarityMatrix_Backup clone() {

    TermSimilarityMatrix_Backup result = new TermSimilarityMatrix_Backup();

    for ( Integer aRowIdx : this.keySet() ) {

      HashMap<Integer,Double> curRow = this.get(aRowIdx);

      HashMap<Integer,Double> newRow = new HashMap<Integer,Double>( curRow.size() );
      result.put(aRowIdx, newRow);

      for ( int curColIdx : curRow.keySet() ) {
        newRow.put( curColIdx, new Double(curRow.get(curColIdx)) );
      }

    }

    result.supportedTerms = (GoTermList) this.supportedTerms.clone();
    return result;

  }


  /** This is not meant to be called from outside the class. */
  private TermSimilarityMatrix_Backup() {
    // do nothing
  }
  



  
}
