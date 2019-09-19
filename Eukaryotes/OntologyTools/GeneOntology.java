package OntologyTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

  
/**
 * A class representing the whole Gene Ontology.
 * 
 * @author Fran Supek (fran.supek[AT]irb.hr)
 */
public class GeneOntology extends OntologyTools.Graph<GOTerm> {

  /** For serialization. */
  //static final long serialVersionUID = 1;

  /** Used when no Gene Ontology OBO-XML file is specified in the constructor. */
  private static final String defaultOboXmlFile =
          "go_201404-termdb.obo-xml.gz";
          //"go_20141215-termdb.obo-xml.gz";

  /**
   * If the Gene Ontology OBO-XML file is not in the working directory, look
   * for it here.
   */
  public static final String defaultOboXmlDirectory =
            "Z:/databases/GeneOntology/";
          //"/home/vedrana/genomics/databases/GeneOntology/";
  
  /** As given in the header of the OBO-XML file. */
  private String date;

  public String getDate() {
    return date;
  }



  /**
   * Reads in the structure of the whole Gene Ontology from an OBO-XML
   * file.<p>
   * 
   * The files are downloadable from http://archive.geneontology.org/latest/
   * and are named "go_??????-termdb.obo-xml.gz".
   *
   * @throws java.io.IOException
   */
  public GeneOntology(String fileName) throws IOException {
      
      this( fileName.toLowerCase().endsWith(".gz")
            ?
              new GZIPInputStream(new FileInputStream( fileName ) )
            :
              new FileInputStream(fileName) );
    
  }

  
  /**
   * Reads in the structure of the whole Gene Ontology from the default OBO-XML
   * file.
   *
   * @throws IOException
   */
  public GeneOntology() throws IOException {
    this( defaultOboXmlDirectory + defaultOboXmlFile );
  }

  
  /**
   * Reads in the structure of the whole Gene Ontology from
   * an InputStream. <p>
   *
   * Covenient for use with getServletContext().getResourceAsStream("/yourfilename.txt")
   *
   * The files are downloadable from http://archive.geneontology.org/latest/
   * and are named "go_??????-termdb.obo-xml.gz".
   *
   * @throws java.io.IOException
   */ 
  public GeneOntology ( InputStream stream ) throws IOException {

    this.date = "date of OBO-XML file is unknown";

    BufferedReader in = new BufferedReader( new InputStreamReader(stream) );
    in.mark( 80000000 ); // 80 megabytes.... ouch
    
    String line;
    
    int curGoId = -1, goAltId = -1, parent = -1;
   
    // first pass - read in all GO terms and their names, ids, and namespace
    GOTerm curGoTerm = null;
    while ( ( line = in.readLine() ) != null ) {
      
      line = line.trim();
      String lineLow = line.toLowerCase();
      
      if ( lineLow.equals("<typedef>") ) break;
      
      if ( lineLow.startsWith("<id>go:") ) {
        curGoId = Integer.parseInt( line.substring(7, 14) );
        curGoTerm = new GOTerm(curGoId);
        this.put(curGoId, curGoTerm);
        curGoTerm.geneOntology = this;
         
      }
      
      if ( lineLow.startsWith("<name>") ) {
        if ( curGoTerm == null )
          continue;
        if (curGoTerm.getName() == null)
          curGoTerm.setName( line.substring(6, line.length()-7) );
      }
      
      if ( lineLow.startsWith("<namespace>") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.namespace = GoNamespace.valueOf(
                line.substring(11, line.length()-12 ).toUpperCase() );
      }
      
      if ( lineLow.startsWith("<alt_id>go:") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        goAltId = Integer.parseInt( line.substring(11, 18) );
        curGoTerm.getAltIds().add(goAltId);
        this.put(goAltId, curGoTerm);
      }
      
      if ( lineLow.startsWith("<subset>gosubset_prok</subset>")  ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.inProkaryoticSubset = true;          
      }
      
      if ( lineLow.startsWith("<subset>goslim_generic</subset>") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.inGoSlim = true;          
      }
      
      if ( lineLow.startsWith("<is_obsolete>1") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.obsolete = true;          
      }

      // alternatives for obsolete terms?
      if ( lineLow.startsWith("<replaced_by>") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        String replacementTermString =
                line.replace("<replaced_by>", "").replace("</replaced_by>", "");
        int replacementId = Integer.parseInt( replacementTermString.substring(3) );
        if (curGoTerm.replaceByTermId == 0)
          curGoTerm.replaceByTermId = replacementId;
      }


      if ( lineLow.startsWith("<defstr>") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.setDescription(line.replace("<defstr>", "").replace("</defstr>", ""));
      }

      if ( lineLow.startsWith("<synonym_text>") ) {
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        curGoTerm.getAltNames().add(line.replace("<synonym_text>", "").replace("</synonym_text>", ""));
      }
      
      
      
    }
    in.reset();

    // second pass - read in parents for each term
    curGoTerm = null;

    while ( ( line = in.readLine() ) != null ) {
      
      line = line.trim().toLowerCase();
      
      if ( line.equals("<typedef>") ) break;
      
      if ( line.startsWith("<id>go:") ) {
        curGoId = Integer.parseInt( line.substring(7, 14) );
        curGoTerm = this.get(curGoId);
        parent = curGoId;
      }      
      
      if ( line.startsWith("<is_a>go:")  // is-a relationship 
              || line.startsWith("<to>go:") ) {  // part-of relationship
        
        if ( curGoTerm == null )
          throw new IllegalArgumentException("OBO-XML file has weird contents!");
        
        int where = line.indexOf("go:") + 3;
        curGoId = Integer.parseInt( line.substring(where, where+7) );
        
        curGoTerm.getParents().add( this.get(curGoId) );
        // adding children
        this.get(curGoId).getChildren().add(this.get(parent));
      }
    }

    // third pass - just read in the date of the OBO-XML file
    in.reset();

    while ( ( line = in.readLine() ) != null ) {
      line = line.trim().toLowerCase();
      if ( line.startsWith("<date>") ) {
        this.date = line.replace("<date>", "").replace("</date>", "");
        break;
      }
      if ( line.equals("</header>") ) break;
    }


    in.close();

    // now, kill obsolete terms
    Set<Integer> allIds = new HashSet<Integer>();

    for ( int goId : this.keySet() ){
      GOTerm curTerm = this.get(goId);
      if ( curTerm.obsolete )
        allIds.add(goId);
    }
    
    for ( int goId : allIds ) {
      GOTerm curTerm = this.get(goId);
      if ( curTerm.obsolete && curTerm.replaceByTermId != 0 ) {
        this.put( goId, this.get(curTerm.replaceByTermId) );
      }

    }

  }
  

  
  /** Returns the name of a specific GO term. */
  public String translateId( int id ){
    if ( this.containsKey(id) )
      return this.get(id).getName();
    else
      return "Unknown GO category";
  }

  
  
  /**
   * A utility function that takes a file with GO IDs (one per line) and splits
   * t into three files, one for each of the three gene ontologies.
   * 
   * @param filename
   */
  public void splitGOtermFile ( String filename ) throws Exception {
    BufferedReader in = new BufferedReader( new FileReader( filename ) );
    String line;
    
    PrintWriter molFunPw = new PrintWriter( new FileWriter(
            "molecularFunction_" + filename, false ), true );
    PrintWriter biolProcPw = new PrintWriter( new FileWriter(
            "biologicalProcess_" + filename, false ), true );
    PrintWriter celCompPw = new PrintWriter( new FileWriter(
            "cellularComponent_" + filename, false ), true );
    
    while ( ( line = in.readLine() ) != null ) {
      int curGoId = Integer.parseInt( line.replace("GO:", "") );
      int curTopNodeId = get(curGoId).getTopNode().getId();
      PrintWriter aPW;
      
      if ( curTopNodeId == 3674  )       // mol. function
        aPW = molFunPw;
      else if ( curTopNodeId == 5575  )  // component
        aPW = celCompPw;
      else                               // biol. process
        aPW = biolProcPw;
        
      aPW.printf("GO:%07d\n", curGoId);
      
    }
    in.close();
    molFunPw.close(); celCompPw.close(); biolProcPw.close();
    
  }

  
  /**
   * Returns a List with all GO terms inside. Can be limited to GO Slim, or
   * prokaryotic subset only. If not limited, the resulting list will be 
   * quite large. Never outputs obsolete GO categories.
   * 
   * @param goSlimOnly
   * @param prokaryoticSubsetOnly
   * @return
   */
  public Collection<GOTerm> getAllTerms( boolean goSlimOnly,
          boolean prokaryoticSubsetOnly ) {

    Set<GOTerm> result = new LinkedHashSet<GOTerm>();
    
    for ( int curGoTermId : this.keySet() ) {
      GOTerm curGOTerm = this.get(curGoTermId);
      if ( goSlimOnly && !curGOTerm.isGoSlim() )
        continue;
      if ( prokaryoticSubsetOnly && !curGOTerm.isInProkaryoticSubset() )
        continue;
      if ( curGOTerm.isObsolete() )
        continue;
      result.add(curGOTerm);
    }
    return result;
    
  }
  

    public static GeneOntology deserializeFromFile ( File file )
          throws IOException, ClassNotFoundException {

    InputStream is = new FileInputStream(file);
    ObjectInput oi = new ObjectInputStream(is);
    Object newObj = oi.readObject();
    if ( ! (newObj instanceof GeneOntology) )
      throw new IllegalArgumentException("Wrong class in file;");
    oi.close();
    return (GeneOntology)newObj;

  }


  public void serializeToFile( File file )
          throws FileNotFoundException, IOException {

    OutputStream os = new FileOutputStream(file);
    ObjectOutput oo = new ObjectOutputStream(os);
    oo.writeObject(this);
    oo.close();

  }

  /**
   * Clears the "properties" HashMap of all constituent GO Terms.
   */
  public void clearGoTermProperties() {
    for ( GOTerm goTerm : this.values() ) {
      goTerm.getProperties().clear();
    }

  }


  public void addBasicKeywordsToAllTerms() {
    for (GOTerm term : this.values()) {
      term.getKeywords();
    }
  }


  public void addKeywordsFromUniprotKeywords( InputStream keywlistFileStream ) throws IOException {

    BufferedReader in = new BufferedReader( new InputStreamReader(keywlistFileStream) );

    String line;

    StringBuilder curText = new StringBuilder();
    while ( ( line = in.readLine() ) != null ) {
      
      if ( line.startsWith("DE") ){
        curText.append(line.replaceFirst("^DE  ", ""));
        if ( !line.endsWith("-") )
          curText.append(" ");
        continue;
      }

      if ( line.startsWith("ID") ) {
        curText.append(line.replaceFirst("^ID  ", ""));
        continue;
      }

      if ( line.startsWith("SY") ) {
        curText.append(line.replaceFirst("^SY  ", ""));
        continue;
      }

      if ( line.startsWith("//") ) {
        curText = new StringBuilder();
        continue;
      }

      if ( line.startsWith("GO") ) {
        int termId = Integer.parseInt(line.substring(8, 15));

        Set<String> newKeywords = new HashSet<String>();
        String oneBigString = curText.toString().toLowerCase();
        StringTokenizer tok = new StringTokenizer(oneBigString, ",:=;. ");
        while (tok.hasMoreTokens()) {
          String nextTok = tok.nextToken();
          if ( nextTok.startsWith("(") )
            nextTok = nextTok.substring(1);
          if ( nextTok.endsWith(")") )
            nextTok = nextTok.substring(0, nextTok.length()-1);
          // the word "protein" is very overrepresented in Uniprot definitions
          if ( !nextTok.equals("protein") )
            newKeywords.add(nextTok); // a Set does not allow duplicates
        }

        this.get(termId).keywords.addAll(newKeywords);


        curText = new StringBuilder();
        continue;
      }

    }

    in.close();


  }



  
  /**
   * A demonstration function for the GeneOntology class.
   * 
   * @param args
   * @throws java.lang.Exception
   */
  public static void main(String[] args) throws Exception {
   
    // make sure that OBO-XML file is in the right directory before testing...
    
    Long millis;
    GeneOntology myGO;

    millis = System.currentTimeMillis();
    myGO= new GeneOntology("Z:/databases/GeneOntology/go_200911-termdb.obo-xml");
    System.out.println( System.currentTimeMillis() - millis );

    /*millis = System.currentTimeMillis();
    myGO= new GeneOntology("Z:/databases/GeneOntology/go_200911-termdb.obo-xml");
    System.out.println( System.currentTimeMillis() - millis );

    millis = System.currentTimeMillis();
    myGO = GeneOntology.deserializeFromFile(new File("C:/Users/Fran/Genomes/Revigo/RevigoWepApp/web/data/go_200911-termdb.ser"));
    System.out.println( System.currentTimeMillis() - millis );

    millis = System.currentTimeMillis();
    myGO = GeneOntology.deserializeFromFile(new File("C:/Users/Fran/Genomes/Revigo/RevigoWepApp/web/data/go_200911-termdb.ser"));
    System.out.println( System.currentTimeMillis() - millis );*/

    //myGO.serializeToFile(new File("C:/Users/Fran/Genomes/Revigo/RevigoWepApp/web/data/go_200911-termdb.ser"));

    // (1) translate some GO IDs to see if this works...
    
    System.out.println( "GO ID  6662 is: " + myGO.translateId( 6662 ));
    System.out.println( "GO ID 44238 is: " + myGO.translateId( 44238 ));    
    System.out.println( "GO ID    27 is: " + myGO.translateId( 27 ));
    System.out.println( "GO ID 60245 is: " + myGO.translateId( 60245 ));

    // (2) output all the IDs to a file
    /*PrintWriter aPW = new PrintWriter( new FileWriter(
            "all_GO_IDs_200710-termdb.obo-xml.txt", false ), true );
    
    GOTerm[] allGoTerms = new GOTerm[0];
    for ( GOTerm curGOterm : myGO.values() ) {
      aPW.printf( "GO:%07d/n", curGOterm.getId() );
    }
    
    aPW.close(); */
    
    // (3) try the splitGOtermFile function
    // myGO.splitGOtermFile( "all_GO_IDs_200710-termdb.obo-xml.txt" );
    
    // (4) try the getAllGoTerms function
    /*Collection<GOTerm> subset = myGO.getAllTerms(false, false);
    for ( GOTerm curGoTerm : subset ) {
      System.out.println(curGoTerm + ":" + curGoTerm.getKeywords().toString());
    }*/

    Collection<GOTerm> parents = myGO.get(7067).getAllParents();
    for ( GOTerm curPar : parents ) {
      System.out.println( curPar );

    }



  }


}
