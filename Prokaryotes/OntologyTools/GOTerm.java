package OntologyTools;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

/**
 * Class representing a single Gene Ontology term.
 * 
 * @author Fran Supek (fran.supek[AT]irb.hr)
 */
public class GOTerm extends hr.irb.directedAcyclicGraph.Node {

  /**
   * Is a term part of the Gene Ontology prokaryotic subset? Default = false.
   * see http://www.geneontology.org/GO.slims.shtml#prok
   */ 
  protected boolean inProkaryoticSubset = false;
  /**
   * Is a term part of the GOslim-generic subset? Default = false.
   * see http://www.geneontology.org/GO.slims.shtml#prok
   */ 
  protected boolean inGoSlim;
  
  /**
   * Is term declared obsolete? Default = false;
   */ 
  protected boolean obsolete = false;
  /**
   * For obsolete terms, here is a suitable replacement. (Note: just one 
   * replacement term is stored here while the OBO-XML file may offer more 
   * than one.
   */
  protected int replaceByTermId = 0;

  public boolean isGoSlim() {
    return this.inGoSlim;
  }

  public boolean isInProkaryoticSubset() {
    return this.inProkaryoticSubset;
  }
  
  public boolean isObsolete() {
    return this.obsolete;
  }
  
  
  /** The namespace this term belongs to. */
  protected GoNamespace namespace;

  public GoNamespace getNamespace() {
    return namespace;
  }


  /** The GeneOntology object this term belongs to. */
  protected GeneOntology geneOntology = null;

  public GeneOntology getGeneOntology() {
    return geneOntology;
  }


  /**
   * Retrieves a shorter version of the node's name, where: <ul>
   * <li> "metabolic process" is replaced by "metabolism"
   * <li> "catabolic process" is replaced by "catabolism"
   * <li> "biosynthetic process" is replaced by "biosynthesis"
   * </ul>
   * <p>
   * These replacements concern only terms from the Biological Function
   * namespace. They produce names that are legal synonyms for the original
   * names.
   *
   * @return
   */
  public String getAbbreviatedName() {

    if (this.getNamespace() == GoNamespace.BIOLOGICAL_PROCESS) {
      String newName = new String(this.getName());
      newName = newName.replace("metabolic process", "metabolism")
              .replace("catabolic process", "catabolism")
              .replace("biosynthetic process", "biosynthesis");
      return newName;
    } else
      return this.getName();
    
  }


  /** Returns GO Term ID in the format "GO:0006915". */
  public String getFormattedId() {
    return String.format("GO:%07d", this.getId());
  }
  
  public String getFormattedId2() {
    return String.format("GO%07d", this.getId());
  }


  protected Set<String> keywords = null;


  /**
   * Provides a set of all words (in lowercase) used in any of the term names,
   * or in the term's definition.
   */
  public Set<String> getKeywords() {

    if (this.keywords == null)
      this.keywords = makeKeywordsFromDescription();
    return this.keywords;

  }


  public Set<String> makeKeywordsFromDescription() {

    if (this.keywords != null)
      return this.keywords;

    StringBuilder sb = new StringBuilder(this.getName());
    for (String curAltName : this.getAltNames()) {
      sb.append(" ").append( curAltName );
     }
    sb.append(" ").append(this.getDescription());

    Set<String> result = new HashSet<String>();
    String oneBigString = sb.toString().toLowerCase();
    StringTokenizer tok = new StringTokenizer(oneBigString, ",:=;. ");
    while (tok.hasMoreTokens()) {
      String nextTok = tok.nextToken();
      if ( nextTok.startsWith("(") )
        nextTok = nextTok.substring(1);
      if ( nextTok.endsWith(")") )
        nextTok = nextTok.substring(0, nextTok.length()-1);
      result.add(nextTok); // a Set does not allow duplicates
    }

    this.keywords = result;
    return result;

  }


  
  /**
   * A unique ID (integer) must be specified for the GOTerm at the time of
   * construction.
   * 
   * @param IDparam The main unique ID for the node.
   */
  protected GOTerm (int idParam) {
    super(idParam);
  }    

  
  /** Returns true if the GOTerm belongs to the "Molecular Function" ontology */
  public boolean isMolecularFunction() {
    return (this.namespace == GoNamespace.MOLECULAR_FUNCTION);
  }
  
  /** Returns true if the GOTerm belongs to the "Biological Process" ontology */
  public boolean isBiologicalProcess() {
    return (this.namespace == GoNamespace.BIOLOGICAL_PROCESS);
  }

  /** Returns true if the GOTerm belongs to the "Cellular Component" ontology */
  public boolean isCellularComponent() {
    return (this.namespace == GoNamespace.CELLULAR_COMPONENT);
  }

  @Override
  public Set<GOTerm> getParents() {
    return (Set<GOTerm>) super.getParents();
  }


  /**
   * Returns all parents of the GO Term, all of their parents and so on,
   * searching recursivelly through the graph. Output does not
   * contain duplicate nodes.<p/>
   *
   * The list DOES NOT include the Term itself, only its parents!
   *
   * @return ArrayList of Nodes.
   */
  @Override
  public Set<GOTerm> getAllParents() {
    return (Set<GOTerm>) super.getAllParents();
  }


  /**
   * Returns all of the parent terms that this GOTerm has in common with
   * another GOTerm.
   */
  public Set<GOTerm> getAllCommonParents( GOTerm anotherTerm ) {
    
    return (Set<GOTerm>) super.getAllCommonParents(anotherTerm);

  }


  @Override
  public Set<GOTerm> getChildren() {
    return (Set<GOTerm>) super.getChildren();
  }


  /**
   * Returns all children of the GO Term, all of their children and so on,
   * searching recursivelly through the graph. Output does not
   * contain duplicate nodes.<p/>
   *
   * The list DOES NOT include the Term itself, only its children!
   *
   * @return ArrayList of Nodes.
   */
  @Override
  public Set<GOTerm> getAllChildren() {
    return (Set<GOTerm>) super.getAllChildren();
  }

  /**
   * Returns all the children of any parent of the given Term. When looking
   * for children, looks only one level deep. The returned collection does not
   * include the given term.
   */
  @Override
  public Set<GOTerm> getSiblings() {
    return (Set<GOTerm>) super.getSiblings();
  }


}
