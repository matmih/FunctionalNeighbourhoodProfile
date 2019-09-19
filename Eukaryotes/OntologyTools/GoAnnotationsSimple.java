package OntologyTools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

/**
 * Loads the GO annotations from the EBI Uniprot-GOA project
 * http://www.ebi.ac.uk/GOA/proteomes
 * in a simple object, linking UniProt IDs to a set of GOs.
 * 
 * Does not store annotation type (IEA vs non-IEA).
 * 
 * @author Fran Supek
 */
public class GoAnnotationsSimple {
 
  
  /**
   * Stores the GO terms for each protein, keyed by the "long" UniProtID, such 
   * as "BRO1_YEAST". Necessarily includes all parent terms for a given term as 
   * propagateUpwards() is invoked at construction time.
   */
  public final Map<String, Set<GOTerm>> data = new LinkedHashMap<>();

  /**
   * Links eg P48582 to BRO1_YEAST, which can then be looked up in the "data"
   * as a key. When looking up gene function, if it found nothing using the "data"
   * tries to use this field to translate.
   */
  public final Map<String, String> shortToLongUniprotId = new HashMap<>();
  
  
  /** Stores a reference to the GeneOntology object supplied at creation time. */
  public final GeneOntology go;
  
  
  /**
   * Loads an annotation file. Can be gzipped. All annotations are automatically
   * propagated upwards.
   * 
   * @param inFile
   */
  public GoAnnotationsSimple(String inFile, GeneOntology go) throws IOException {
    
    BufferedReader in;
    if ( inFile.endsWith(".gz") ) {
      GZIPInputStream gzip = new GZIPInputStream( new FileInputStream(inFile) );
      in = new BufferedReader(new InputStreamReader(gzip));
    } else {
      in = new BufferedReader( new FileReader( inFile ) );
    }
    String line;
    
    while ( ( line = in.readLine() ) != null ) {
      
      if ( line.isEmpty() || line.startsWith("!") ) { 
        continue;
      }
      
      // cols[1] is the short UniProtID, e.g. P48582  
      // cols[4] is the GO ID
      // cols[10] contains identifiers separated by pipes, where the first one (in form "LCE6A_HUMAN") is interesting to us
      String[] cols = line.split("\t");
      int goId;
      try {
        goId = Integer.parseInt(cols[4].substring(3));
      } catch (NumberFormatException numberFormatException) {
        continue; // for some reason, this is unparseable
      }
      GOTerm term = go.get(goId);
      if ( term==null )  // could not find it, maybe the GO version is ancient?
        continue; 
      String proteinId = cols[10].split("\\|",2)[0]; // like "3HAO_YEAST"

      Set<GOTerm> annots;
      if ( data.containsKey( proteinId )  ) {
        annots = data.get( proteinId );
      } else {
        annots = new HashSet<>();
        data.put( proteinId, annots );
      }
      annots.add(term);  // note: this does not necessarily include all parent terms!
      shortToLongUniprotId.put( cols[1], proteinId );          
    }
    
    in.close();
    this.go = go;
    propagateUpwards();
  }
  
  /**
   * Loads a map of GO-annotated genes. All annotations are automatically
   * propagated upwards.
   * 
   * @param gos    A map containing GO-annotated genes
   * @param go     Gene Ontology
   * @author Vedrana Vidulin <vedrana.vidulin[AT]irb.hr>
   */
  public GoAnnotationsSimple(Map<String, Set<Integer>> gos, GeneOntology go)
  {
    for (String gene_ID : gos.keySet())
    {
        Set<Integer> goAnnots = gos.get(gene_ID);
        
        Set<GOTerm> annots = new HashSet<>();
        if (data.containsKey(gene_ID))
            annots = data.get(gene_ID);
        
        for (Integer g : goAnnots)
        {
            
            GOTerm term = go.get(g);
            
            if(g== 3674 && term==null)
                System.out.println("Contained function but null term!");
            else if(g == 3674 && term!=null){
                System.out.println(term.name);
            }
            
            if (term == null) //Could not find it, maybe the GO version is ancient?
                continue;
            
            annots.add(term); //Note: this does not necessarily include all parent terms!
        }
        
        if (annots.size() > 0)
            data.put(gene_ID, annots);
    }
    
    this.go = go;
    propagateUpwards();
    
    /*System.out.println("-----------------------------------------------------------------------------------------");
    for (String gene_ID : gos.keySet())
    {
        System.out.print("gene_ID: " + gene_ID + " - Original: " + gos.get(gene_ID) + " - With Parents: ");
        
        Set<GOTerm> annots = data.get(gene_ID);
        
        for (GOTerm g : annots)
            System.out.print(g.getFormattedId() + ", ");
        
        System.out.println();
    }
    System.out.println("-----------------------------------------------------------------------------------------");*/
    
    //System.out.println(data);
  }
  
  /**
   * For all annotations stored in field "data", makes sure that all of their 
   * GO parents are included as well. Is run automatically at construction time.
   */
  public void propagateUpwards() {
    for ( Set<GOTerm> annotsForGene : data.values() ) { // we don't really care about the gene here
      Set<GOTerm> newAnnots = new HashSet<>();
      for ( GOTerm anAnnot : annotsForGene ) {
        newAnnots.addAll( anAnnot.getAllParents() );
      }
      annotsForGene.addAll(newAnnots);
    }
  }
  
  
  /**
   * Counts how many annotations are there per GO term. Not cached, avoid
   * repeated invocations.
   * 
   * @return A Map linking the GOTerms to the integer count of annotations.
   * Divide by the data.size() to get the relative frequency. GO terms absent
   * from this Map have a count of 0.
   */
  public Map<GOTerm, Integer> countAnnots() {
    Map<GOTerm, Integer> result = new HashMap<>();
    for ( Set<GOTerm> annotsForGene : data.values() ) { // we don't really care about the gene here
      for ( GOTerm term : annotsForGene ) {
        if ( result.containsKey(term) ) {
          result.put( term, result.get(term) + 1 );
        } else {
          result.put( term, 1 );
        }
      }  // all annotations for current gene
    }  // gene by gene
    return result;    
  }
  
  
  /**
   * Filters out annotations that occur infrequently.
   * 
   * @param keepThreshold A GO term must be annotated to at least this many
   * genes to be retained (in any gene).
   * 
   * @return A Set of removed GO terms (they were removed from ALL genes).
   */
  public Set<GOTerm> removeRareAnnotations( int keepThreshold ) {
    Map<GOTerm, Integer> termCounts = this.countAnnots();
    Set<GOTerm> toRemove = new HashSet<>(); // will remove annotations with these
                                            // terms across ALL genes
    for ( GOTerm term : termCounts.keySet() ) {
      if ( termCounts.get(term) < keepThreshold ) {
        toRemove.add(term);
      }
    }
    removeAnnotsFromAllGenes(toRemove);
    return toRemove;
  }
  
  
  
  /**
   * Removes annotations to a set of GO terms from all genes. 
   * After removing, propagateUpwards() is invoked so some terms may be added
   * back (if their children were in the set).
   */ 
  public void removeAnnotsFromAllGenes( Set<GOTerm> termsToRemove ){
    for ( Set<GOTerm> annotsForGene : data.values() ) { 
      Iterator<GOTerm> it = annotsForGene.iterator();
      while (it.hasNext()) {
        if ( termsToRemove.contains( it.next() ) ) {
          it.remove();
        }
      }
    }  // gene by gene (or more precisely, their GO annotations)
    propagateUpwards();
  }
  
  
  
  public int getNumGenes() {
    return data.size();
  }
  
  public Set<GOTerm> getGeneAnnots(String geneId) {
    return data.get(geneId);
  }
  
  
  /**
   * Returns all annotations in a CLUS-HMC compatible format. <p>
   * 
   * The three root terms (GO:0008150 biological_process, 
   * GO:0003674 molecular_function and GO:0005575 cellular_component) are never
   * output in these annotations.
   * 
   * @param geneId By default, takes the long Uniprot ID, e.g. BRO1_YEAST,
   * but will also attempt to translate the short Uniprot IDs (like P48582)
   * by using the shortToLongUniprotId dictionary.
   * 
   * @return Something like "GO00059P4858275@GO0006066@GO0007047@GO0006412@GO0006464".
   * If no annotations, returns the string "root".
   */
  public String getGeneAnnotsAsClusString(String geneId ) {
    StringBuilder result = new StringBuilder();
    
    Set<GOTerm> annots = getGeneAnnots(geneId);
    if ( annots == null && this.shortToLongUniprotId.containsKey(geneId) ) {
      annots = getGeneAnnots( this.shortToLongUniprotId.get(geneId) ); 
    }
    if ( annots != null ) {
      for ( GOTerm annot : annots ) {
        if ( annot.isTopLevel() )
          continue;
        result.append( annot.getFormattedId().replaceAll(":","") ).append("@");
      }
    }
    if ( result.length() > 0 )
      result.deleteCharAt(result.length()-1); // remove terminal @
    else
      result.append("root");  // we found nothing, so we say "root" (who knows how CLUS will treat this...)
    return result.toString();
  }
  
  
  
  
  /**
   * Writes out the whole GO hierarchy in the format required for the header
   * of CLUS-HMC. <p>
   * 
   * Importantly, writes out only the terms actually present in the annotations.
   * It will thus not bother CLUS with e.g. terms with few annotations (if they 
   * were manually filtered before), or with namespaces that weren't loaded. <p>
   * 
   * The three root terms (GO:0008150 biological_process, 
   * GO:0003674 molecular_function and GO:0005575 cellular_component) are all
   * merged together and simply output as "root" here. Thus, if the annotations
   * contain different namespaces, they will be combined into one.
   * 
   * @param geneId
   * @return 
   */
  
  public HashSet<Integer> usedGOs=new HashSet<>();
  
  public String getGoHierarchyAsClusString() {
    StringBuilder output = new StringBuilder();
    
    // add all used terms - their parents were also automatically included at construction time. 
    
    // the list should be created top-to-bottom (so CLUS can reconstruct the hierarchy).
    Set<GOTerm> termsOfInterest = this.countAnnots().keySet();
    
    Set<String> alreadyPrinted = new HashSet<>();  // stores strings in form "root/GO0051234" or "GO0051179/GO0051234" ... so as not to print them twice (this will occur because of multi-parent relationships)
    
    List<GOTerm> lastExamined = Arrays.asList( go.get(8150), go.get(3674), go.get(5575) );
    while ( ! lastExamined.isEmpty() ) {
      List<GOTerm> nextGen = new ArrayList<>();
      for ( GOTerm parentTerm : lastExamined ) {
        for ( GOTerm term : parentTerm.getChildren() ) {
        
        
          if ( ! termsOfInterest.contains(term) ) {
            continue; //will not be printed, and its descendants won't be marked for examining at the next level
          }
          String toAppend;
          if ( parentTerm.isTopLevel() )
            toAppend = "root/" ;  // all three root terms are pooled together in a common "root"
          else
            toAppend = String.valueOf(parentTerm.getId()) + "/";
            //toAppend = parentTerm.getFormattedId().replaceFirst(":", "") + "/";
          System.out.println("tID: "+term.getId());
          System.out.println("pID: "+parentTerm.getId());
          toAppend += term.getId();
          usedGOs.add(term.getId());
          usedGOs.add(parentTerm.getId());
          //toAppend += term.getFormattedId().replaceFirst(":", "");

          if ( alreadyPrinted.contains(toAppend) ) {  // we already had this combination of parent + child
            continue; // will not be printed, and its descendants won't be marked for examining at the next level
          }
          output.append(toAppend).append(",");
          alreadyPrinted.add(toAppend);
          nextGen.add(term);
        
        }
      } // GO term by GO term (from last layer)
      
      lastExamined = nextGen;  // move down one level
      
    }  // level by level
    
    output.deleteCharAt(output.length()-1);  // remove the final comma
    return output.toString();
  }
  
  
  
  
  /**
   * A testing routine.
   */
  public static void main(String[] args) throws IOException
  {
        GeneOntology go = new GeneOntology("D:\\Programi\\GeneFunctionPrediction\\resources\\go_daily-termdb.obo-xml.gz");
        
        //Map<Integer, List<Integer>> genes = new TreeMap<Integer, List<Integer>>();
        
        Set<Integer> gos = new HashSet<Integer>();
        gos.add(47244);
        gos.add(9058);
        //genes.put(10000000, gos);
        
        //gos = new ArrayList<Integer>();
        gos.add(16740);
        //genes.put(10000001, gos);
        
        //gos = new ArrayList<Integer>();
        gos.add(4872);
        gos.add(5886);
        gos.add(6810);
        gos.add(9279);
        gos.add(5215);
        //genes.put(10000005, gos);
        
        //gos = new ArrayList<Integer>();
        gos.add(4872);
        gos.add(5886);
        gos.add(6810);
        gos.add(9279);
        gos.add(5215);
        //genes.put(10000008, gos);
        
        //System.out.println("Genes: " + genes);
        System.out.println("GO annotations: " + gos);
        
        //GoAnnotationsSimple annots = new GoAnnotationsSimple(gos, go);
        //System.out.println(annots.getGoHierarchyAsClusString());
      
    /*GeneOntology go = new GeneOntology(
            "C:\\Users\\fsupek\\Genomes\\CAFA2013\\go_201312-termdb.obo-xml.gz" 
    );
    GoAnnotationsSimple annots = new GoAnnotationsSimple(
            //"C:\\Users\\fsupek\\Genomes\\CAFA2013\\geneAnnots_2.1.2014\\gene_association.goa_human.gz",
            //"C:\\Users\\fsupek\\Genomes\\CAFA2013\\UniprotGOA_10Dec2013\\25.H_sapiens.goa",
            "C:\\Users\\fsupek\\Genomes\\CAFA2013\\UniprotGOA_10Dec2013\\71242.S_cerevisiae_ATCC_204508.goa",
            go);
    
    annots.removeRareAnnotations(1000);
    //System.out.println( annots.countAnnots().toString().replaceAll(", ",",\n") );
    System.out.println( annots.getGeneAnnots("BRAF_HUMAN") );
    System.out.println( annots.getGeneAnnots("P53_HUMAN") );
    System.out.println( annots.getGeneAnnotsAsClusString("BRAF_HUMAN") );
    System.out.println( annots.getGeneAnnotsAsClusString("P53_HUMAN") );
    System.out.println( annots.getGoHierarchyAsClusString() );
    //System.out.println( annots.data.toString() );*/
  }
  
  
  
}
