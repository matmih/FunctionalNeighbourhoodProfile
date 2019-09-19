package OntologyTools;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Fran Supek
 */
public class GoTermList extends ArrayList<GOTerm> {

  public GoTermList(Collection<? extends GOTerm> c) {
    super(c);
  }


  public GoTermList() {
    super();
  }


  public GoTermList(int initialCapacity) {
    super(initialCapacity);
  }


  private class GoTermSortByDispensability implements Comparator<GOTerm> {

    private final GoTermProperties termProperties;

    public GoTermSortByDispensability(GoTermProperties termProperties) {
      this.termProperties = termProperties;
    }

    @Override
    public int compare(GOTerm o1, GOTerm o2) {

      Double d1 = this.termProperties.get(o1, "dispensability");
      Double d2 = this.termProperties.get(o2, "dispensability");

      if ( d1 == null || d2 == null )
        throw new IllegalArgumentException("GOTerms in Comparator do not have" +
                " 'dispensability' property defined.");

      return d1.compareTo(d2);

    }

  }


  private class GoTermSortByRepresentative implements Comparator<GOTerm> {

    private final GoTermProperties termProperties;

    public GoTermSortByRepresentative(GoTermProperties termProperties) {
      this.termProperties = termProperties;
    }

    @Override
    public int compare(GOTerm o1, GOTerm o2) {

      Double repId1 = termProperties.get(o1, "representative");
      Double repId2 = termProperties.get(o2, "representative");

      if ( repId1 == null || repId2 == null )
        throw new IllegalArgumentException("GOTerms in Comparator do not have" +
                " 'representative' property defined.");

      // return 0 if equal, -1 if o1 is before o2, or +1 if o1 is after o2
      if ( repId1.equals(repId2) ) {

        // here, the representative is first among equals
        if ( repId1.intValue() == o1.getId() )
          return -1;
        else if ( repId2.intValue() == o2.getId() )
          return 1;
        else
          return 0;

      } else { // if terms have different representatives, the one with the less
               // dispensable representative goes first

        GOTerm o1_r = o1.getGeneOntology().get( repId1.intValue() );
        GOTerm o2_r = o2.getGeneOntology().get( repId2.intValue() );

        if ( termProperties.get(o1_r,"dispensability")
                < termProperties.get(o2_r,"dispensability") ) {
          return -1;
        } else if ( termProperties.get(o1_r, "dispensability")
                > termProperties.get(o2_r, "dispensability") ) {
          return 1;
        } else {    // if both representatives are equally dispensible, the one
                    // with lower GO number goes first
          return repId1.compareTo(repId2);
        }
        
      }

    }
    
  }




  public void sortByDispensability(GoTermProperties termProps) {
    Collections.sort(this, new GoTermSortByDispensability(termProps));
  }

  public void sortByRepresentative(GoTermProperties termProps) {
    Collections.sort(this, new GoTermSortByRepresentative(termProps));
  }


  public void findClustersAndSortByThem( GoTermProperties termProperties, double dispensabilityCutoff ) {

    // first make a set of "untouchables" - those GO terms that are below
    // dispensability cutoff
    Set<GOTerm> representatives = new HashSet<GOTerm>();
    for ( GOTerm term : this ) {
      if ( termProperties.get(term, "dispensability") < dispensabilityCutoff )
        representatives.add(term);
    }

    // now, reconstruct the whole chain, who gets dispensed by whom
    Map<GOTerm, GOTerm> dispensedBy = new HashMap<GOTerm, GOTerm>();
    for ( GOTerm term : this ) {
      dispensedBy.put(term, term.getGeneOntology().get( termProperties.get(term,"dispensedBy").intValue() )  );
    }

    // now determine the representatives using the last two arrays
    for ( GOTerm term : this ) {
      GOTerm representativeTerm = term;  // first, the term tries to represent itsels
      while ( ! representatives.contains(representativeTerm) ) {
        representativeTerm = dispensedBy.get(representativeTerm);
      }
      termProperties.put(term, "representative", (double) representativeTerm.getId() );
    }

    // finally, sort list by representatives
    sortByRepresentative(termProperties);

  }



}
