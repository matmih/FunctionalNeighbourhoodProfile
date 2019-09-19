package OntologyTools;

import java.io.Serializable;
import java.util.List;

/**
 *
 * @author Fran Supek
 */
public enum GoNamespace implements Serializable {

  BIOLOGICAL_PROCESS, CELLULAR_COMPONENT, MOLECULAR_FUNCTION, MIXED_NAMESPACE;

  /**
   * Returns if this Gene Ontology Namespace is the parent of the supplied 
   * GOTerm.
   */
  public boolean isParentOf(GOTerm term) {

    switch (this) {
      case MOLECULAR_FUNCTION:
        return term.isMolecularFunction();
      case CELLULAR_COMPONENT:
        return term.isCellularComponent();
      case BIOLOGICAL_PROCESS:
        return term.isBiologicalProcess();
      case MIXED_NAMESPACE:
        return true;
      default:
        throw new IllegalArgumentException("Namespace not supported.");
    }

  }

  @Override
  public String toString() {
    switch (this) {
      case MOLECULAR_FUNCTION:
        return "Molecular Function";
      case CELLULAR_COMPONENT:
        return "Cellular Component";
      case BIOLOGICAL_PROCESS:
        return "Biological Process";
      case MIXED_NAMESPACE:
        return "Mixed Namespace";
      default:
        throw new IllegalArgumentException("Namespace not supported.");
    }
  }


  /**
   * Determine the namespace of supplied terms, set to null if list has mixed
   * namespaces, or is empty.
   */
  public static GoNamespace findNamespaceOfList(List<GOTerm> termList) {

    GoNamespace namespaceOfFirstTerm;
    if (termList.size() == 0)
      return null;
    else
      namespaceOfFirstTerm = termList.get(0).getNamespace();

    for (GOTerm curGoTerm : termList) {
      if (curGoTerm.getNamespace() != namespaceOfFirstTerm)
        return GoNamespace.MIXED_NAMESPACE;
    }
    return namespaceOfFirstTerm;

  }


}
