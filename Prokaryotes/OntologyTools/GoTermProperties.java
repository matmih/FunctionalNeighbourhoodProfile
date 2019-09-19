/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package OntologyTools;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Stores a Map of arbitraty String-Double key-value pairs for GO Terms.
 * Empty at creation.
 *
 * @author Fran
 */
public class GoTermProperties implements Serializable {


  private Map<GOTerm, HashMap<String, Double>> data
          = new HashMap<GOTerm, HashMap<String,Double>>();


  private GeneOntology myGo = null;

  /**
   * Please provide a GO object at initialization time, for convenience in
   * finding GO terms by ID later.
   */
  public GoTermProperties(GeneOntology go) {
    myGo = go;
  }


  public void clearAll() {
    data.clear();
  }

  public Double get( GOTerm term, String propName ){
    return data.get(term).get(propName);
  }

  public Double get( int termId, String propName ){
    return data.get( myGo.get(termId) ).get(propName);
  }

  public boolean hasTermWithProperty( GOTerm term, String propName ) {
    if ( ! data.containsKey(term) )
      return false;
    else
      return data.get(term).containsKey(propName);
  }

  public boolean hasTermWithProperty( int termId, String propName ) {
    return this.hasTermWithProperty( myGo.get(termId), propName );
  }


  
  public boolean hasTerm( GOTerm term ){
    if ( ! data.containsKey(term) )
      return false;
    else
      return data.get(term) != null;
  }


  public boolean hasTerm( int termId ){
    return this.hasTerm( myGo.get(termId) );
  }

  
  public void put( GOTerm term, String propName, double value ) {

    if ( data.containsKey(term) ) {
      data.get(term).put(propName, value);
    } else {
      HashMap<String, Double> newMap = new HashMap<String, Double>();
      newMap.put(propName, value);
      data.put(term, newMap);
    }

  }

  public void put( int termId, String propName, double value ) {

    GOTerm term = myGo.get(termId);
    this.put( term, propName, value );

  }


  public Map<String,Double> getAllProps( GOTerm term ) {
    return data.get(term);
  }

  public Map<String,Double> getAllProps( int termId ) {
    return data.get(myGo.get(termId));
  }


}
