package OntologyTools;

import java.util.Collection;
import java.util.HashMap;
import java.util.Random;

/**
 * A directed acyclic graph. Essentially a Map containing Nodes, indexed by 
 * an Integer ID.
 * 
 * @author Fran Supek
 */
public class Graph<V extends Node> extends HashMap<Integer,V> {

  //private static final long serialVersionUID = 1L;
  
  /** Returns a random Node. */
  public V getRandomNode() {
    Random r = new Random();
    return this.get( this.keySet().toArray(new Integer[0])[r.nextInt( this.size() )] );
  }
  
  /** Empties the Properties Maps of each Node. */
  public void clearAllProperties() {
    for ( Node curNode : this.values() )
      curNode.properties.clear();
  }
  
  public Collection<V> getAllNodes() {
    return this.values();
  }

  @Override
  public V remove(Object key) {
    return super.remove(key);
  }
  
  
  
  
}
