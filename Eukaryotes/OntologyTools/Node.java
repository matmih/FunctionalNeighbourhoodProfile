package OntologyTools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A node in a directed acyclic graph.
 * 
 * The node must have at least one unique ID (of type int) specified at time
 * of creation, but may also have many alternate unique IDs.
 * 
 * The node may have a name (of type String), and it may have a number of 
 * alternate names.
 * 
 * The node may have zero, one or many parent nodes.
 * 
 * @author Fran Supek (fran.supek[AT]irb.hr)
 */
public class Node implements Serializable {

  String description = null;

  public String getDescription() {
    return description;
  }


  public void setDescription(String description) {
    this.description = description;
  }


  boolean traversable;

  public void setTraversable() {
    this.traversable = true;
  }


  public boolean getTraversable() {
    return traversable;
  }


  /**
   * A unique ID (integer) must be specified for the node at the time of
   * construction.
   *
   * @param IDparam The main unique ID for the node.
   */
  public Node(int idParam) {
    super();
    id = idParam;
  }


  /** The main unique ID for the node. */
  protected int id;

  /** Returns the main unique ID of the node. */
  public int getId() {
    return id;
  }


  /** List of alternate IDs for the node. */
  protected List<Integer> altIds = new ArrayList<Integer>();

  /** Returns a list of alternate IDs for the node. */
  public List<Integer> getAltIds() {
    return altIds;
  }


  /** The main name of the node. */
  protected String name = null;

  /** Sets the name of the node. */
  public void setName(String name) {
    this.name = name;
  }


  /** Returns the name of the node. */
  public String getName() {
    return name;
  }



  /** A list of all parents of the node. */
  protected Set<Node> parents = new HashSet<Node>();

  /**
   * Returns a reference to a list of all parents of the node (modifications
   * will be reflected in the original Node object).
   */
  public Set<? extends Node> getParents() {
    return parents;
  }

  /**
   * Creates a new Set containing references to all parents of the Node,
   * all of their parents and so on, searching recursively through the graph.
   *
   * @return Set of Nodes.
   */
  public Set<? extends Node> getAllParents() {
    Set<Node> myParentList = new HashSet<Node>();
    getAllParentsHelper(myParentList);
    return myParentList;
  }

  /** A helper function for getAllParents() recursive search. */
  private final void getAllParentsHelper(Set<Node> setParam) {
    for (Node curNode : parents) {
      setParam.add(curNode);
      curNode.getAllParentsHelper(setParam);
    }
  }
  
  
  /**
   * Creates a new Set containing references to all parents of the Node,
   * all of their parents and so on, searching recursively through the graph.
   * For each parent, also returns the distance to the current node.<p>
   * 
   * If null is supplied as "usePropAsWeight", the default behavior is to 
   * find the topological distance, i.e. every edge counts as 1. Otherwise, 
   * the value of the property for each node is used as the weight of the edge.
   * If the property is not defined for some of the nodes, throws an error.<p>
   *
   * @param inclusive The list of a node's parents includes itself, with a 
   * distance of 0.
   * 
   * @return Set of Nodes.
   */
  public Map<? extends Node, Double> getAllParentsWithDist( String usePropAsWeight, boolean inclusive ) {
    Map<Node, Double> myParentList = new HashMap<Node, Double>();
    getAllParentsWithDistHelper(myParentList, 0, 0.0, usePropAsWeight);
    if ( inclusive )
      myParentList.put(this, 0.0);
    return myParentList;
  }

  /** A helper function for getAllParents() recursive search. */
  private void getAllParentsWithDistHelper(Map<Node, Double> nodesAndDists, 
          int curTopoDist, double curPropDist, String usePropAsWeight) {
    curTopoDist++;
    if ( usePropAsWeight != null ) {
      if ( this.getProperties().get(usePropAsWeight) == null )
        throw new IllegalArgumentException("Node " + this.toString() + 
                " does not have the property " + usePropAsWeight + " defined." );
      curPropDist += this.getProperties().get(usePropAsWeight);
    }
    for (Node curParent : this.parents) {  // there may be more than one!
      if ( usePropAsWeight == null )
        nodesAndDists.put(curParent, (double)curTopoDist );
      else
        nodesAndDists.put(curParent, curPropDist);
      curParent.getAllParentsWithDistHelper(nodesAndDists, curTopoDist, curPropDist, usePropAsWeight);
    }
  }  
  

  
  /**
   * Returns a new set with references to all of the parent Nodes that this Node
   * has in common with another Node.
   */
  public Set<? extends Node> getAllCommonParents( Node anotherTerm ) {

    Set<Node> result = new HashSet<Node>();
    Set<Node> temp = new HashSet<Node>(this.getAllParents());

    for ( Node term : anotherTerm.getAllParents() ) {
      if ( temp.contains(term) )
        result.add(term);
    }

    return result;

  }  

  
  /**
   * Returns a new set with references to all of the parent Nodes that this Node
   * has in common with another Node.
   * 
   * In the special case when one term is a parent of another, returns the 
   * former term.
   */
  public Set<? extends Node> getAllCommonParentsInclusive( Node anotherTerm ) {

    Set<Node> result = new HashSet<Node>();
    Set<Node> temp = new HashSet<Node>(this.getAllParents());
    temp.add(this);
    
    Set<Node> temp2 = new HashSet<Node>(anotherTerm.getAllParents());
    temp2.add(anotherTerm);

    for ( Node term : temp2 ) {
      if ( temp.contains(term) )
        result.add(term);
    }

    return result;

  }  
  
  
  
  /**
   * Returns the deepest of the parent terms that this Node has in common with
   * another Node.<p>
   * 
   * @param inclusive If "true", in the special case when one parent node is a 
   * parent of another, the former node is returned. If "false", their parent
   * is sought.
   */
  public Node getFirstCommonParent( Node anotherNode, boolean inclusive ) {

    Set<Node> result = new HashSet<Node>();
    Set<Node> allCommonParents;
    if (inclusive ){
      allCommonParents = (Set<Node>) this.getAllCommonParentsInclusive(anotherNode);
    } else {
      allCommonParents = (Set<Node>) this.getAllCommonParents(anotherNode);
    }

    // now, the first (=deepest) common node is the one in the Set above that
    // is not a parent to any other node in the same Set (=has no children in
    // the set)
    Set<Node> allParentsOfAllCommonParents = new HashSet<Node>();
    for ( Node n : allCommonParents ) {
      allParentsOfAllCommonParents.addAll( n.getAllParents() );
    }
    for ( Node n : allCommonParents ) {
      if ( !allParentsOfAllCommonParents.contains(n) )
        result.add(n);
    }
    if ( result.isEmpty() ) {
      throw new IllegalArgumentException("No deepest common parents for " + this.toString()
              + " and " + anotherNode.toString() + " Very strange." );
    } else if ( result.size()>=2 ) {
      throw new IllegalArgumentException("More than one deepest common parent for" + this.toString()
              + " and " + anotherNode.toString() + " Very strange." );
    }
    
    return result.iterator().next();

  }  
  
    
  /**
   * Returns the deepest of the parent terms that this Node has. Can return 
   * the Node itself, if it has no parents.<p>
   */
  public Node getDeepestParent() {

    Set<Node> result = new HashSet<Node>();
    Set<Node> allParents;
    allParents = (Set<Node>) this.getAllParents();
    if ( allParents.isEmpty() )
      allParents.add(this); // special case
    for ( Node n : allParents ) {
      if ( n.parents.isEmpty() )
        result.add(n);
    }
    
    if ( result.isEmpty() ) {
      throw new IllegalArgumentException("No deepest parent for " + this.toString()
              + " Very strange." );
    } else if ( result.size()>=2 ) {
      throw new IllegalArgumentException("More than one deepest parent for " + this.toString()
              + " Very strange." );
    }
    
    return result.iterator().next();

  }    
  
  
  
  /**
   * Find distance to another node.<p>
   * 
   * "calcDiff" may be true, false or null, meaning the function will return:
   * the difference of the nodes' distances to the common ancestor (if true),
   * the sum (if false), or only the distance of the current node to the
   * common ancestor (if null).<p>
   * 
   * By default, if you supply null to "usePropAsWeight", the topological
   * distance will be calculated i.e. every edge counts as 1. If you supply 
   * a name of a property instead, this property will be used to find the
   * distance between all adjacent nodes. (If "usePropAsWeight" is given and 
   * some of the nodes don't have this property defined, an Exception will be 
   * thrown).
   */
  protected double topologicalDistanceTo( Node anotherNode, Boolean calcDiff, String usePropAsWeight ) {
    Node firstCommonParent = this.getFirstCommonParent(anotherNode, true);
    double dist1 = this.topologicalDistanceToParent(firstCommonParent, usePropAsWeight);
    if ( calcDiff==null )
      return dist1;
    else {
      double dist2 = anotherNode.topologicalDistanceToParent(firstCommonParent, usePropAsWeight);
      if ( calcDiff )
        return dist1 - dist2;
      else
        return dist1 + dist2;
    }
  }
  
  /**
   * A helper function for topologicalDistanceTo().
   */
  private double topologicalDistanceToParent( Node parentNode, String usePropAsWeight ) {

    // this is VERY inefficient! why "getAllParents" with dist, when we could have only a single one?
    Map<Node, Double> parentsWithDists = (Map<Node, Double>) this.getAllParentsWithDist( usePropAsWeight, true );
    return parentsWithDists.get(parentNode);

  }

  
  
  /**
   * Find distance to the root (actually, it's the distance to the deepest 
   * parent of the node).<p>
   * 
   * By default, if you supply null to "usePropAsWeight", the topological
   * distance will be calculated i.e. every edge counts as 1. If you supply 
   * a name of a property instead, this property will be used to find the
   * distance between all adjacent nodes. (If "usePropAsWeight" is given and 
   * some of the nodes don't have this property defined, an Exception will be 
   * thrown).
   * 
   * @return 
   */
  protected double topologicalDistanceToRoot( String usePropAsWeight ) {
    return topologicalDistanceToParent( this.getDeepestParent(), usePropAsWeight );
  }  
  
  
  
  
  /** A list of all children of the node. */
  protected Set<Node> children = new HashSet<Node>();

  /**
   * Returns a reference to a list of all children of the node (modifications
   * will be reflected in the original Node object).
   */
  public Set<? extends Node> getChildren() {
    return children;
  }

  /**
   * Creates a new Set containing references to all children of the Node, all of
   * their children and so on, searching recursively through the graph.
   */
  public Set<? extends Node> getAllChildren() {
    Set<Node> myChildrenList = new HashSet<Node>();
    getAllChildrenHelper(myChildrenList);
    return myChildrenList;
  }

  /** A helper function for getAllChildren() recursive search. */
  private final void getAllChildrenHelper(Set<Node> setParam) {
    for ( Node curNode : children ) {
      setParam.add( curNode );
      curNode.getAllChildrenHelper(setParam);
    }
  }

  /**
   * Creates a new Set containing references to all the children of any parent 
   * of the given Node. When looking for children, looks only one level deep. 
   * The returned collection does not include the given node.
   */
  public Set<? extends Node> getSiblings() {
    Set<Node> mySibs = new HashSet<Node>();
    for ( Node n : parents ) {
      mySibs.addAll( n.children );
    }
    mySibs.remove(this);
    return mySibs;
  }


  /**
   * Gives back the name of the Node, but with any characters illegal
   * in filenames, commas and spaces replaced with "_", and double quotes with
   * single quotes.
   */
  public String getSafeName() {
    String result = getName();
    if (result != null) {
      result = result.replace('\\', '_');
      result = result.replace('/', '_');
      result = result.replace('*', '_');
      result = result.replace('?', '_');
      result = result.replace('"', '\'');
      result = result.replace('<', '_');
      result = result.replace('>', '_');
      result = result.replace('|', '_');
      result = result.replace(':', '_');
      result = result.replace(' ', '_');
      result = result.replace(',', '_');
      result = result.trim();
    }
    return result;
  }


  /** List of alternate names for the node. */
  protected List<String> altNames = new ArrayList<String>();

  /** Returns a list of alternate names for the node. */
  public List<String> getAltNames() {
    return altNames;
  }


  /** This speeds up topmost parent lookups on very deep graphs. */
  protected Node cachedTopLevelNode = null;

  /**
   * Finds the topmost parent of the given node (i.e. the one that has no
   * parents.)
   */
  public Node getTopNode() {
    if (cachedTopLevelNode == null) {
      Node curNode = this;
      while (curNode.parents.size() > 0) {
        curNode = (Node) curNode.parents.toArray()[0];
      }
      cachedTopLevelNode = curNode;
    }
    return cachedTopLevelNode;
  }


  /** Returns true if the node is a top-level node (if it has no parents). */
  public boolean isTopLevel() {
    return this.parents.isEmpty();
  }
  
  /** Returns true if the node is a terminal node (if it has no children). */
  public boolean isTerminal() {
    return this.children.isEmpty();
  }
  

  /**
   * Used to store arbitrary properties for a node. Quite inefficient
   * because values are stored in HashMaps indexed by Strings. Use sparingly.
   */
  protected final Map<String, Double> properties = new LinkedHashMap<String, Double>();

  public Map<String, Double> getProperties() {
    return properties;
  }


  /**
   * Two nodes are equal if they have equal unique IDs. Note that lists of
   * alternate IDs (field: altIds) are not taken into account by this
   * equals() implementation.
   *
   * @param obj A node to compare.
   * @return true if two nodes have the same main unique ID.
   */
  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    final Node other = (Node) obj;
    if (this.id != other.id) {
      return false;
    }
    return true;
  }


  /** hashCode is based only on the unique ID of a node. */
  @Override
  public int hashCode() {
    int hash = 3;
    hash = 89 * hash + this.id;
    return hash;
  }


  @Override
  public String toString() {
    return this.id + ":" + this.name;
  }


  public boolean isTraversable() {
    return traversable;
  }


  public void setTraversable(boolean traversable) {
    this.traversable = traversable;
  }


}
