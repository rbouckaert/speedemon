package delimitation.snapper;


import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import delimitation.ClusterTreeSetAnalyser;

/*
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */

@Description("Implements branch exchange operations. There is a NARROW and WIDE variety. " +
        "The narrow exchange is very similar to a rooted-beast.tree nearest-neighbour " +
        "interchange but with the restriction that node height must remain consistent.")
public class Exchange extends TreeOperator {
    final public Input<Boolean> isNarrowInput = new Input<>("isNarrow", "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

    @Override
    public void initAndValidate() {
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);

        double logHastingsRatio = 0;

        if (isNarrowInput.get()) {
            logHastingsRatio = narrow(tree);
        } else {
            logHastingsRatio = wide(tree);
        }

        return logHastingsRatio;
    }

    private int isg(final Node n) {
      return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     */
    public double narrow(final Tree tree) {

        final int internalNodes = tree.getInternalNodeCount();
        
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        int attempt = 0;
        do {
	        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf() ||
	        		grandParent.getHeight() <= ClusterTreeSetAnalyser.EPSILON) {
	            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
	        	attempt++;
	        	if (attempt > 20) {
	        		return Double.NEGATIVE_INFINITY;
	        	}
	            
	        }
	        
	        if (grandParent.getLeft().getHeight() <= ClusterTreeSetAnalyser.EPSILON) {
	        	// swap left & right
	        	Node left = grandParent.getLeft();
	        	grandParent.removeChild(left);
	        	grandParent.addChild(left);
	        	attempt++;
	        	if (attempt > 20) {
	        		return Double.NEGATIVE_INFINITY;
	        	}
	        }
        } while (grandParent.getLeft().getHeight() <= ClusterTreeSetAnalyser.EPSILON);
        

        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        if( parentIndex.isLeaf() ) {
            // tree with dated tips
            return Double.NEGATIVE_INFINITY;
        }

        int validGP = 0;
        {
            for(int i = internalNodes + 1; i < 1 + 2*internalNodes; ++i) {
                validGP += isg(tree.getNode(i));
            }
        }

        final int c2 = sisg(parentIndex) + sisg(uncle);

        final Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
        exchangeNodes(i, uncle, parentIndex, grandParent);

        final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

        return Math.log((float)validGP/validGPafter);
    }

    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     * @param tree
     */
    public double wide(final Tree tree) {

        final int nodeCount = tree.getNodeCount();

        Node i = tree.getRoot();

        while (i.isRoot()) {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        }

        Node j = i;
        while (j.getNr() == i.getNr() || j.isRoot()) {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
        }

        final Node p = i.getParent();
        final Node jP = j.getParent();

        if ((p != jP) && (i != jP) && (j != p)
                && (j.getHeight() < p.getHeight())
                && (i.getHeight() < jP.getHeight())) {
            exchangeNodes(i, j, p, jP);

            // All the nodes on the path from i/j to the common ancestor of i/j parents had a topology change,
            // so they need to be marked FILTHY.
            if( markCladesInput.get() ) {
                Node iup = p;
                Node jup = jP;
                while (iup != jup) {
                    if( iup.getHeight() < jup.getHeight() ) {
                        assert !iup.isRoot();
                        iup = iup.getParent();
                        iup.makeDirty(Tree.IS_FILTHY);
                    } else {
                        assert !jup.isRoot();
                        jup = jup.getParent();
                        jup.makeDirty(Tree.IS_FILTHY);
                    }
                }
            }
            return 0;
        }

        // Randomly selected nodes i and j are not valid candidates for a wide exchange.
        // reject instead of counting (like we do for narrow).
        return Double.NEGATIVE_INFINITY;
    }


    /* exchange sub-trees whose root are i and j */

    protected void exchangeNodes(Node i, Node j,
                                 Node p, Node jP) {
        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i
    }
}
