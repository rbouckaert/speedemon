package speedemon.snapper;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;
import speedemon.ClusterCounter;
import speedemon.ClusterTreeSetAnalyser;

@Description("Reversible jump a tree to collapse/expand its leaf branches")
public class ClusterOperator extends Operator {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	private TreeInterface tree;
	private int [] map;
	private boolean [] clustered;
	
	
	private double [] logHR;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		map = new int[tree.getLeafNodeCount()];
		clustered = new boolean[tree.getNodeCount()];

		
		int n = tree.getLeafNodeCount();
		
		// treecount[k] = number of trees with k taxa 
		BigInteger [] treecount = new BigInteger[n+1];
		treecount[0] = BigInteger.ZERO;
		treecount[1] = BigInteger.ONE;
		long t = 1l;
		for (int k = 0; k < n-1; k++) {
			if (t < 0) {
				throw new IllegalArgumentException("Too many taxa leading to overflow");
			}
			// # (unranked) tree topologies
			// t = t*(2*(k+2)-3);
			
			// # (ranked) tree topologies
			t = ((k+1) * (k+2)/2) * t;
			BigInteger s2 = sterling2(n, k+2);
			treecount[k+2] = s2.multiply(new BigInteger(t+""));
			// treecount[k+2] = new BigInteger(t+"");
		}
		logHR = new double[n+1];
		for (int k = 1; k < n; k++) {
			BigDecimal d = new BigDecimal(treecount[k+1]).divide(new BigDecimal(treecount[k]), MathContext.DECIMAL32);
			logHR[k] = Math.log(d.doubleValue());
		}
		
		System.out.println(Arrays.toString(treecount));
		System.out.println(Arrays.toString(logHR));
	}

	
	private static Map<String,BigInteger> COMPUTED = new HashMap<>();
	// from https://rosettacode.org/wiki/Stirling_numbers_of_the_second_kind#Java
    private static final BigInteger sterling2(int n, int k) {
        String key = n + "," + k;
        if ( COMPUTED.containsKey(key) ) {
            return COMPUTED.get(key);
        }
        if ( n == 0 && k == 0 ) {
            return BigInteger.valueOf(1);
        }
        if ( (n > 0 && k == 0) || (n == 0 && k > 0) ) {
            return BigInteger.ZERO; 
        }
        if ( n == k ) {
            return BigInteger.valueOf(1);
        }
        if ( k > n ) {
            return BigInteger.ZERO;
        }
        BigInteger result = BigInteger.valueOf(k).multiply(sterling2(n-1, k)).add(sterling2(n-1, k-1));
        COMPUTED.put(key, result);
        return result;
    }
    
    
	@Override
	public double proposal() {
        final TreeInterface tree = treeInput.get(this);

        List<Node> mergeCandidates = new ArrayList<>();
        List<Node> splitCandidates = new ArrayList<>();
		boolean doSplit = Randomizer.nextBoolean();
		Node [] nodes = tree.getNodesAsArray();


		
		Arrays.fill(clustered, false);
		int  k = ClusterCounter.countClusters(tree, map, clustered);
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			if (nodes[i].getLength() <= ClusterTreeSetAnalyser.EPSILON) {
				clustered[nodes[i].getParent().getNr()] = true;
			}
		}
		
		if (k == tree.getLeafNodeCount() && doSplit) {
			return Double.NEGATIVE_INFINITY;
		} else if (k == 2 && !doSplit) {
			return Double.NEGATIVE_INFINITY;
		}
							
		if (!doSplit) {
			// merge
			for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount() - 1; i++) {
				Node left = nodes[i].getLeft();
				Node right = nodes[i].getRight();
				if ((left.isLeaf() || clustered[left.getNr()]) && 
					(right.isLeaf() || clustered[right.getNr()]) && 
					!clustered[i]) {
					mergeCandidates.add(nodes[i]);
				}
			}

			if (mergeCandidates.size() <= 0) {
				throw new RuntimeException("Programmer error -- should not get here");
				// return Double.NEGATIVE_INFINITY;
			}
			// Node node = mergeCandidates.get(Randomizer.nextInt(mergeCandidates.size()));
			
			// merge lowest ranked candidate
			Node node = mergeCandidates.get(0);
			for (Node c : mergeCandidates) {
				if (c.getHeight() < node.getHeight()) {
					node = c;
				}
			}
			node.setHeight(node.getLeft().getHeight());

			
			clustered[node.getNr()] = true;
			for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount() - 1; i++) {
				node = nodes[i];
				Node left = node.getLeft();
				Node right = node.getRight();
				if ((left.isLeaf() || clustered[left.getNr()]) && 
					(right.isLeaf() || clustered[right.getNr()]) &&
					clustered[node.getNr()] && 
					node.getLength() > ClusterTreeSetAnalyser.EPSILON) {
					splitCandidates.add(nodes[i]);
				}
			}		
			
			// j = number of ways to split
			int j = 0;
			for (Node n : splitCandidates) {
				int m = n.getLeafNodeCount();
				j += Math.pow(2,m-1) - 1;
			}
			j = 1;
			// k * (k-1)/2) = number of ways to merge from k clusters = k choose 2
			return -logHR[k-1];
			// return -logHR[k-1] - (Math.log(j) - Math.log(k * (k-1)/2));
			// return - Math.log(j);
		} else {
			// split
			for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount() - 1; i++) {
				Node node = nodes[i];
				Node left = node.getLeft();
				Node right = node.getRight();
				if ((left.isLeaf() || clustered[left.getNr()]) && 
					(right.isLeaf() || clustered[right.getNr()]) &&
					clustered[node.getNr()] && 
					node.getLength() > ClusterTreeSetAnalyser.EPSILON) {
					splitCandidates.add(nodes[i]);
				}
			}		

			if (splitCandidates.size() == 0) {
				throw new RuntimeException("Programmer error -- should not get here");
				// return Double.NEGATIVE_INFINITY;
			}
			Node node = splitCandidates.get(Randomizer.nextInt(splitCandidates.size()));

			double hMax = tree.getRoot().getHeight();
			for (Node n : tree.getInternalNodes()) {
				if (n.getHeight() > ClusterTreeSetAnalyser.EPSILON && n.getHeight() < hMax) {
					hMax = n.getHeight();
				}
			}			
			//double h = node.getLeft().getHeight() + Randomizer.nextDouble() * node.getLength();
			double h = node.getLeft().getHeight() + Randomizer.nextDouble() * hMax;
			node.setHeight(h);
			
			// j = number of ways to split
			int j = 0;
			for (Node n : splitCandidates) {
				int m = n.getLeafNodeCount();
				j += Math.pow(2,m-1) - 1;
			}
			j = 1;
			// (k+1) * k/2) = number of ways to merge from k+1 clusters = k+1 choose 2
			return logHR[k];
			// return logHR[k] + (Math.log(j) - Math.log((k+1) * k/2));
			// return + Math.log(j);
		}
		

	}

}
