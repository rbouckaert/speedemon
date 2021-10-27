package delimitation;


import java.io.PrintStream;
import java.util.Arrays;

import beast.app.beauti.Beauti;
import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import biceps.YuleSkyline;

@Description("Mixture of skyline version of Yule tree prior that integrates out birth rate parameters"
		+ " under a gamma prior"
		+ " and spike distribution on internal node heights")
public class YuleSkylineCollapse extends YuleSkyline {

    final public Input<RealParameter> collapseHeightInput = new Input<>("epsilon", "collapse height value below wich taxa are considered to be the same species.", Validate.REQUIRED);
    final public Input<RealParameter> collapseWeightInput =  new Input<>("weight", "mixture weight between Yule and spike density.", Validate.REQUIRED);

    private RealParameter weight;
    private RealParameter epsilon;
    private TreeInterface tree;

    @Override
    public void initAndValidate() {
    	if (Beauti.isInBeauti()) {
    		return;
    	}
    	super.initAndValidate();
    	epsilon = collapseHeightInput.get();
    	ClusterTreeSetAnalyser.EPSILON = epsilon.getValue();
    	weight = collapseWeightInput.get();
		tree = treeInput.get() == null ?
				treeIntervalsInput.get().treeInput.get():
					treeInput.get();
    }


    
	@Override
	public double calculateLogP() {
		double epsilon = this.epsilon.getValue();

		logP = super.calculateLogPbyEqualEpochs(epsilon);
		
		double w = this.weight.getValue();
		
		int k = 0; // number of node heights >= epsilon
		int n = tree.getInternalNodeCount(); // number of internal nodes
		for (int i = 0; i < n; i++) {
			if (tree.getNode(n + 1 + i).getHeight() >= epsilon) {
				k++;
			}
		}

		logP += k * Math.log(1-w) + (n-k) * Math.log(w/epsilon);
	
		return logP;
    }

	

	public static int countClusters(TreeInterface tree, int [] map, boolean [] done, double epsilon) {
		Node [] nodes = tree.getNodesAsArray();
		
		Arrays.fill(map, -1);
		Arrays.fill(done, false);

		int k = 0;
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			if (!done[i]) {
				if (nodes[i].getLength() >= epsilon) {
					map[i] = k;
					done[i] = true;
				} else {
					// nodes[i] is part of a cluster
					Node node = nodes[i];
					while (!node.getParent().isRoot() && node.getLength() <= epsilon) {
						node = node.getParent();
					}
					visit(node,k,map,done);
				}
				k++;
			}
		}
		return k;
	}
	
    private static void visit(Node node, int k, int[] map, boolean[] done) {
		if (node.isLeaf()) {
			int i = node.getNr();
			map[i] = k;
			done[i] = true;
		} else {
			for (Node child : node.getChildren()) {
				visit(child, k, map, done);
			}
		}
	}

    
    
	@Override
	public void init(PrintStream out) {
		super.init(out);
		
		out.append("clusterCount\t");
	}
	
	@Override
	public void log(long sampleNr, PrintStream out) {
		super.log(sampleNr, out);
		
		int [] map = new int[tree.getLeafNodeCount()];
		boolean [] done = new boolean[tree.getLeafNodeCount()];
		int clusterCount = countClusters(tree, map, done, epsilon.getValue());
		out.append(clusterCount + "\t");
	}



	
	
	
}
