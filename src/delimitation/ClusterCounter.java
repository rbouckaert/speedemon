package delimitation;

import java.io.PrintStream;
import java.util.Arrays;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Allows logging of the number of clusters")
public class ClusterCounter extends BEASTObject implements Loggable {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    
	private TreeInterface tree;
	private int [] map;
	private boolean [] done;
	
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		map = new int[tree.getLeafNodeCount()];
		done = new boolean[tree.getLeafNodeCount()];
	}

	@Override
	public void init(PrintStream out) {
		out.append("ClusterCount\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		int k = countClusters(tree, map, done);
		out.append(k + "\t");
	}
	
	@Override
	public void close(PrintStream out) {
	}

	
	public static int countClusters(TreeInterface tree, int [] map, boolean [] done) {
		Node [] nodes = tree.getNodesAsArray();
		
		Arrays.fill(map, -1);
		Arrays.fill(done, false);

		int k = 0;
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			if (!done[i]) {
				if (nodes[i].getHeight() > ClusterTreeSetAnalyser.EPSILON) {
					map[i] = k;
					done[i] = true;
				} else {
					// nodes[i] is part of a cluster
					Node node = nodes[i];
					while (!node.getParent().isRoot() && node.getHeight() <= ClusterTreeSetAnalyser.EPSILON) {
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
}
