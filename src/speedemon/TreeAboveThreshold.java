package speedemon;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

public class TreeAboveThreshold extends Tree {

	
	/**
	 * Filter out some of the taxa so that all internal nodes
	 * are above epsilon.
	 * This assumes the tree does not have dated tips
	 * @param tree: tree to be filtered, it is assumed all taxa
	 * are in the taxon set of this tree.
	 * @param epsilon: threshold level below which nodes are removed
	 */
	public void filterTree(TreeInterface tree, double epsilon) {
		for (Node node : m_nodes) {
			node.removeAllChildren(false);
			node.setParent(null);
		}
		
		leafNodeCount = 0;
		copyFrom(tree.getRoot(), epsilon);
		
		// reorder m_nodes array
		nodeCount = leafNodeCount;
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			if (m_nodes[i].getChildCount() != 0) {
				// swap m_nodes[i] with latest unused node in m_nodes
				Node tmp = m_nodes[nodeCount];
				m_nodes[nodeCount] = m_nodes[i];
				m_nodes[i] = tmp;
				m_nodes[nodeCount].setNr(nodeCount);
				nodeCount++;
			}
		}
	}

	/**
	 * Copy part of the tree that is above epsilon
	 */
	private Node copyFrom(Node node, double epsilon) {
		if (node.isLeaf()) {
			leafNodeCount++;
			m_nodes[leafNodeCount-1].setHeight(node.getHeight());			
			return m_nodes[leafNodeCount-1];
		} else {
			Node left = copyFrom(node.getLeft(), epsilon);
			if (node.getHeight() <= epsilon) {
				return left;
			}
			Node right = copyFrom(node.getLeft(), epsilon);
			
			Node thisNode = m_nodes[node.getNr()];
			thisNode.addChild(left);
			thisNode.addChild(right);
			return thisNode;
		}
	}
	
	
	
}
