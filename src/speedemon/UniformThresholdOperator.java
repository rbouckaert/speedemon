package speedemon;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.Uniform;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;


/*
 * Selects a node whose parent or child is on other side of threshold (at height epsilon) and moves the node across the threshold to a position uniformly at random
 * This will add/remove a cluster
 */
public class UniformThresholdOperator extends Uniform {
	
	

	final public Input<Function> epsilonInput = new Input<>("epsilon", "the threshold parameter.", Validate.REQUIRED);
	
	
	
	  /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        Tree tree = treeInput.get();
        double epsilon = epsilonInput.get().getArrayValue();

        
        // Select a node above or below threshold?
        boolean startAboveThreshold = Randomizer.nextBoolean();
       
        
        // Get candidate nodes
        List<Node> eligibleNodesBefore = getCandidates(tree, epsilon, startAboveThreshold);
        int nOptionsBefore = eligibleNodesBefore.size();
    	
        
       // System.out.println("nOptionsBefore " + nOptionsBefore + " " + startAboveThreshold);
        
    	
    	// If not possible, then reject proposal
    	if (nOptionsBefore == 0) return Double.NEGATIVE_INFINITY;
    	
    	
    	
    	
    	
    	int num = Randomizer.nextInt(eligibleNodesBefore.size());
    	Node node = eligibleNodesBefore.get(num);
    	
    	
    	// Moving up or down?
        double upperBefore, lowerBefore;
        if (startAboveThreshold) {
        	upperBefore = epsilon;
        	lowerBefore = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        }else {
        	upperBefore = node.getParent().getHeight();
        	lowerBefore = epsilon;
        }
    	

        
       

        // Move value into range uniformly at random
        final double newValue = (Randomizer.nextDouble() * (upperBefore - lowerBefore)) + lowerBefore;
        node.setHeight(newValue);
        
        
        
        // Hastings ratio contribution from uniform height sampling
        double upperAfter, lowerAfter;
        if (startAboveThreshold) {
        	upperAfter = node.getParent().getHeight();
        	lowerAfter = epsilon;
        }else {
        	upperAfter = epsilon;
        	lowerAfter = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        }
        double logHR = Math.log(upperBefore-lowerBefore) - Math.log(upperAfter-lowerAfter);
        
        
        // Hastings ratio contribution from node sampling
        List<Node> eligibleNodesAfter = getCandidates(tree, epsilon, !startAboveThreshold);
        int nOptionsAfter = eligibleNodesAfter.size();
        logHR += Math.log(nOptionsBefore) - Math.log(nOptionsAfter);
        
        
       // System.out.println("nOptionsAfter " + nOptionsAfter);
        
        
        
        return logHR;
    }

    
    
    /**
     * Get a list of candidate nodes above/below boundary
     * @param movingDown
     * @return
     */
    public List<Node> getCandidates(Tree tree, double epsilon, boolean startAboveThreshold){
    	
    	
    	  
        // Get list of internal nodes above threshold, whose children are all below threshold
    	List<Node> eligibleNodes = new ArrayList<>();
    	for (Node parent : tree.getInternalNodes()) {
    		
    		if (parent.getHeight() < epsilon) continue;
    		
    		
    		boolean allChildrenBelowEpsilon = true;
    		for (Node c : parent.getChildren()) {
    			if (c.getHeight() > epsilon) {
    				allChildrenBelowEpsilon = false;
    			}else if (!startAboveThreshold && !c.isLeaf()) {
    				eligibleNodes.add(c);
    			}
    		}
    		
    		
    		if (startAboveThreshold && allChildrenBelowEpsilon && !parent.isRoot()) eligibleNodes.add(parent);
    		
    		
    	}
    	
    	
    	return eligibleNodes;
    	
    }
    
    
    public List<StateNode> listStateNodes() {
        // pick up all inputs that are stateNodes that are estimated
        final List<StateNode> list = new ArrayList<>();
        list.add(treeInput.get());
        return list;
    }
    
	
	

}
