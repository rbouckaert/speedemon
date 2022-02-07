package speedemon;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.BactrianNodeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


/*
 * Selects a node whose parent or child is on other side of threshold (at height epsilon) and moves the node across the threshold
 * This will add/remove a cluster
 */
public class BactrianThresholdOperator extends BactrianNodeOperator {
	
	

	final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "the threshold parameter.", Validate.REQUIRED);
	
	
	
	  /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        Tree tree = treeInput.get(this);
        double epsilon = epsilonInput.get(this).getValue();

        
        // Select a node above or below threshold?
        boolean startAboveThreshold = Randomizer.nextBoolean();
       
        
        // Get candidate nodes
        List<Node> eligibleNodesBefore = getCandidates(tree, epsilon, startAboveThreshold);
        int nOptionsBefore = eligibleNodesBefore.size();
    	
    	
    	// If not possible, then reject proposal
    	if (nOptionsBefore == 0) return Double.NEGATIVE_INFINITY;
    	
    	
    	
    	
    	
    	int num = Randomizer.nextInt(eligibleNodesBefore.size());
    	Node node = eligibleNodesBefore.get(num);
    	
    	
    	// Moving up or down?
        double upper, lower;
        if (startAboveThreshold) {
        	upper = epsilon;
        	lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        }else {
        	upper = node.getParent().getHeight();
        	lower = epsilon;
        }
    	

        
        // Scale into range, using Bactrian distribution and a logit transformation
        double scale = kernelDistribution.getScaler(0, Double.NaN, scaleFactor);

        // transform value
        double value = node.getHeight();
        double y = (upper - value) / (value - lower);
        y *= scale;
        double newValue = (upper + lower * y) / (y + 1.0);
        
        if (newValue < lower || newValue > upper) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        node.setHeight(newValue);
        
        
        
        // Hastings ratio contribution from random walk
        double logHR = Math.log(scale) + 2.0 * Math.log((newValue - lower)/(value - lower));
        
        
        // Hastings ratio contribution from node sampling
        List<Node> eligibleNodesAfter = getCandidates(tree, epsilon, !startAboveThreshold);
        int nOptionsAfter = eligibleNodesAfter.size();
        logHR += Math.log(nOptionsBefore) - Math.log(nOptionsAfter);
        
        
        
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
    	for (Node node : tree.getInternalNodes()) {
    		
    		if (node.getHeight() < epsilon) continue;
    		
    		
    		boolean allChildrenBelowEpsilon = true;
    		for (Node c : node.getChildren()) {
    			if (c.getHeight() > epsilon) {
    				allChildrenBelowEpsilon = false;
    			}else if (!startAboveThreshold && !c.isLeaf()) {
    				eligibleNodes.add(c);
    			}
    		}
    		
    		
    		if (startAboveThreshold && allChildrenBelowEpsilon) eligibleNodes.add(node);
    		
    		
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
