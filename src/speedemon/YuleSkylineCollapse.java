package speedemon;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.beauti.Beauti;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import biceps.YuleSkyline;

@Description("Mixture of skyline version of Yule tree prior that integrates out birth rate parameters"
		+ " under a gamma prior"
		+ " and spike distribution on internal node heights")
@Citation(value="Jordan Douglas and Remco Bouckaert. Quantitatively defining species boundaries with more efficiency and more biological realism. Communications Biology 5, 755 (2022)", DOI="110.1038/s42003-022-03723-z")
public class YuleSkylineCollapse extends YuleSkyline {

    final public Input<Function> collapseHeightInput = new Input<>("epsilon", "collapse height value below wich taxa are considered to be the same species.", Validate.REQUIRED);
    final public Input<RealParameter> collapseWeightInput =  new Input<>("weight", "mixture weight between Yule and spike density.", Validate.REQUIRED);

    private RealParameter weight;
    private Function epsilon;
    private TreeInterface tree;

    @Override
    public void initAndValidate() {
    	if (ProgramStatus.name.equals("BEAUti")) {
    		return;
    	}
    	super.initAndValidate();
    	epsilon = collapseHeightInput.get();
    	ClusterTreeSetAnalyser.EPSILON = epsilon.getArrayValue();
    	weight = collapseWeightInput.get();
		tree = treeInput.get() == null ?
				treeIntervalsInput.get().treeInput.get():
					treeInput.get();
    }


    
	@Override
	public double calculateLogP() {
		double epsilon = this.epsilon.getArrayValue();

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
					while (!node.isRoot() && node.getParent().getHeight() <= epsilon) {
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
		int clusterCount = countClusters(tree, map, done, epsilon.getArrayValue());
		out.append(clusterCount + "\t");
	}




	
    @Override
    public List<String> getConditions() {
    	List<String> conditions = new ArrayList<>();
    	conditions.add(treeInput.get().getID());
    	conditions.add(birthRateShapeInput.get().getID());
    	conditions.add(birthRateRateInput.get().getID());
    	conditions.add(collapseWeightInput.get().getID());
    	// conditions.add(collapseHeightInput.get().getID());
    	return conditions;
    }
	
	@Override
	public void sample(State state, Random random) {
		


        if (sampledFlag) return;
        sampledFlag = true;
        
    	
		if (useEqualEpochs) {
			//throw new IllegalArgumentException(YuleSkylineCollapse.class.getCanonicalName() + ": please ensure that non equal epochs are used (set equalEpochs to false)");
		}
		
		if (linkedMeanInput.get()) {
			throw new IllegalArgumentException(YuleSkylineCollapse.class.getCanonicalName() + ": please ensure that linkedMean is false");
		}
		
		if (groupCount != 1) {
			throw new IllegalArgumentException(YuleSkylineCollapse.class.getCanonicalName() + ": please ensure that groupCount is 1");
		}
		
		
		
        
        if (!isPrepared) {
            prepare();
        }

        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        Tree tree = (Tree) treeInput.get();
        double w = collapseWeightInput.get().getValue();
        double epsilon = collapseHeightInput.get().getArrayValue();

        
        // Simulate tree conditional on new parameters
        List<Node> activeLineages = new ArrayList<>();
        for (Node oldLeaf : tree.getExternalNodes()) {
            Node newLeaf = new Node(oldLeaf.getID());
            newLeaf.setNr(oldLeaf.getNr());
            newLeaf.setHeight(0.0);
            activeLineages.add(newLeaf);

        }
        
        
        
        // How many of the internal nodes will be below epsilon? Binomial(n-1, w) distribution. Assuming binary tree
        List<Double> collapseHeights = new ArrayList<>();
        int n = activeLineages.size();
        for (int i = 0; i < n-1; i ++) {
        	if (random.nextDouble() < w) {
        	
        		// Sample a collapse height
        		double h = random.nextDouble() * epsilon;
        		collapseHeights.add(h);
        		
        	}
        	
        }
	
	
        int nextNr = activeLineages.size();
	
		// Create collapse epoch
		Collections.sort(collapseHeights);
		for (int i = 0; i < collapseHeights.size(); i ++) {
			
			
			double t = collapseHeights.get(i);
			int k = activeLineages.size();
			
			// Sample 2 nodes
			Node node1 = activeLineages.get(random.nextInt(k));
            Node node2;
            do {
                node2 = activeLineages.get(random.nextInt(k));
            } while (node2.equals(node1));
            
            
            
            // Join them
            Node newParent = new Node();
            newParent.setNr(nextNr++);
            newParent.setHeight(t);
            newParent.addChild(node1);
            newParent.addChild(node2);

            activeLineages.remove(node1);
            activeLineages.remove(node2);
            activeLineages.add(newParent);
  
			
		}

        

		// Sample Yule birth rates
		double[][] rates = this.sampleBirthRates();
	    double[] birthRates = rates[0];
	    
	    // Epoch 0
	    int epochNumber = 0;
	    double birthRate = birthRates[epochNumber];
        int groupSize = (int)groupSizes.getArrayValue(epochNumber);
	    
	    
        
		// Sample from Yule beginning at time epsilon
        double t = epsilon;
        while (activeLineages.size() > 1) {
            int k = activeLineages.size();
            
            
            // Proceed to the next epoch
            if (k <= groupSize) {
            	epochNumber++;
            	birthRate = birthRates[epochNumber];
                groupSize = (int)groupSizes.getArrayValue(epochNumber);
            }
    
    		
  
    		// Sample from Yule (exponential distribution)
			double a = birthRate * k;
			double t_ = -Math.log(random.nextDouble())/a;
            t += t_;

            Node node1 = activeLineages.get(random.nextInt(k));
            Node node2;
            do {
                node2 = activeLineages.get(random.nextInt(k));
            } while (node2.equals(node1));

            Node newParent = new Node();
            newParent.setNr(nextNr++);
            newParent.setHeight(t);
            newParent.addChild(node1);
            newParent.addChild(node2);

            activeLineages.remove(node1);
            activeLineages.remove(node2);
            activeLineages.add(newParent);
        }

        tree.assignFromWithoutID(new Tree(activeLineages.get(0)));
		

		
		
	}

	
	
}
