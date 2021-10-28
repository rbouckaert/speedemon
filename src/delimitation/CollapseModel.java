package delimitation;

import java.io.PrintStream;

import beast.app.beauti.Beauti;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

@Description("Mixture of a custom tree prior to be used above threshold "
		   + "and spike distribution on internal node heights")
public class CollapseModel extends SpeciesTreeDistribution {
	final public Input<RealParameter> collapseHeightInput = new Input<>("epsilon", "collapse height value below wich taxa are considered to be the same species.", Validate.REQUIRED);
    final public Input<RealParameter> collapseWeightInput =  new Input<>("weight", "mixture weight between Yule and spike density.", Validate.REQUIRED);
	final public Input<SpeciesTreeDistribution> treepriorInput = new Input<>("treePrior", "Tree prior to be used on tree filtered to be above threshold", Validate.REQUIRED);

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

		TreeAboveThreshold treeAboveThreshold = new TreeAboveThreshold();
		treeAboveThreshold.assignFrom((Tree) tree);
		treeAboveThreshold.filterTree(tree, epsilon);
		logP = treepriorInput.get().calculateTreeLogLikelihood(treeAboveThreshold);
		
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

	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {
		return this.calculateLogP();
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
		int clusterCount = YuleSkylineCollapse.countClusters(tree, map, done, epsilon.getValue());
		out.append(clusterCount + "\t");
	}
    


}
