package speedemon;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.beauti.Beauti;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

@Description("Mixture of birth-death skyline model "
		   + "and spike distribution on internal node heights")
public class BirthDeathSkylineCollapseModel extends BirthDeathSkylineModel {
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

		TreeAboveThreshold treeAboveThreshold = new TreeAboveThreshold();
		treeAboveThreshold.assignFrom((Tree) tree);
		treeAboveThreshold.filterTree(tree, epsilon);
		logP = super.calculateTreeLogLikelihood(treeAboveThreshold);
		
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
    
	
	
	
	public double calculateTreeLogLikelihood(Tree tree, double epsilon) {

        logP = 0.;

        int nTips = tree.getLeafNodeCount();

        if (preCalculation(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // number of lineages at each time ti
        int[] n = new int[totalIntervals];

        int index = 0;
        if (times[index] < 0.)
            index = index(0.);

        double x0 = 0.;
        double temp = 0.;

        switch (conditionOn) {
            case NONE:
                temp = log_q(index, times[index], x0);
                break;
            case SURVIVAL:
                temp = p0(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                if (conditionOnRootInput.get()) {
                    temp = log_q(index, times[index], x0) - 2 * Math.log(1 - temp) - Math.log(birth[index]);
                } else {
                    temp = log_q(index, times[index], x0) - Math.log(1 - temp);
                }
                break;
            case RHO_SAMPLING:
                temp = p0hat(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                if (conditionOnRootInput.get()) {
                    temp = log_q(index, times[index], x0) - 2 * Math.log(1 - temp);
                } else {
                    temp = log_q(index, times[index], x0) - Math.log(1 - temp);
                }
                break;
            default:
                break;
        }

        logP = temp;
        if (Double.isInfinite(logP))
            return logP;

        if (taxonInput.get() != null) {
            if (taxonAge > origin.get().getValue()) {
                return Double.NEGATIVE_INFINITY;
            }
            double x = times[totalIntervals - 1] - taxonAge;
            index = index(x);
            if (SATaxonInput.get().getValue() == 0) {
                logP += Math.log(p0(index, times[index], x));
            } else {
                logP += Math.log(1-p0(index, times[index], x));
            }

            return logP;
        }

        if (printTempResults) System.out.println("first factor for origin = " + temp);

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[totalIntervals - 1] - tree.getNode(nTips + i).getHeight();
            index = index(x);
            if (!(tree.getNode(nTips + i)).isFake()) {
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; interval = " + i);
                if (Double.isInfinite(logP))
                    return logP;
            }
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            if (!isRhoTip[i] || m_rho.get() == null) {
                double y = times[totalIntervals - 1] - tree.getNode(i).getHeight();
                index = index(y);

                if (!(tree.getNode(i)).isDirectAncestor()) {
                    if (!SAModel) {
                        temp = Math.log(psi[index]) - log_q(index, times[index], y);
                    } else {
                        temp = Math.log(psi[index] * (r[index] + (1 - r[index]) * p0(index, times[index], y))) - log_q(index, times[index], y);
                    }
                    logP += temp;
                    if (printTempResults) System.out.println("2nd PI = " + temp);
                    if (psi[index] == 0 || Double.isInfinite(logP))
                        return logP;
                } else {
                    if (r[index] != 1) {
                        logP += Math.log((1 - r[index])*psi[index]);
                        if (Double.isInfinite(logP)) {
                            return logP;
                        }
                    } else {
                        //throw new Exception("There is a sampled ancestor in the tree while r parameter is 1");
                        System.out.println("There is a sampled ancestor in the tree while r parameter is 1");
                        System.exit(0);
                    }
                }
            }
        }

        // last product term in f[T], factorizing from 1 to m //
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
            int[] k = {0};
            if (!SAModel) {
                n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree));
            } else {
                n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree, k));
            }
            if (n[j] > 0) {
                temp = n[j] * (log_q(j, times[j], time) + Math.log(1 - rho[j-1]));
                logP += temp;
                if (printTempResults)
                    System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }

            if (SAModel && j>0 && N != null) { // term for sampled leaves and two-degree nodes at time t_i
                logP += k[0] * (log_q(j, times[j], time) + Math.log(1-r[j])) + //here g(j,..) corresponds to q_{i+1}, r[j] to r_{i+1},
                        (N[j-1]-k[0])*(Math.log(r[j]+ (1-r[j])*p0(j, times[j], time))); //N[j-1] to N_i, k[0] to K_i,and thus N[j-1]-k[0] to M_i
                if (Double.isInfinite(logP)) {
                    return logP;
                }
            }

            if (rho[j] > 0 && N[j] > 0) {
                temp = N[j] * Math.log(rho[j]);    // term for contemporaneous sampling
                logP += temp;
                if (printTempResults)
                    System.out.println("3rd factor (Nj loop) = " + temp + "; interval = " + j + "; N[j] = " + N[j]);
                if (Double.isInfinite(logP))
                    return logP;

            }
        }

        if (SAModel) {
            int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
            logP +=  Math.log(2)*internalNodeCount;
        }

        return logP;
    }
    
}
