package delimitation.snapper;

import java.util.Arrays;

import beast.core.Description;
import beast.evolution.tree.TreeInterface;
import delimitation.ClusterCounter;
import snapper.SnapperTreeLikelihood;


@Description("As SnapperTreeLikelihood but the likelihood correction is updated when "
		+ "branches to leafs collapse to zero, making them form a cluster")
// TODO: only recalculate when `map` changes
// TODO: make tree traversal more efficient
public class ClusterableSnapperTreeLikelihood extends SnapperTreeLikelihood {

	public ClusterableSnapperTreeLikelihood() throws Exception {
		super();
	}

	private TreeInterface tree;
	private int [] map, sites, counts;
	private boolean [] done;
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		tree = treeInput.get();
		map = new int[tree.getLeafNodeCount()];
		done = new boolean[tree.getLeafNodeCount()];
		
		sites = new int[tree.getLeafNodeCount()];
		counts = new int[tree.getLeafNodeCount()];
	}
	
		
	@Override
	public double calculateLogP() {
		m_fLogLikelihoodCorrection = calcCorrection();
		logP = super.calculateLogP();
		return logP;
	}


	private double calcCorrection() {
		int k = ClusterCounter.countClusters(tree, map, done);
		
		
		// calculate Likelihood Correction. 
		// When the assignment of individuals to populations/species is fixed, the allele counts in each population are sufficient 
		// statistics for the species tree parameters. However when testing species assignments this is no longer the case.
		// To address this we multiply the likelihood computed from allele counts by the probability of observing
		// the given sequences given those allele counts (and the species assignments).
		m_fLogLikelihoodCorrection = 0;
			// RRB: note that increasing the number of constant sites
			// does not change the m_fLogLikelihoodCorrection since the
			// contribution of constant sites is zero. This means,
			// m_fLogLikelihoodCorrection does not need to be recalculated
			// when ascSiteCount changes.
			// DJB: This is true, but only until we start looking at non-constant sites being ascertained.

			
			int numPatterns = m_data2.getPatternCount();
	    	for (int i = 0; i < numPatterns; i++) {
	            int [] thisSite = m_data2.getPattern(i);  //count of red alleles for this site
	            int [] lineageCounts = m_data2.getPatternLineagCounts(i); //count of total lineages for this site
	            
	            
				Arrays.fill(sites, 0);
				Arrays.fill(counts, 0);
	            for (int j = 0; j < thisSite.length; j++) {
	            	sites[map[j]] += thisSite[j];
	            	counts[map[j]] += lineageCounts[j];
	            }
	            for (int j = 0; j < k; j++) {
	            	m_fLogLikelihoodCorrection -= logBinom(sites[j], counts[j]) * m_data2.getPatternWeight(i);
	            }
	    	}
		return 0;
	}
	



	

	private double logBinom(int k, int n) {
    	double f = 0;
    	for (int i = k + 1; i <= n; i++) {
    		f += Math.log(i) - Math.log(n - i + 1);
    	}
		return f;
	}

}
