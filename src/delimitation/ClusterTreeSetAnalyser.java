package delimitation;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.*;

import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.FastTreeSet;
import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.app.util.TreeFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import delimitation.snapper.ClusterOperator;

@Description("Analyses the tree set, and compares with a tree if provided (e.g. the original used to simulate data from)\n" +
		   "-tree <newick tree>: tree in newick format on command line\n" +
		   "-b <burnin percentage>: percentage (so a number between 0 and 100) of tree at the start that are discarded, default 10%\n" +
		   "<tree set file>: name of the file containing tree set in NEXUS format\n" +
		   "Outputs:\n"+
		   "(1) Size of the 95% credible set.\n"+
		   "(2) Probabilities of the top 20 (?) trees\n"+
		   "(3) Whether or not the true tree is in the credible set."
		   )
public class ClusterTreeSetAnalyser extends beast.core.Runnable {
	public Input<TreeFile> treesInput = new Input<>("trees", "NEXUS file containing a tree set", new TreeFile("[[none]]"));
	public Input<TreeFile> originalInput = new Input<>("original", "tree to check if it is in the 95% credible set");
	public Input<OutFile> outputInput = new Input<>("out", "output file. Print to stdout if not specified", new OutFile("[[none]]"));
	public Input<Boolean> verboseInput = new Input<>("verbose", "print out extra information while processing", true);
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<Double> epsilonInput = new Input<>("epsilon", "threshold below which taxa are deemed to be clustered", EPSILON);
	public Input<Boolean> printTailInput = new Input<>("printTail", "print out clades and clusters outside 95% credible set", true);
	public Input<Integer> maxInput = new Input<>("max", "maximum number of trees/clusters to print out (if non-positive, "
			+ "all available data will be printed). Ignored if printTail=true", -1);

	final static private String SPECIES_SEPARATOR = "___";
	static public double EPSILON = 0;
	
	@Override
	public void initAndValidate() {
	}

	private double epsilon = EPSILON;
	private int max;
	
	@Override
	public void run() throws Exception {
		if (!verboseInput.get()) {
			Log.setLevel(Log.Level.error);
		}
		epsilon = epsilonInput.get();
		max = maxInput.get();
		if (max <= 0) {
			max = Integer.MAX_VALUE;
		}

		
		Map<String,Integer> treeMap = new HashMap<>();
		Map<String,Integer> clusterMap = new HashMap<>();
        FastTreeSet trees = new TreeAnnotator().new FastTreeSet(treesInput.get().getAbsolutePath(), burnInPercentageInput.get());
        trees.reset();
        double treeCount = 0;
        while (trees.hasNext()) {
            Tree tree = trees.next();
            Node root = tree.getRoot();
            root.sort();
            String s = getShortTopology(root);
            if (treeMap.containsKey(s)) {
            	treeMap.put(s, treeMap.get(s) + 1);
            } else {
            	treeMap.put(s, 1);
            }
            
            s = mapTree2Cluster(s);
            
            if (clusterMap.containsKey(s)) {
            	clusterMap.put(s, clusterMap.get(s) + 1);
            } else {
            	clusterMap.put(s, 1);
            }
            treeCount++;
        }
		
		
        String [] clusterKeys = clusterMap.keySet().toArray(new String[] {});
        Arrays.sort(clusterKeys, (a,b) -> {
        	final int i1 = clusterMap.get(a);
        	final int i2 = clusterMap.get(b);
        	if (i1 > i2) {
        		return -1;
        	} else if (i1 == i2) {
        		return 0;
        	}
        	return 1;
        });
        
        String [] treeKeys = treeMap.keySet().toArray(new String[] {});
        Arrays.sort(treeKeys, (a,b) -> {
        	final int i1 = treeMap.get(a);
        	final int i2 = treeMap.get(b);
        	if (i1 > i2) {
        		return -1;
        	} else if (i1 == i2) {
        		return 0;
        	}
        	return 1;
        });

        
        // determine size of 95% credible set
        int threshold = (int)(1 + treeCount * 95 / 100);
        int coverageSetCount = 0;
        int k = 0;
        while (coverageSetCount < threshold && k < treeKeys.length) {
        	coverageSetCount += treeMap.get(treeKeys[k++]);
        }
        
        
        // check if original tree is in 95% credible set
        boolean originalIsInCredibleSet = false;
        boolean hasOriginal = false;
        if (originalInput.get() != null && !originalInput.get().getName().equals("[[none]]")) {
        	hasOriginal = true;
            FastTreeSet originalTrees = new TreeAnnotator().new FastTreeSet(originalInput.get().getAbsolutePath(), 0);
            originalTrees.reset();
            Tree originalTree = originalTrees.next();
            originalTree.getRoot().sort();
            String originalTopology = getShortTopology(originalTree.getRoot());
            if (treeMap.containsKey(originalTopology)) {
            	if (treeKeys.length == 1) {
            		originalIsInCredibleSet = true;
            	} else {
            		int count = treeMap.get(originalTopology);
            		if (count >= treeMap.get(treeKeys[k-1])) {
            			originalIsInCredibleSet = true;
            		}
            	}
            } else {
            	originalIsInCredibleSet = false;
            }
        	
        }
        
        
		PrintStream out = System.out;
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
        	out = new PrintStream(outputInput.get());
        }
        
        out.println("The 95 precent credible set consists of " + k + " out of " + treeKeys.length + " topologies");
        if (hasOriginal) {
        	out.println("The original tree is " + (originalIsInCredibleSet ? "indeed" : "not") + " in the 95% credible set");
        }
        
        DecimalFormat f = new DecimalFormat("#.##");
        int i = 0;
        double cum = 0;
        out.println("\nsupport\t#taxa\ttopology");
        for (String s : treeKeys) {
        	if (i <= k && i < max  || printTailInput.get()) {
        		double percent = 100*treeMap.get(s)/treeCount;
        		cum += percent;
        		out.println(f.format(percent) + "%\t" + taxonCount(s) + "\t" + s.replaceAll(SPECIES_SEPARATOR, " + "));
            	if (i==k && printTailInput.get()) {
            		out.println("=== 95% credible set ===");
            	}
        	}
        	i++;
        }
        out.println(f.format(cum) +"% displayed");
        
        
        coverageSetCount = 0;
        k = 0;
        while (coverageSetCount < threshold && k < clusterKeys.length) {
        	coverageSetCount += clusterMap.get(clusterKeys[k++]);
        }


        i = 0;
        cum = 0;
        out.println("\nsupport\t#taxa\tclusters");
        for (String s : clusterKeys) {
        	if (i <= k && i < max || printTailInput.get()) {
        		double percent = 100 * clusterMap.get(s)/treeCount;
        		cum += percent;
        		out.println(f.format(percent) + "%\t" + taxonCount2(s) + "\t" + s.replaceAll(SPECIES_SEPARATOR, " + "));
            	if (i==k && printTailInput.get()) {
            		out.println("=== 95% credible set ===");
            	}
        	}
        	i++;
        }
        out.println(f.format(cum) +"% displayed");

        
        // collect info on individual clusters
        Map<String,Integer> singleClusterCount = new HashMap<>();
        for (String s : clusterKeys) {
    		int percent = clusterMap.get(s);
    		for (String singleCluster: s.split(",")) {
    			if (!singleClusterCount.containsKey(singleCluster)) {
    				singleClusterCount.put(singleCluster, percent);
    			} else {
    				singleClusterCount.put(singleCluster, percent + singleClusterCount.get(singleCluster));
    			}
    		}
        }
        
        String [] singleClusterKeys = singleClusterCount.keySet().toArray(new String[] {});
        Arrays.sort(singleClusterKeys, (a,b) -> {
        	final int i1 = singleClusterCount.get(a);
        	final int i2 = singleClusterCount.get(b);
        	if (i1 > i2) {
        		return -1;
        	} else if (i1 == i2) {
        		return 0;
        	}
        	return 1;
        });
        
        out.println("\nsupport\tcount\tcluster");
        for (String s : singleClusterKeys) {
    		double percent = 100*singleClusterCount.get(s)/treeCount;
    		out.println(f.format(percent) + "%\t" + singleClusterCount.get(s) + "\t" + s.replaceAll(SPECIES_SEPARATOR, " + "));
        }

        // collect info on pairs of taxa
        Map<String,Integer> pairClusterCount = new HashMap<>();
        for (String s : clusterKeys) {
    		int percent = clusterMap.get(s);
    		for (String singleCluster: s.split(",")) {
    			String [] taxa = singleCluster.split(SPECIES_SEPARATOR);
    			for (int d = 0; d < taxa.length; d++) {
    				for (int e = d+1; e < taxa.length; e++) {
    					String pair = taxa[d] + SPECIES_SEPARATOR + taxa[e];
    	    			if (!pairClusterCount.containsKey(pair)) {
    	    				pairClusterCount.put(pair, percent);
    	    			} else {
    	    				pairClusterCount.put(pair, percent + pairClusterCount.get(pair));
    	    			}
    				}
    			}
    			
    		}
        }
        
        String [] pairClusterKeys = pairClusterCount.keySet().toArray(new String[] {});
        Arrays.sort(pairClusterKeys, (a,b) -> {
        	final int i1 = pairClusterCount.get(a);
        	final int i2 = pairClusterCount.get(b);
        	if (i1 > i2) {
        		return -1;
        	} else if (i1 == i2) {
        		return 0;
        	}
        	return 1;
        });
        
        out.println("\nsupport\tcount\tpair");
        for (String s : pairClusterKeys) {
    		double percent = 100 * pairClusterCount.get(s)/treeCount;
    		out.println(f.format(percent) + "%\t" + pairClusterCount.get(s) + "\t" + s.replaceAll(SPECIES_SEPARATOR, " + "));
        }
       
        
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
        	out.close();
        }
	}

	
	private String mapTree2Cluster(String s) {
		s = s.replaceAll("[\\(\\),]"," ").trim();
		String [] strs = s.split(" +");
		Arrays.sort(strs);
		s = Arrays.toString(strs);
		s = s.replaceAll("[\\[\\] ]","");
		return s;
	}


	/** number of taxa = number of open braces in newick string + 1 **/
	private int taxonCount(String newick) {
		int k = 0;
		for (int i = 0; i < newick.length(); i++) {
			if (newick.charAt(i) == '(') {
				k++;
			}
		}
		return k + 1;
	}

	private int taxonCount2(String clusterString) {
		int k = 0;
		for (int i = 0; i < clusterString.length(); i++) {
			if (clusterString.charAt(i) == ',') {
				k++;
			}
		}
		return k + 1;
	}


	private String getShortTopology(Node node) {
		if (node.isLeaf()) {
			return node.getID(); 
		} else {
			if (node.getHeight() <= epsilon) {
				return collectTaxa(node);
			} else {
				return "(" + getShortTopology(node.getLeft()) + "," + 
					getShortTopology(node.getRight())  
				+")";
			}
		}
	}
	
	private String collectTaxa(Node node) {
		Set<String> taxa = new HashSet<>();
		collectTaxa(node, taxa);
		String [] taxaNames = taxa.toArray(new String[] {});
		Arrays.sort(taxaNames);
		return Arrays.toString(taxaNames).replaceAll(", ", SPECIES_SEPARATOR).replaceAll("\\[", "").replaceAll("\\]", "");
	}


	private void collectTaxa(Node node, Set<String> taxa) {
		if (node.isLeaf()) {
			taxa.add(node.getID());
		} else {
			for (Node child : node.getChildren()) {
				collectTaxa(child, taxa);
			}
		}
		
	}


	public static void main(String [] args) throws Exception {
		new Application(new ClusterTreeSetAnalyser(), "Clustered Tree Set Analyser", args);
	} // main



}

