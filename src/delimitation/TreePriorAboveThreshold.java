package delimitation;

// tree prior that assumes the tree is cut off at threshold epsilon
// and only calculates the contribution above that threshold
public interface TreePriorAboveThreshold {
	
	/**
	 * like calculateLogP, but only considering the part of the tree above threshold epsilon
	 * @param epsilon: threshold value
	 * @return
	 */
	double calculateLogPabove(double epsilon);

}
