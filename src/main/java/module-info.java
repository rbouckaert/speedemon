open module speedemon {
	requires beast.pkgmgmt;
    requires beast.base;
    requires static javafx.controls;
    requires static beast.fx;
    requires static org.apache.commons.statistics.distribution;
    requires static org.apache.commons.math4.legacy;
    

    exports speedemon;


    provides beast.base.core.BEASTInterface with
        
        speedemon.BirthDeathSkylineCollapseModel,
        speedemon.BirthDeathSkylineModel,
        speedemon.ClusterCounter,
        speedemon.ClusterOperator,
        speedemon.ClusterTreeSetAnalyser,
        speedemon.CollapseModel,
        speedemon.TreeAboveThreshold,
        speedemon.UniformThresholdOperator,
        speedemon.EpochTreeDistribution,
        speedemon.YuleSkyline,
        speedemon.YuleSkylineCollapse;
    
    

}
