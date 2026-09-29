open module speedemon {
	requires beast.pkgmgmt;
    requires beast.base;
    requires static javafx.controls;
    requires static beast.fx;
    requires static biceps;
    requires static org.apache.commons.statistics.distribution;
    requires static org.apache.commons.math4.legacy;
    

    exports speedemon;
    exports speedemon.inputeditor;


    provides beast.base.core.BEASTInterface with
        
        speedemon.BirthDeathSkylineCollapseModel,
        speedemon.BirthDeathSkylineModel,
        speedemon.ClusterCounter,
        speedemon.ClusterOperator,
        speedemon.ClusterTreeSetAnalyser,
        speedemon.CollapseModel,
        speedemon.TreeAboveThreshold,
        speedemon.UniformThresholdOperator,
        speedemon.YuleSkylineCollapse;
    
    
    provides beastfx.app.inputeditor.InputEditor with
    	speedemon.inputeditor.ConstantInputEditor;
 

}
