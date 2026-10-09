function analysisOptions = defaultAnalysisOptions()
% defaultAnalysisOptions Return default project analysis options.
%
% Outputs:
%   analysisOptions - Struct controlling analysis-level load scaling.
%
% Notes:
%   loadFactor = 1.0 preserves the limit-load reference comparison.
%   loadFactor = 1.5 should be used for final-design ultimate-load
%   structural screening.

analysisOptions = struct();
analysisOptions.loadFactor = 1.0;
analysisOptions.loadFactorPurpose = "limit_load_reference";
analysisOptions.applyLoadFactorTo = ["Vx", "Vz", "Mx", "Mz", "My"];
analysisOptions.suppressPracticalGeometryWarnings = false;
analysisOptions.usePanelIdealisation = true;
analysisOptions.panelIdealisationMethod = "modified_megson_z_ratio";
analysisOptions.enablePanelBuckling = false;
analysisOptions.panelBuckling.Kc = 6;
analysisOptions.panelBuckling.Kb = 32;
analysisOptions.panelBuckling.Ks = 8;
analysisOptions.panelBuckling.method = "constant_high_aspect_ratio_coefficients";
analysisOptions.panelBuckling.assumption = ...
    "Representative high-aspect-ratio coefficients between simply supported and clamped edge restraint; approved for project-level screening.";
analysisOptions.panelBuckling.enableCombinedRatios = true;
analysisOptions.panelBuckling.combinedMethod = "combined_loads_interaction";
analysisOptions.panelBuckling.combinedNotes = ...
    "S+C uses Rs^1.5 + Rc < 1; B+S uses Rb^2 + Rs^2 < 1.";
analysisOptions.panelBuckling.useBayAveragedRecovery = true;
analysisOptions.panelBuckling.recoveryMethod = ...
    "bay_averaged_four_stress_and_shear_flow";
analysisOptions.panelBuckling.enableCurvedPanelCorrection = true;
analysisOptions.panelBuckling.curvedPanelCorrectionMethod = ...
    "cylindrical_compression_addition";
analysisOptions.panelBuckling.curvedPanelApplyTo = "compression_only";
analysisOptions.panelBuckling.notes = ...
    "Panel buckling uses bay-averaged compressive stress from four surrounding boom/stringer stresses and average shear flow between rib stations.";
analysisOptions.notes = [
    "loadFactor = 1.0 preserves the limit-load reference comparison."
    "loadFactor = 1.5 should be used for final-design ultimate-load structural screening."
    "The load factor applies to structural loads used for stress, shear-flow, von Mises, yield, and column buckling checks."
    "Weight and divergence are geometry/material checks and are not multiplied by this structural load factor."
    "suppressPracticalGeometryWarnings = false leaves broad practical geometry warnings visible for manual checks."
    "usePanelIdealisation = false preserves the boom-only section model."
    "usePanelIdealisation = true assigns panel/web area into adjacent boom areas for the final-design section model."
    "Modified Megson panel idealisation is the final-design panel idealisation method; equal-split panel area remains available for regression comparison."
    "enablePanelBuckling = false preserves the G3 spanwise checks without panel buckling."
    "enablePanelBuckling = true adds separate Rc, Rb, and Rs panel buckling screens."
    "Panel buckling uses constant project-approved high-aspect-ratio coefficients unless overridden."
    "Panel buckling combined ratios use S+C and B+S rules when enabled."
    "Panel buckling pass/fail uses bay-averaged stress recovery when enabled."
    "Curved airfoil skin/nose panels use the cylindrical compression addition when enableCurvedPanelCorrection = true."
    ];
end
