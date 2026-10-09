function analysisOptions = defaultFinalAnalysisOptions()
% defaultFinalAnalysisOptions Return final-design analysis options.

analysisOptions = defaultAnalysisOptions();
analysisOptions.loadFactor = 1.5;
analysisOptions.loadFactorPurpose = "ultimate_load_final_design";
analysisOptions.usePanelIdealisation = true;
analysisOptions.panelIdealisationMethod = "modified_megson_z_ratio";
analysisOptions.enablePanelBuckling = true;
analysisOptions.panelBuckling.enableCombinedRatios = true;
analysisOptions.panelBuckling.useBayAveragedRecovery = true;
analysisOptions.panelBuckling.enableCurvedPanelCorrection = true;
analysisOptions.suppressPracticalGeometryWarnings = true;
analysisOptions.finalWorkflow = true;
analysisOptions.notes(end + 1, 1) = ...
    "defaultFinalAnalysisOptions enables the final-design analysis settings.";
analysisOptions.notes(end + 1, 1) = ...
    "Final Monte Carlo suppresses broad practical spar-location warnings while retaining hard geometry validation.";
end
