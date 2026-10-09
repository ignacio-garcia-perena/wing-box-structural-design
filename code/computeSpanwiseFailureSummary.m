function summary = computeSpanwiseFailureSummary(stationResults, weight, divergence, analysisOptions, bayResults)
% computeSpanwiseFailureSummary Aggregate station results.
%
% Aggregates currently implemented checks. Panel buckling and combined
% combined panel ratios are included only when enabled in analysisOptions.
% Crippling and final certification margin stack-up are intentionally
% excluded.

if nargin < 4 || isempty(analysisOptions)
    analysisOptions = defaultAnalysisOptions();
end
if nargin < 5
    bayResults = [];
end

validateInputs(stationResults, weight, divergence);

nStations = numel(stationResults);

[yieldMax, yieldStationIdx] = maxStationMetric(stationResults, "yield");
[columnMax, columnStationIdx] = maxStationMetric(stationResults, "column");
[vonMisesMax, vonMisesStationIdx] = maxStationMetric(stationResults, "vonMises");
[shearMax, shearStationIdx] = maxStationMetric(stationResults, "shear");
[bendingMax, bendingStationIdx] = maxStationMetric(stationResults, "bending");

yieldEntry = stationResults(yieldStationIdx).yieldScreening.globalMaxRatio;
columnEntry = stationResults(columnStationIdx).columnBuckling.globalMaxRatio;
vonMisesEntry = stationResults(vonMisesStationIdx).vonMises.panel.globalMax;
shearEntry = stationResults(shearStationIdx).shearStress.panelShear.globalMaxAbs;
bendingEntry = stationResults(bendingStationIdx).stresses.bending.globalMaxAbs;

yieldPass = yieldMax <= 1.0;
columnPass = columnMax <= 1.0;
divergencePass = divergence.pass;
panelBucklingEnabled = isfield(analysisOptions, 'enablePanelBuckling') && ...
    analysisOptions.enablePanelBuckling;
bayAveragedPanelBucklingEnabled = panelBucklingEnabled && ...
    isfield(analysisOptions, 'panelBuckling') && ...
    isfield(analysisOptions.panelBuckling, 'useBayAveragedRecovery') && ...
    logical(analysisOptions.panelBuckling.useBayAveragedRecovery);

if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    validateBayResults(bayResults);
    [panelRcMax, panelRcBayIdx] = maxBayMetric(bayResults, "panelRc");
    [panelRbMax, panelRbBayIdx] = maxBayMetric(bayResults, "panelRb");
    [panelRsMax, panelRsBayIdx] = maxBayMetric(bayResults, "panelRs");

    panelRcEntry = bayResults(panelRcBayIdx).panelBuckling.maxRc;
    panelRbEntry = bayResults(panelRbBayIdx).panelBuckling.maxRb;
    panelRsEntry = bayResults(panelRsBayIdx).panelBuckling.maxRs;
    combinedEnabled = isPanelCombinedEnabled(analysisOptions);
    if combinedEnabled
        [panelRscMax, panelRscBayIdx] = maxBayMetric(bayResults, "panelRsc");
        [panelRbsMax, panelRbsBayIdx] = maxBayMetric(bayResults, "panelRbs");
        [panelCombinedMax, panelCombinedBayIdx] = maxBayMetric(bayResults, "panelCombined");
        panelRscEntry = bayResults(panelRscBayIdx).panelBuckling.maxR_SC;
        panelRbsEntry = bayResults(panelRbsBayIdx).panelBuckling.maxR_BS;
        panelCombinedEntry = bayResults(panelCombinedBayIdx).panelBuckling.maxCombinedActive;
        panelCombinedPass = panelCombinedMax <= 1.0;
    else
        panelRscMax = NaN;
        panelRbsMax = NaN;
        panelCombinedMax = NaN;
        panelRscEntry = makeDisabledPanelEntry("R_SC");
        panelRbsEntry = makeDisabledPanelEntry("R_BS");
        panelCombinedEntry = makeDisabledPanelEntry("R_combined_active");
        panelRscBayIdx = 1;
        panelRbsBayIdx = 1;
        panelCombinedBayIdx = 1;
        panelCombinedPass = true;
    end

    panelBucklingPass = panelRcMax <= 1.0 && panelRbMax <= 1.0 && ...
        panelRsMax <= 1.0 && panelCombinedPass;
elseif panelBucklingEnabled
    [panelRcMax, panelRcStationIdx] = maxStationMetric(stationResults, "panelRc");
    [panelRbMax, panelRbStationIdx] = maxStationMetric(stationResults, "panelRb");
    [panelRsMax, panelRsStationIdx] = maxStationMetric(stationResults, "panelRs");

    panelRcEntry = stationResults(panelRcStationIdx).panelBuckling.maxRc;
    panelRbEntry = stationResults(panelRbStationIdx).panelBuckling.maxRb;
    panelRsEntry = stationResults(panelRsStationIdx).panelBuckling.maxRs;
    combinedEnabled = isPanelCombinedEnabled(analysisOptions);
    if combinedEnabled
        [panelRscMax, panelRscStationIdx] = maxStationMetric(stationResults, "panelRsc");
        [panelRbsMax, panelRbsStationIdx] = maxStationMetric(stationResults, "panelRbs");
        [panelCombinedMax, panelCombinedStationIdx] = maxStationMetric(stationResults, "panelCombined");
        panelRscEntry = stationResults(panelRscStationIdx).panelBuckling.maxR_SC;
        panelRbsEntry = stationResults(panelRbsStationIdx).panelBuckling.maxR_BS;
        panelCombinedEntry = stationResults(panelCombinedStationIdx).panelBuckling.maxCombinedActive;
        panelCombinedPass = panelCombinedMax <= 1.0;
    else
        panelRscMax = NaN;
        panelRbsMax = NaN;
        panelCombinedMax = NaN;
        panelRscEntry = makeDisabledPanelEntry("R_SC");
        panelRbsEntry = makeDisabledPanelEntry("R_BS");
        panelCombinedEntry = makeDisabledPanelEntry("R_combined_active");
        panelRscStationIdx = 1;
        panelRbsStationIdx = 1;
        panelCombinedStationIdx = 1;
        panelCombinedPass = true;
    end

    panelBucklingPass = panelRcMax <= 1.0 && panelRbMax <= 1.0 && ...
        panelRsMax <= 1.0 && panelCombinedPass;
else
    panelRcMax = NaN;
    panelRbMax = NaN;
    panelRsMax = NaN;
    panelRcEntry = makeDisabledPanelEntry("Rc");
    panelRbEntry = makeDisabledPanelEntry("Rb");
    panelRsEntry = makeDisabledPanelEntry("Rs");
    panelRcStationIdx = 1;
    panelRbStationIdx = 1;
    panelRsStationIdx = 1;
    combinedEnabled = false;
    panelRscMax = NaN;
    panelRbsMax = NaN;
    panelCombinedMax = NaN;
    panelRscEntry = makeDisabledPanelEntry("R_SC");
    panelRbsEntry = makeDisabledPanelEntry("R_BS");
    panelCombinedEntry = makeDisabledPanelEntry("R_combined_active");
    panelRscStationIdx = 1;
    panelRbsStationIdx = 1;
    panelCombinedStationIdx = 1;
    panelCombinedPass = true;
    panelBucklingPass = true;
end

overallPass = yieldPass && columnPass && divergencePass && panelBucklingPass;

summary = struct();
summary.nStations = nStations;
summary.station_y_in = [stationResults.y_in];
summary.station_y_ft = [stationResults.y_ft];
summary.weight_lb = weight.fullspan.total_including_ribs_excluding_misc_lb;
summary.weight_excluding_ribs_misc_lb = weight.fullspan.total_excluding_ribs_and_misc_lb;
summary.weight_including_ribs_excluding_misc_lb = weight.fullspan.total_including_ribs_excluding_misc_lb;
summary.divergence_ratio = divergence.V_divergence_over_VD;

summary.max_yield_ratio = yieldMax;
summary.maxYield = addStationMetadata(yieldEntry, stationResults(yieldStationIdx));
summary.max_column_buckling_ratio = columnMax;
summary.maxColumnBuckling = addStationMetadata(columnEntry, stationResults(columnStationIdx));
summary.max_von_mises_psi = vonMisesMax;
summary.maxVonMises = addStationMetadata(vonMisesEntry, stationResults(vonMisesStationIdx));
summary.max_panel_shear_stress_psi = shearMax;
summary.maxPanelShearStress = addStationMetadata(shearEntry, stationResults(shearStationIdx));
summary.max_bending_stress_psi = bendingMax;
summary.maxBendingStress = addStationMetadata(bendingEntry, stationResults(bendingStationIdx));

summary.panelBuckling.enabled = panelBucklingEnabled;
summary.panelBuckling.use_bay_averaged_recovery = bayAveragedPanelBucklingEnabled;
if bayAveragedPanelBucklingEnabled
    summary.panelBuckling.recovery_method = string(analysisOptions.panelBuckling.recoveryMethod);
    summary.panelBuckling.n_bays = numel(bayResults);
else
    summary.panelBuckling.recovery_method = "station_endpoint_envelope";
    summary.panelBuckling.n_bays = max(0, nStations - 1);
end
summary.panelBuckling.max_Rc = panelRcMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxRc = addBayMetadata(panelRcEntry, bayResults(panelRcBayIdx));
else
    summary.panelBuckling.maxRc = addStationMetadata(panelRcEntry, stationResults(panelRcStationIdx));
end
summary.panelBuckling.max_Rb = panelRbMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxRb = addBayMetadata(panelRbEntry, bayResults(panelRbBayIdx));
else
    summary.panelBuckling.maxRb = addStationMetadata(panelRbEntry, stationResults(panelRbStationIdx));
end
summary.panelBuckling.max_Rs = panelRsMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxRs = addBayMetadata(panelRsEntry, bayResults(panelRsBayIdx));
else
    summary.panelBuckling.maxRs = addStationMetadata(panelRsEntry, stationResults(panelRsStationIdx));
end
summary.panelBuckling.combined_enabled = combinedEnabled;
summary.panelBuckling.max_R_SC = panelRscMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxR_SC = addBayMetadata(panelRscEntry, bayResults(panelRscBayIdx));
else
    summary.panelBuckling.maxR_SC = addStationMetadata(panelRscEntry, stationResults(panelRscStationIdx));
end
summary.panelBuckling.max_R_BS = panelRbsMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxR_BS = addBayMetadata(panelRbsEntry, bayResults(panelRbsBayIdx));
else
    summary.panelBuckling.maxR_BS = addStationMetadata(panelRbsEntry, stationResults(panelRbsStationIdx));
end
summary.panelBuckling.max_combined_active = panelCombinedMax;
if panelBucklingEnabled && bayAveragedPanelBucklingEnabled
    summary.panelBuckling.maxCombinedActive = addBayMetadata(panelCombinedEntry, ...
        bayResults(panelCombinedBayIdx));
else
    summary.panelBuckling.maxCombinedActive = addStationMetadata(panelCombinedEntry, ...
        stationResults(panelCombinedStationIdx));
end
summary.panelBuckling.combined_pass = panelCombinedPass;
summary.panelBuckling.pass = panelBucklingPass;
if panelBucklingEnabled && isfield(analysisOptions, 'panelBuckling')
    summary.panelBuckling.coefficients = analysisOptions.panelBuckling;
else
    summary.panelBuckling.coefficients = struct();
end

summary.yield_pass = yieldPass;
summary.column_buckling_pass = columnPass;
summary.divergence_pass = divergencePass;
summary.panel_buckling_pass = panelBucklingPass;
summary.overall_pass = overallPass;
summary.method = ...
    "Maxima aggregated over rib stations, load cases, and existing panel/boom entries.";
if bayAveragedPanelBucklingEnabled
    summary.method = ...
        "Yield/column/von Mises maxima aggregated over rib stations; panel buckling final pass/fail aggregated over rib bays using bay-averaged recovery.";
end
if panelBucklingEnabled
    if combinedEnabled
        summary.exclusions = [
            "crippling"
            "general crippling interaction margin stack-up"
            "final margins"
            "optimisation"
            ];
    else
        summary.exclusions = [
            "Combined panel buckling ratios"
            "crippling"
            "general crippling interaction margin stack-up"
            "final margins"
            "optimisation"
            ];
    end
else
    summary.exclusions = [
        "panel buckling"
        "crippling"
        "Combined panel buckling ratios"
        "general crippling interaction margin stack-up"
        "final margins"
        "optimisation"
        ];
end
end

function validateInputs(stationResults, weight, divergence)
if nargin < 1 || isempty(stationResults) || ~isstruct(stationResults)
    error('computeSpanwiseFailureSummary:MissingStationResults', ...
        'stationResults struct array is required.');
end
requiredStationFields = ["y_in", "y_ft", "stresses", "shearStress", ...
    "vonMises", "yieldScreening", "columnBuckling"];
missing = setdiff(requiredStationFields, string(fieldnames(stationResults)));
if ~isempty(missing)
    error('computeSpanwiseFailureSummary:MissingStationFields', ...
        'stationResults is missing fields: %s', strjoin(missing, ', '));
end
if nargin < 2 || isempty(weight) || ~isfield(weight, 'fullspan') || ...
        ~isfield(weight.fullspan, 'total_including_ribs_excluding_misc_lb')
    error('computeSpanwiseFailureSummary:InvalidWeight', ...
        'weight.fullspan.total_including_ribs_excluding_misc_lb is required.');
end
if nargin < 3 || isempty(divergence) || ...
        ~isfield(divergence, 'V_divergence_over_VD') || ...
        ~isfield(divergence, 'pass')
    error('computeSpanwiseFailureSummary:InvalidDivergence', ...
        'divergence ratio and pass flag are required.');
end
end

function [value, stationIdx] = maxStationMetric(stationResults, metricName)
values = zeros(numel(stationResults), 1);
for iStation = 1:numel(stationResults)
    switch string(metricName)
        case "yield"
            values(iStation) = stationResults(iStation).yieldScreening.globalMaxRatio.yield_ratio;
        case "column"
            values(iStation) = stationResults(iStation).columnBuckling.globalMaxRatio.buckling_ratio;
        case "vonMises"
            values(iStation) = stationResults(iStation).vonMises.panel.globalMax.sigma_vm_psi;
        case "shear"
            values(iStation) = stationResults(iStation).shearStress.panelShear.globalMaxAbs.abs_tau_psi;
        case "bending"
            values(iStation) = stationResults(iStation).stresses.bending.globalMaxAbs.abs_sigma_psi;
        case "panelRc"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxRc.ratio;
        case "panelRb"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxRb.ratio;
        case "panelRs"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxRs.ratio;
        case "panelRsc"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxR_SC.ratio;
        case "panelRbs"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxR_BS.ratio;
        case "panelCombined"
            requirePanelBucklingEnabled(stationResults(iStation), metricName);
            values(iStation) = stationResults(iStation).panelBuckling.maxCombinedActive.ratio;
        otherwise
            error('computeSpanwiseFailureSummary:UnknownMetric', ...
                'Unknown metric %s.', metricName);
    end
end
if any(~isfinite(values))
    error('computeSpanwiseFailureSummary:InvalidMetricValues', ...
        'Metric %s contains non-finite values.', metricName);
end
[value, stationIdx] = max(values);
end

function validateBayResults(bayResults)
if isempty(bayResults) || ~isstruct(bayResults)
    error('computeSpanwiseFailureSummary:MissingBayResults', ...
        'Bay-averaged panel buckling requires non-empty bayResults.');
end
requiredBayFields = ["bay_index", "y_start_in", "y_end_in", "y_mid_in", ...
    "bay_length_in", "panelBuckling"];
missing = setdiff(requiredBayFields, string(fieldnames(bayResults)));
if ~isempty(missing)
    error('computeSpanwiseFailureSummary:MissingBayFields', ...
        'bayResults is missing fields: %s', strjoin(missing, ', '));
end
for iBay = 1:numel(bayResults)
    if ~isfield(bayResults(iBay).panelBuckling, 'enabled') || ...
            ~bayResults(iBay).panelBuckling.enabled
        error('computeSpanwiseFailureSummary:BayPanelBucklingDisabled', ...
            'Bay %d does not contain enabled panel buckling results.', iBay);
    end
end
end

function [value, bayIdx] = maxBayMetric(bayResults, metricName)
values = zeros(numel(bayResults), 1);
for iBay = 1:numel(bayResults)
    pb = bayResults(iBay).panelBuckling;
    switch string(metricName)
        case "panelRc"
            values(iBay) = pb.maxRc.ratio;
        case "panelRb"
            values(iBay) = pb.maxRb.ratio;
        case "panelRs"
            values(iBay) = pb.maxRs.ratio;
        case "panelRsc"
            values(iBay) = pb.maxR_SC.ratio;
        case "panelRbs"
            values(iBay) = pb.maxR_BS.ratio;
        case "panelCombined"
            values(iBay) = pb.maxCombinedActive.ratio;
        otherwise
            error('computeSpanwiseFailureSummary:UnknownBayMetric', ...
                'Unknown bay metric %s.', metricName);
    end
end
if any(~isfinite(values))
    error('computeSpanwiseFailureSummary:InvalidBayMetricValues', ...
        'Bay metric %s contains non-finite values.', metricName);
end
[value, bayIdx] = max(values);
end

function enabled = isPanelCombinedEnabled(analysisOptions)
enabled = isfield(analysisOptions, 'panelBuckling') && ...
    isfield(analysisOptions.panelBuckling, 'enableCombinedRatios') && ...
    logical(analysisOptions.panelBuckling.enableCombinedRatios);
end

function entry = addStationMetadata(entry, stationResult)
entry.station_index = stationResult.station_index;
entry.station_y_in = stationResult.y_in;
entry.station_y_ft = stationResult.y_ft;
end

function entry = addBayMetadata(entry, bayResult)
entry.bay_index = bayResult.bay_index;
entry.bay_start_y_in = bayResult.y_start_in;
entry.bay_end_y_in = bayResult.y_end_in;
entry.bay_mid_y_in = bayResult.y_mid_in;
entry.bay_length_in = bayResult.bay_length_in;
end

function requirePanelBucklingEnabled(stationResult, metricName)
if ~isfield(stationResult, 'panelBuckling') || ...
        ~isfield(stationResult.panelBuckling, 'enabled') || ...
        ~stationResult.panelBuckling.enabled
    error('computeSpanwiseFailureSummary:PanelBucklingDisabled', ...
        'Metric %s requires enabled panel buckling results.', metricName);
end
end

function entry = makeDisabledPanelEntry(ratioName)
entry = struct();
entry.ratio_name = string(ratioName);
entry.ratio = NaN;
entry.case_index = NaN;
entry.case_name = "";
entry.panel_index = NaN;
entry.panel_name = "";
entry.panel_group = "";
entry.applied_stress_psi = NaN;
entry.critical_stress_psi = NaN;
entry.a_in = NaN;
entry.b_in = NaN;
entry.t_in = NaN;
entry.aspect_ratio_a_over_b = NaN;
end
