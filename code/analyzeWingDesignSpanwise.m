function spanwiseResult = analyzeWingDesignSpanwise(design, loads, designOptions, analysisOptions)
% analyzeWingDesignSpanwise Evaluate current checks at each rib station.
%
% This spanwise analysis reuses the section stress modules by passing a
% local root-shaped load struct for each rib station. The section is
% constant along the span for the current optimisation model.

if nargin < 4 || isempty(analysisOptions)
    analysisOptions = defaultAnalysisOptions();
end

validateInputs(design, loads, designOptions, analysisOptions);

wing = setGeometry(design, designOptions, analysisOptions);
section = computeSectionUnitLoads(wing, analysisOptions);
weightSection = sectionForPhysicalWeight(section);
weight = computeWeight(wing, weightSection);
divergence = computeDivergence(wing, section);
spanLoads = interpolateSpanwiseLoads(loads, wing, analysisOptions);

stationResults = repmat(makeStationResult(), spanLoads.nStations, 1);
for iStation = 1:spanLoads.nStations
    localLoads = makeLocalLoads(loads, spanLoads, iStation);

    stresses = computeStresses(localLoads, section);
    shearFlows = computeShearFlows(localLoads, section);
    shearStress = computePanelShearStress(shearFlows);
    vonMises = computeVonMises(stresses, shearStress);
    yieldScreening = computeYieldScreening(vonMises, wing);
    columnBuckling = computeColumnBuckling(wing, section, stresses, ...
        wing.ribs.spacing_in, wing.columnBuckling.K_column);
    if isfield(analysisOptions, 'enablePanelBuckling') && analysisOptions.enablePanelBuckling
        panelBuckling = computePanelBuckling(wing, shearStress, vonMises, analysisOptions);
    else
        panelBuckling = makeDisabledPanelBuckling();
    end

    stationResults(iStation).station_index = iStation;
    stationResults(iStation).y_in = spanLoads.stations_in(iStation);
    stationResults(iStation).y_ft = spanLoads.stations_ft(iStation);
    stationResults(iStation).loads = localLoads;
    stationResults(iStation).stresses = stresses;
    stationResults(iStation).shearFlows = shearFlows;
    stationResults(iStation).shearStress = shearStress;
    stationResults(iStation).vonMises = vonMises;
    stationResults(iStation).yieldScreening = yieldScreening;
    stationResults(iStation).columnBuckling = columnBuckling;
    stationResults(iStation).panelBuckling = panelBuckling;
end

bayResults = makeBayResults(stationResults, wing, section, analysisOptions);
summary = computeSpanwiseFailureSummary(stationResults, weight, divergence, analysisOptions, bayResults);

spanwiseResult = struct();
spanwiseResult.design = design;
spanwiseResult.wing = wing;
spanwiseResult.section = section;
spanwiseResult.weight = weight;
spanwiseResult.divergence = divergence;
spanwiseResult.spanLoads = spanLoads;
spanwiseResult.stationResults = stationResults;
spanwiseResult.bayResults = bayResults;
spanwiseResult.summary = summary;
spanwiseResult.analysisOptions = analysisOptions;
spanwiseResult.status = "spanwise_rib_station_analysis";
if isfield(analysisOptions, 'enablePanelBuckling') && analysisOptions.enablePanelBuckling
    spanwiseResult.status = "spanwise_analysis_with_panel_buckling_screening";
    if isfield(analysisOptions, 'panelBuckling') && ...
            isfield(analysisOptions.panelBuckling, 'enableCombinedRatios') && ...
            analysisOptions.panelBuckling.enableCombinedRatios
        spanwiseResult.status = ...
            "spanwise_analysis_with_panel_buckling_and_combined_panel_ratios";
    end
    if useBayAveragedPanelBuckling(analysisOptions)
        spanwiseResult.status = ...
            "spanwise_analysis_with_bay_averaged_panel_buckling";
    end
    spanwiseResult.exclusions = [
        "crippling"
        "general crippling interaction margin stack-up"
        "final margins"
        "optimisation"
        "spanwise-varying section properties"
        ];
else
    spanwiseResult.exclusions = [
        "panel buckling"
        "crippling"
        "combined panel ratios"
        "general crippling interaction margin stack-up"
        "final margins"
        "optimisation"
        "spanwise-varying section properties"
        ];
end
end

function validateInputs(design, loads, designOptions, analysisOptions)
if nargin < 1 || isempty(design) || ~isstruct(design)
    error('analyzeWingDesignSpanwise:MissingDesign', ...
        'design struct is required.');
end
if nargin < 2 || isempty(loads) || ~isstruct(loads)
    error('analyzeWingDesignSpanwise:MissingLoads', ...
        'loads struct from criticalLoads is required.');
end
if nargin < 3 || isempty(designOptions) || ~isstruct(designOptions)
    error('analyzeWingDesignSpanwise:MissingDesignOptions', ...
        'designOptions struct from loadDesignOptions is required.');
end
if nargin < 4 || isempty(analysisOptions) || ~isstruct(analysisOptions)
    error('analyzeWingDesignSpanwise:MissingAnalysisOptions', ...
        'analysisOptions struct is required.');
end
if ~isfield(analysisOptions, 'loadFactor') || ...
        ~(isnumeric(analysisOptions.loadFactor) && isscalar(analysisOptions.loadFactor) && ...
        isfinite(analysisOptions.loadFactor) && analysisOptions.loadFactor > 0)
    error('analyzeWingDesignSpanwise:InvalidLoadFactor', ...
        'analysisOptions.loadFactor must be a positive finite scalar.');
end
end

function localLoads = makeLocalLoads(loads, spanLoads, iStation)
localLoads = struct();
localLoads.caseNames = spanLoads.caseNames;
localLoads.nCases = spanLoads.nCases;
localLoads.root.Vx_lb = spanLoads.Vx_lb(:, iStation);
localLoads.root.Vz_lb = spanLoads.Vz_lb(:, iStation);
localLoads.root.Mx_lbin = spanLoads.Mx_lbin(:, iStation);
localLoads.root.Mz_lbin = spanLoads.Mz_lbin(:, iStation);
localLoads.root.My_lbin = spanLoads.My_lbin(:, iStation);
localLoads.root.Mx_lbft = localLoads.root.Mx_lbin / 12.0;
localLoads.root.Mz_lbft = localLoads.root.Mz_lbin / 12.0;
localLoads.root.My_lbft = localLoads.root.My_lbin / 12.0;
localLoads.units.force = "lb";
localLoads.units.length = "in";
localLoads.units.moment = "lb-in";
localLoads.station.y_in = spanLoads.stations_in(iStation);
localLoads.station.y_ft = spanLoads.stations_ft(iStation);
localLoads.station.index = iStation;
localLoads.loadFactorApplied = spanLoads.loadFactorApplied;
localLoads.loadFactorPurpose = spanLoads.loadFactorPurpose;
localLoads.source = "interpolateSpanwiseLoads rib-station values";
if isfield(loads, 'aero')
    localLoads.aero = loads.aero;
end
if isfield(loads, 'notes')
    localLoads.notes = loads.notes;
end
end

function result = makeStationResult()
result = struct();
result.station_index = NaN;
result.y_in = NaN;
result.y_ft = NaN;
result.loads = struct();
result.stresses = struct();
result.shearFlows = struct();
result.shearStress = struct();
result.vonMises = struct();
result.yieldScreening = struct();
result.columnBuckling = struct();
result.panelBuckling = struct();
end

function panelBuckling = makeDisabledPanelBuckling()
panelBuckling = struct();
panelBuckling.enabled = false;
panelBuckling.status = "disabled";
panelBuckling.pass = true;
panelBuckling.note = "Set analysisOptions.enablePanelBuckling = true to compute panel buckling screens.";
end

function bayResults = makeBayResults(stationResults, wing, section, analysisOptions)
nBays = max(0, numel(stationResults) - 1);
bayResults = repmat(makeBayResult(), nBays, 1);
if nBays == 0
    return;
end
for iBay = 1:nBays
    bayResults(iBay).bay_index = iBay;
    bayResults(iBay).station_start_index = stationResults(iBay).station_index;
    bayResults(iBay).station_end_index = stationResults(iBay + 1).station_index;
    bayResults(iBay).y_start_in = stationResults(iBay).y_in;
    bayResults(iBay).y_end_in = stationResults(iBay + 1).y_in;
    bayResults(iBay).y_mid_in = 0.5 * (bayResults(iBay).y_start_in + bayResults(iBay).y_end_in);
    bayResults(iBay).bay_length_in = bayResults(iBay).y_end_in - bayResults(iBay).y_start_in;
    if useBayAveragedPanelBuckling(analysisOptions)
        bayResults(iBay).panelBuckling = computePanelBucklingBayAveraged( ...
            wing, section, stationResults(iBay), stationResults(iBay + 1), analysisOptions);
    else
        bayResults(iBay).panelBuckling = makeDisabledBayPanelBuckling();
    end
end
end

function tf = useBayAveragedPanelBuckling(analysisOptions)
tf = isfield(analysisOptions, 'enablePanelBuckling') && analysisOptions.enablePanelBuckling && ...
    isfield(analysisOptions, 'panelBuckling') && ...
    isfield(analysisOptions.panelBuckling, 'useBayAveragedRecovery') && ...
    logical(analysisOptions.panelBuckling.useBayAveragedRecovery);
end

function result = makeBayResult()
result = struct();
result.bay_index = NaN;
result.station_start_index = NaN;
result.station_end_index = NaN;
result.y_start_in = NaN;
result.y_end_in = NaN;
result.y_mid_in = NaN;
result.bay_length_in = NaN;
result.panelBuckling = struct();
end

function panelBuckling = makeDisabledBayPanelBuckling()
panelBuckling = struct();
panelBuckling.enabled = false;
panelBuckling.status = "disabled";
panelBuckling.pass = true;
panelBuckling.note = "Bay-averaged panel buckling recovery is disabled.";
end

function sectionOut = sectionForPhysicalWeight(sectionIn)
sectionOut = sectionIn;
if ~isfield(sectionIn, 'idealization') || ~isfield(sectionIn.idealization, 'enabled') || ...
        ~sectionIn.idealization.enabled
    return;
end

sectionOut.sparCaps = resetBoomArrayToOriginalArea(sectionOut.sparCaps);
sectionOut.stringers.top = resetBoomArrayToOriginalArea(sectionOut.stringers.top);
sectionOut.stringers.bottom = resetBoomArrayToOriginalArea(sectionOut.stringers.bottom);
sectionOut.stringers.noseTop = resetBoomArrayToOriginalArea(sectionOut.stringers.noseTop);
sectionOut.stringers.noseBottom = resetBoomArrayToOriginalArea(sectionOut.stringers.noseBottom);
end

function booms = resetBoomArrayToOriginalArea(booms)
for iBoom = 1:numel(booms)
    if isfield(booms(iBoom), 'original_area_in2') && ...
            isfinite(booms(iBoom).original_area_in2) && booms(iBoom).original_area_in2 > 0
        booms(iBoom).area_in2 = booms(iBoom).original_area_in2;
    end
end
end
