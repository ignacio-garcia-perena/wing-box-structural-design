function adaptiveResults = runAdaptiveFinalMonteCarlo(loads, designOptions, analysisOptions, optimisationOptions)
% runAdaptiveFinalMonteCarlo Run staged adaptive Monte Carlo optimisation search.
%
% Each stage uses the same final structural analysis as
% runFinalMonteCarlo. Passing designs from a stage are sorted by weight and
% the elite subset refines candidate lists for the next stage.

validateInputs(loads, designOptions, analysisOptions, optimisationOptions);

initialRanges = optimisationOptions.ranges;
currentRanges = initialRanges;
numStages = optimisationOptions.adaptive.numStages;
nPerStage = optimisationOptions.adaptive.nPerStage;

allDesigns = table();
bestByStage = table();
stageRows = repmat(makeStageSummaryRow(), numStages, 1);
rangeHistoryRows = repmat(makeRangeHistoryRow(), 0, 1);
stageResults = cell(numStages, 1);
completedStages = 0;
stoppedEarly = false;
stopReason = "";

fprintf('\nAdaptive final Monte Carlo is ON\n');
fprintf('Stages: %d\n', numStages);
fprintf('N_PER_STAGE: %s\n', mat2str(nPerStage));
fprintf('Elite fraction: %.3f, minimum elite count: %d\n\n', ...
    optimisationOptions.adaptive.eliteFraction, ...
    optimisationOptions.adaptive.minEliteCount);

for iStage = 1:numStages
    stageOptions = optimisationOptions;
    stageOptions.ranges = currentRanges;
    stageOptions.numDesigns = nPerStage(iStage);

    fprintf('\nAdaptive Monte Carlo stage %d/%d\n', iStage, numStages);
    fprintf('N = %d\n', stageOptions.numDesigns);
    printRangeSummary(stageOptions.ranges, "Active ranges:");

    stageResult = runFinalMonteCarlo(loads, designOptions, analysisOptions, stageOptions);
    completedStages = iStage;
    stageResults{iStage} = compactStageResult(stageResult);

    stageAll = addvars(stageResult.allDesigns, ...
        repmat(iStage, height(stageResult.allDesigns), 1), ...
        'Before', 1, 'NewVariableNames', 'stage');
    allDesigns = appendTables(allDesigns, stageAll);

    nPassing = height(stageResult.passingDesigns);
    if nPassing > 0
        stageBest = addvars(stageResult.bestDesign, iStage, ...
            'Before', 1, 'NewVariableNames', 'stage');
        bestByStage = appendTables(bestByStage, stageBest);
    end

    nElite = eliteCount(nPassing, optimisationOptions.adaptive);
    stageRows(iStage) = makeStageSummaryRowFromResult(iStage, ...
        stageOptions.numDesigns, stageResult, nElite);

    fprintf('Stage %d summary: passing = %d/%d (%.2f%%), best weight = %.4f lb, elite count = %d\n', ...
        iStage, nPassing, height(stageResult.allDesigns), ...
        100 * safeDivide(nPassing, height(stageResult.allDesigns)), ...
        stageRows(iStage).best_weight_including_ribs_excluding_misc_lb, nElite);

    if iStage < numStages
        if nPassing == 0
            stoppedEarly = true;
            stopReason = "No passing designs in stage " + string(iStage) + ...
                "; adaptive search stopped and active final candidate artifacts were not overwritten.";
            fprintf('\n%s\n', stopReason);
            break;
        end
        eliteDesigns = stageResult.passingDesigns(1:nElite, :);
        refineOptions = struct();
        refineOptions.rangePaddingSteps = optimisationOptions.adaptive.rangePaddingSteps;
        refineOptions.minSectionCandidates = optimisationOptions.adaptive.minSectionCandidates;
        [currentRanges, newHistory] = refineFinalMonteCarloRanges( ...
            currentRanges, eliteDesigns, initialRanges, refineOptions, iStage);
        rangeHistoryRows = [rangeHistoryRows; newHistory(:)]; %#ok<AGROW>
        printRangeSummary(currentRanges, "Refined ranges for Stage " + string(iStage + 1) + ":");
    end
end

stageSummary = struct2table(stageRows(1:completedStages));
if isempty(rangeHistoryRows)
    rangeHistory = emptyRangeHistoryTable();
else
    rangeHistory = struct2table(rangeHistoryRows);
end

passingAll = allDesigns(allDesigns.overall_pass == true & ...
    allDesigns.error_flag == false, :);
if height(passingAll) > 0
    passingAll = sortrows(passingAll, ...
        'weight_including_ribs_excluding_misc_lb', 'ascend');
    bestOverallDesign = passingAll(1, :);
else
    bestOverallDesign = passingAll;
end

finalStageResult = stageResults{completedStages};
finalStageHasPassing = ~isempty(finalStageResult) && ...
    height(finalStageResult.passingDesigns) > 0 && ~stoppedEarly;
if finalStageHasPassing
    stage4BestDesign = finalStageResult.bestDesign;
else
    stage4BestDesign = table();
end
globalHasPassing = height(bestOverallDesign) > 0;

adaptiveResults = struct();
adaptiveResults.allDesigns = allDesigns;
adaptiveResults.stageSummary = stageSummary;
adaptiveResults.bestByStage = bestByStage;
adaptiveResults.rangeHistory = rangeHistory;
adaptiveResults.stageResults = stageResults(1:completedStages);
adaptiveResults.initialRanges = initialRanges;
adaptiveResults.finalRanges = currentRanges;
adaptiveResults.analysisOptions = analysisOptions;
adaptiveResults.optimisationOptions = optimisationOptions;
adaptiveResults.completedStages = completedStages;
adaptiveResults.stoppedEarly = stoppedEarly;
adaptiveResults.stopReason = stopReason;
adaptiveResults.bestOverallDesign = bestOverallDesign;
adaptiveResults.finalStageResult = finalStageResult;
adaptiveResults.stage4BestDesign = stage4BestDesign;
adaptiveResults.bestDesign = table();
adaptiveResults.bestDesignStruct = [];
adaptiveResults.bestSpanwiseResult = [];
adaptiveResults.candidateUsedForUpdate = "none";

if globalHasPassing
    adaptiveResults.bestDesign = bestOverallDesign;
    adaptiveResults.bestDesignStruct = designStructForBestRow(bestOverallDesign, stageResults);
    adaptiveResults.candidateUsedForUpdate = "global_best_across_all_stages";
    fprintf('\nRe-running detailed spanwise analysis for global best passing design across all stages...\n');
    adaptiveResults.bestSpanwiseResult = analyzeWingDesignSpanwise( ...
        adaptiveResults.bestDesignStruct, loads, designOptions, analysisOptions);
end

adaptiveResults.summary = makeAdaptiveSummary(adaptiveResults);

if globalHasPassing
    updateOnlyIfImproved = getFieldWithDefault(optimisationOptions, ...
        'updateFinalOnlyIfImproved', false);
    adaptiveResults.finalCandidateUpdate = updateFinalCandidateArtifactsIfImproved( ...
        adaptiveResults, analysisOptions, updateOnlyIfImproved);
else
    fprintf('\nNo final-stage passing design is available; active final candidate artifacts were not overwritten.\n');
    updateOnlyIfImproved = getFieldWithDefault(optimisationOptions, ...
        'updateFinalOnlyIfImproved', false);
    adaptiveResults.finalCandidateUpdate = updateFinalCandidateArtifactsIfImproved( ...
        adaptiveResults, analysisOptions, updateOnlyIfImproved);
end
end

function design = designStructForBestRow(bestRow, stageResults)
stageNumber = bestRow.stage(1);
designIndex = bestRow.index(1);
if ~isfinite(stageNumber) || ~isfinite(designIndex) || ...
        stageNumber < 1 || stageNumber > numel(stageResults)
    error('runAdaptiveFinalMonteCarlo:InvalidGlobalBestRow', ...
        'Global best row does not identify a valid stage/design index.');
end
stageResult = stageResults{stageNumber};
if isempty(stageResult) || ~isfield(stageResult, 'designs') || ...
        designIndex < 1 || designIndex > numel(stageResult.designs) || ...
        isempty(stageResult.designs{designIndex})
    error('runAdaptiveFinalMonteCarlo:MissingGlobalBestDesign', ...
        'Could not recover design struct for global best row.');
end
design = stageResult.designs{designIndex};
end

function compact = compactStageResult(stageResult)
compact = struct();
compact.allDesigns = stageResult.allDesigns;
compact.passingDesigns = stageResult.passingDesigns;
compact.bestDesign = stageResult.bestDesign;
compact.bestDesignStruct = stageResult.bestDesignStruct;
compact.failureCounts = stageResult.failureCounts;
compact.summary = stageResult.summary;
compact.rngState = stageResult.rngState;
compact.designs = stageResult.designs;
end

function validateInputs(loads, designOptions, analysisOptions, optimisationOptions)
if nargin < 1 || isempty(loads) || ~isstruct(loads)
    error('runAdaptiveFinalMonteCarlo:MissingLoads', 'loads struct is required.');
end
if nargin < 2 || isempty(designOptions) || ~isstruct(designOptions)
    error('runAdaptiveFinalMonteCarlo:MissingDesignOptions', ...
        'designOptions struct is required.');
end
if nargin < 3 || isempty(analysisOptions) || ~isstruct(analysisOptions)
    error('runAdaptiveFinalMonteCarlo:MissingAnalysisOptions', ...
        'analysisOptions struct is required.');
end
if nargin < 4 || isempty(optimisationOptions) || ~isstruct(optimisationOptions)
    error('runAdaptiveFinalMonteCarlo:MissingOptimisationOptions', ...
        'optimisationOptions struct is required.');
end
if ~isfield(optimisationOptions, 'adaptive')
    error('runAdaptiveFinalMonteCarlo:MissingAdaptiveOptions', ...
        'optimisationOptions.adaptive is required.');
end
required = ["nPerStage", "numStages", "eliteFraction", "minEliteCount", ...
    "rangePaddingSteps", "minSectionCandidates"];
missing = setdiff(required, string(fieldnames(optimisationOptions.adaptive)));
if ~isempty(missing)
    error('runAdaptiveFinalMonteCarlo:MissingAdaptiveFields', ...
        'optimisationOptions.adaptive is missing fields: %s', strjoin(missing, ', '));
end
nPerStage = optimisationOptions.adaptive.nPerStage;
if numel(nPerStage) ~= optimisationOptions.adaptive.numStages
    error('runAdaptiveFinalMonteCarlo:StageCountMismatch', ...
        'numStages must equal numel(nPerStage).');
end
if any(~isfinite(nPerStage)) || any(nPerStage <= 0) || any(mod(nPerStage, 1) ~= 0)
    error('runAdaptiveFinalMonteCarlo:InvalidNPerStage', ...
        'N per stage values must be positive integers.');
end
end

function nElite = eliteCount(nPassing, adaptiveOptions)
if nPassing <= 0
    nElite = 0;
else
    nElite = max(adaptiveOptions.minEliteCount, ...
        ceil(adaptiveOptions.eliteFraction * nPassing));
    nElite = min(nElite, nPassing);
end
end

function row = makeStageSummaryRow()
row = struct();
row.stage = NaN;
row.n_requested = NaN;
row.n_designs = NaN;
row.n_passing = NaN;
row.pass_rate = NaN;
row.n_errors = NaN;
row.n_yield_fail = NaN;
row.n_column_buckling_fail = NaN;
row.n_divergence_fail = NaN;
row.n_panel_buckling_fail = NaN;
row.best_weight_including_ribs_excluding_misc_lb = NaN;
row.best_index = NaN;
row.best_frontSpar_xc = NaN;
row.best_rearSpar_xc = NaN;
row.best_stringer_section_key = "";
row.best_spar_cap_sec_no = NaN;
row.best_main_skin_t = NaN;
row.best_noseSkin_t = NaN;
row.best_ribSpacing_in = NaN;
row.elite_count = NaN;
row.rngSeed = NaN;
end

function printRangeSummary(ranges, heading)
fprintf('%s\n', heading);
fprintf('  frontSpar_xc      = %s\n', mat2str(ranges.frontSpar_xc));
fprintf('  rearSpar_xc       = %s\n', mat2str(ranges.rearSpar_xc));
fprintf('  main_skin_t       = %s\n', mat2str(ranges.main_skin_t));
fprintf('  noseSkin_t        = %s\n', mat2str(ranges.noseSkin_t));
fprintf('  t1_frontSpar_in   = %s\n', mat2str(ranges.t1_frontSpar_in));
fprintf('  t2_rearSpar_in    = %s\n', mat2str(ranges.t2_rearSpar_in));
fprintf('  numTop            = %s\n', mat2str(ranges.numTop));
fprintf('  numBottom         = %s\n', mat2str(ranges.numBottom));
fprintf('  numNoseTop        = %s\n', mat2str(ranges.numNoseTop));
fprintf('  numNoseBottom     = %s\n', mat2str(ranges.numNoseBottom));
fprintf('  ribSpacing_in     = %s\n', mat2str(ranges.ribSpacing_in));
fprintf('  stringer sections = %s\n', sectionRangeSummary(ranges.stringer_section_key));
fprintf('  spar cap sections = %s\n', sectionRangeSummary(ranges.spar_cap_sec_no));
fprintf('  min separation    = %.3f\n', ranges.minSparSeparation_xc);
end

function text = sectionRangeSummary(values)
if isstring(values) || ischar(values) || iscellstr(values)
    values = string(values);
    values = values(:);
    nValues = numel(values);
    if nValues == 0
        text = "0 candidates";
    elseif nValues <= 10
        text = sprintf('%d candidates [%s]', nValues, strjoin(values', ", "));
    else
        text = sprintf('%d candidates, first %s, last %s', ...
            nValues, values(1), values(end));
    end
else
    values = values(:)';
    nValues = numel(values);
    if nValues == 0
        text = "0 candidates";
    elseif nValues <= 10
        text = sprintf('%d candidates %s', nValues, mat2str(values));
    else
        text = sprintf('%d candidates, first %g, last %g', ...
            nValues, values(1), values(end));
    end
end
end

function row = makeStageSummaryRowFromResult(stageNumber, nRequested, stageResult, nElite)
row = makeStageSummaryRow();
row.stage = stageNumber;
row.n_requested = nRequested;
row.n_designs = height(stageResult.allDesigns);
row.n_passing = height(stageResult.passingDesigns);
row.pass_rate = safeDivide(row.n_passing, row.n_designs);
row.n_errors = stageResult.failureCounts.n_errors(1);
row.n_yield_fail = stageResult.failureCounts.n_yield_fail(1);
row.n_column_buckling_fail = stageResult.failureCounts.n_column_buckling_fail(1);
row.n_divergence_fail = stageResult.failureCounts.n_divergence_fail(1);
row.n_panel_buckling_fail = stageResult.failureCounts.n_panel_buckling_fail(1);
row.elite_count = nElite;
row.rngSeed = stageResult.rngState.Seed;
if height(stageResult.bestDesign) > 0
    best = stageResult.bestDesign;
    row.best_weight_including_ribs_excluding_misc_lb = ...
        best.weight_including_ribs_excluding_misc_lb(1);
    row.best_index = best.index(1);
    row.best_frontSpar_xc = best.frontSpar_xc(1);
    row.best_rearSpar_xc = best.rearSpar_xc(1);
    row.best_stringer_section_key = string(best.stringer_section_key(1));
    row.best_spar_cap_sec_no = best.spar_cap_sec_no(1);
    row.best_main_skin_t = best.main_skin_t(1);
    row.best_noseSkin_t = best.noseSkin_t(1);
    row.best_ribSpacing_in = best.ribSpacing_in(1);
end
end

function row = makeRangeHistoryRow()
row = struct();
row.stage = NaN;
row.variable = "";
row.range_before = "";
row.elite_min = "";
row.elite_max = "";
row.range_after = "";
row.n_parent_values = NaN;
row.n_after_values = NaN;
end

function T = emptyRangeHistoryTable()
T = struct2table(repmat(makeRangeHistoryRow(), 0, 1));
end

function summary = makeAdaptiveSummary(adaptiveResults)
stageSummary = adaptiveResults.stageSummary;
nDesigns = sum(stageSummary.n_designs);
nPassing = sum(stageSummary.n_passing);
if height(adaptiveResults.bestDesign) > 0
    best = adaptiveResults.bestDesign;
    bestWeight = best.weight_including_ribs_excluding_misc_lb(1);
    bestYield = best.max_yield_ratio(1);
    bestColumn = best.max_column_buckling_ratio(1);
    bestDivergence = best.divergence_ratio(1);
    bestCombined = best.max_active_combined_panel_ratio(1);
    bestIndex = best.index(1);
else
    bestWeight = NaN; bestYield = NaN; bestColumn = NaN;
    bestDivergence = NaN; bestCombined = NaN; bestIndex = NaN;
end
lastStage = stageSummary(end, :);
summary = table( ...
    nDesigns, nPassing, safeDivide(nPassing, nDesigns), ...
    sum(stageSummary.n_errors), sum(stageSummary.n_yield_fail), ...
    sum(stageSummary.n_column_buckling_fail), ...
    sum(stageSummary.n_divergence_fail), ...
    sum(stageSummary.n_panel_buckling_fail), ...
    bestIndex, bestWeight, bestYield, bestColumn, bestDivergence, ...
    bestCombined, adaptiveResults.analysisOptions.loadFactor, ...
    string(adaptiveResults.analysisOptions.loadFactorPurpose), ...
    logical(adaptiveResults.analysisOptions.usePanelIdealisation), ...
    string(adaptiveResults.analysisOptions.panelIdealisationMethod), ...
    logical(adaptiveResults.analysisOptions.enablePanelBuckling), ...
    logical(adaptiveResults.analysisOptions.panelBuckling.enableCombinedRatios), ...
    logical(adaptiveResults.analysisOptions.panelBuckling.useBayAveragedRecovery), ...
    logical(adaptiveResults.analysisOptions.panelBuckling.enableCurvedPanelCorrection), ...
    nDesigns, string(adaptiveResults.optimisationOptions.seedMode), ...
    "adaptive_multi_stage", lastStage.rngSeed(1), "", ...
    adaptiveResults.optimisationOptions.optionCounts.n_valid_stringers, ...
    adaptiveResults.optimisationOptions.optionCounts.n_valid_spar_caps, ...
    string(adaptiveResults.optimisationOptions.optionFilters.stringers), ...
    string(adaptiveResults.optimisationOptions.optionFilters.sparCaps), ...
    'VariableNames', { ...
    'n_designs','n_passing','pass_rate','n_errors','n_yield_fail', ...
    'n_column_buckling_fail','n_divergence_fail','n_panel_buckling_fail', ...
    'best_index','best_weight_including_ribs_excluding_misc_lb', ...
    'best_yield_ratio','best_column_buckling_ratio', ...
    'best_divergence_ratio','best_active_combined_panel_ratio', ...
    'loadFactor','loadFactorPurpose','usePanelIdealisation', ...
    'panelIdealisationMethod','enablePanelBuckling', ...
    'enableCombinedPanelRatios','useBayAveragedRecovery', ...
    'enableCurvedPanelCorrection','requested_numDesigns','seedMode', ...
    'rngType','rngSeed','rngState','n_valid_stringer_options', ...
    'n_valid_spar_cap_options','stringer_filter','spar_cap_filter'});
end

function T = appendTables(T, newRows)
if isempty(T)
    T = newRows;
else
    T = [T; newRows];
end
end

function value = safeDivide(num, den)
if den == 0
    value = NaN;
else
    value = num / den;
end
end

function value = getFieldWithDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
