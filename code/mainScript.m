%% mainScript
% Wing-box structural optimisation driver.
%
% This driver loads the project inputs, runs the adaptive Monte Carlo
% wing-box search, reanalyses the selected design, and exports the report
% tables and figures. Detailed analysis settings are defined in
% defaultFinalAnalysisOptions.m; the user controls below define the
% optimisation run.


clear; clc; close all;

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
addpath(scriptFolder);

%% Final Monte Carlo user controls
% true runs the four-stage adaptive Monte Carlo search; false runs a
% single-stage random search.
ADAPTIVE_MONTE_CARLO = true;
% number of random wing designs evaluated in each adaptive stage.
N_PER_STAGE = [5000 5000 5000 5000];
NUM_STAGES = numel(N_PER_STAGE);
ELITE_FRACTION = 0.10;
MIN_ELITE_COUNT = 5;
RANGE_PADDING_STEPS = 1;
MIN_SECTION_CANDIDATES = 5;
% true keeps the best valid candidate found so far and prevents a heavier
% later run from overwriting it.
UPDATE_FINAL_ONLY_IF_IMPROVED = true;
% true uses MATLAB Parallel Computing Toolbox if available; false uses a
% normal serial loop. Set false for maximum compatibility on machines
% without the toolbox.
USE_PARALLEL = true;

% "shuffle" gives a new random search each run.
% "fixed" with SEED set to a number is useful for repeatable verification;
% serial mode is the most reproducible.
% seed mode is not guaranteed to work when using parallel computing;
% set USE_PARALLEL = false if you want to use a seed 
SEED_MODE = "shuffle";
SEED = [];

%% Final Monte Carlo design-variable ranges
FRONT_SPAR_XC_RANGE = [0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30];
% rear spar values start at 0.60c to satisfy the design constraint.
REAR_SPAR_XC_RANGE  = [0.60 0.62 0.66 0.70 0.74 0.78];

MAIN_SKIN_T_RANGE      = [0.034 0.038 0.040 0.042 0.045 0.050 0.056 0.063];
NOSE_SKIN_T_RANGE      = [0.040 0.050 0.056 0.060 0.063 0.067 0.071];
FRONT_SPAR_WEB_T_RANGE = [0.100 0.125 0.135 0.145 0.150 0.156 0.188];
REAR_SPAR_WEB_T_RANGE  = [0.050 0.063 0.070 0.080 0.090 0.100 0.125];

NUM_TOP_RANGE         = [10 12 14 16 18];
NUM_BOTTOM_RANGE      = [10 12 14 16 18];
NUM_NOSE_TOP_RANGE    = [2 3 4 5];
NUM_NOSE_BOTTOM_RANGE = [2 3 4 5];

RIB_SPACING_RANGE = [12 14 16 18 20 22 24];

MIN_SPAR_SEPARATION_XC = 0.20;

fprintf('Wing-box Monte Carlo optimisation\n');
fprintf('Project root: %s\n\n', projectRoot);

criticalLoadsFile = fullfile(projectRoot, 'data', 'critical_load_cases.csv');
spanwiseLoadsFile = fullfile(projectRoot, 'data', 'spanwise_loads.mat');

if ~isfile(criticalLoadsFile)
    error('mainScript:MissingLoadCases', ...
        'Missing aerodynamic load table: %s. Generate the saved load data first.', ...
        criticalLoadsFile);
end
if ~isfile(spanwiseLoadsFile)
    error('mainScript:MissingSpanwiseLoads', ...
        'Missing spanwise load file: %s. Final analysis requires spanwise loads.', ...
        spanwiseLoadsFile);
end

loadCaseTable = readtable(criticalLoadsFile, 'TextType', 'string');
spanwiseData = load(spanwiseLoadsFile);
if ~isfield(spanwiseData, 'spanwiseLoads')
    error('mainScript:InvalidSpanwiseFile', ...
        'Spanwise load file does not contain spanwiseLoads: %s', ...
        spanwiseLoadsFile);
end

loads = criticalLoads(loadCaseTable, spanwiseData.spanwiseLoads);
designOptions = loadDesignOptions();
analysisOptions = defaultFinalAnalysisOptions();
optimisationOptions = defaultFinalOptimisationOptions(designOptions);
panelIdealisationDisplay = "modified_megson_panel_idealisation";
optimisationOptions.numDesigns = N_PER_STAGE(1);
optimisationOptions.seedMode = SEED_MODE;
optimisationOptions.seed = SEED;
optimisationOptions.adaptive.enabled = ADAPTIVE_MONTE_CARLO;
optimisationOptions.adaptive.nPerStage = N_PER_STAGE;
optimisationOptions.adaptive.numStages = NUM_STAGES;
optimisationOptions.adaptive.eliteFraction = ELITE_FRACTION;
optimisationOptions.adaptive.minEliteCount = MIN_ELITE_COUNT;
optimisationOptions.adaptive.rangePaddingSteps = RANGE_PADDING_STEPS;
optimisationOptions.adaptive.minSectionCandidates = MIN_SECTION_CANDIDATES;
optimisationOptions.updateFinalOnlyIfImproved = UPDATE_FINAL_ONLY_IF_IMPROVED;
optimisationOptions.useParallel = USE_PARALLEL;
optimisationOptions.ranges.frontSpar_xc = FRONT_SPAR_XC_RANGE;
optimisationOptions.ranges.rearSpar_xc = REAR_SPAR_XC_RANGE;
optimisationOptions.ranges.main_skin_t = MAIN_SKIN_T_RANGE;
optimisationOptions.ranges.noseSkin_t = NOSE_SKIN_T_RANGE;
optimisationOptions.ranges.t1_frontSpar_in = FRONT_SPAR_WEB_T_RANGE;
optimisationOptions.ranges.t2_rearSpar_in = REAR_SPAR_WEB_T_RANGE;
optimisationOptions.ranges.numTop = NUM_TOP_RANGE;
optimisationOptions.ranges.numBottom = NUM_BOTTOM_RANGE;
optimisationOptions.ranges.numNoseTop = NUM_NOSE_TOP_RANGE;
optimisationOptions.ranges.numNoseBottom = NUM_NOSE_BOTTOM_RANGE;
optimisationOptions.ranges.ribSpacing_in = RIB_SPACING_RANGE;
optimisationOptions.constraints.minSparSeparation_xc = MIN_SPAR_SEPARATION_XC;
optimisationOptions.ranges.minSparSeparation_xc = MIN_SPAR_SEPARATION_XC;

fprintf('Loaded %d aerodynamic load cases.\n', loads.nCases);
fprintf('Loaded %d spanwise source stations.\n', numel(loads.spanwise.y_ft));
fprintf('Valid final stringer options: %d (%s)\n', ...
    optimisationOptions.optionCounts.n_valid_stringers, ...
    optimisationOptions.optionFilters.stringers);
fprintf('Valid final spar cap options: %d (%s)\n\n', ...
    optimisationOptions.optionCounts.n_valid_spar_caps, ...
    optimisationOptions.optionFilters.sparCaps);

fprintf('Final analysis options:\n');
fprintf('- loadFactor: %.3f (%s)\n', ...
    analysisOptions.loadFactor, analysisOptions.loadFactorPurpose);
fprintf('- panel idealisation: %s (%s)\n', ...
    string(analysisOptions.usePanelIdealisation), ...
    panelIdealisationDisplay);
fprintf('- panel buckling enabled: %s\n', string(analysisOptions.enablePanelBuckling));
fprintf('- combined ratios: %s\n', ...
    string(analysisOptions.panelBuckling.enableCombinedRatios));
fprintf('- bay-averaged panel buckling recovery: %s\n', ...
    string(analysisOptions.panelBuckling.useBayAveragedRecovery));
fprintf('- cylindrical curvature correction: %s\n\n', ...
    string(analysisOptions.panelBuckling.enableCurvedPanelCorrection));

fprintf('Optimisation options:\n');
fprintf('- adaptive mode: %s\n', string(ADAPTIVE_MONTE_CARLO));
fprintf('- NUM_STAGES: %d\n', NUM_STAGES);
fprintf('- N_PER_STAGE: %s\n', mat2str(N_PER_STAGE));
fprintf('- single-stage fallback N: %d\n', optimisationOptions.numDesigns);
fprintf('- numDesigns: %d\n', optimisationOptions.numDesigns);
fprintf('- seedMode: %s\n', optimisationOptions.seedMode);
fprintf('- elite fraction: %.3f\n', ELITE_FRACTION);
fprintf('- min elite count: %d\n', MIN_ELITE_COUNT);
fprintf('- range padding steps: %d\n', RANGE_PADDING_STEPS);
fprintf('- min section candidates: %d\n', MIN_SECTION_CANDIDATES);
fprintf('- update final only if improved: %s\n', string(UPDATE_FINAL_ONLY_IF_IMPROVED));
fprintf('- use parallel analysis: %s\n', string(USE_PARALLEL));
fprintf('- objective: %s\n\n', optimisationOptions.objective);
fprintf('Active user design-variable ranges:\n');
fprintf('- frontSpar_xc: %s\n', mat2str(optimisationOptions.ranges.frontSpar_xc));
fprintf('- rearSpar_xc: %s\n', mat2str(optimisationOptions.ranges.rearSpar_xc));
fprintf('- main skin t: %s in\n', mat2str(optimisationOptions.ranges.main_skin_t));
fprintf('- nose skin t: %s in\n', mat2str(optimisationOptions.ranges.noseSkin_t));
fprintf('- front spar web t: %s in\n', mat2str(optimisationOptions.ranges.t1_frontSpar_in));
fprintf('- rear spar web t: %s in\n', mat2str(optimisationOptions.ranges.t2_rearSpar_in));
fprintf('- numTop: %s\n', mat2str(optimisationOptions.ranges.numTop));
fprintf('- numBottom: %s\n', mat2str(optimisationOptions.ranges.numBottom));
fprintf('- numNoseTop: %s\n', mat2str(optimisationOptions.ranges.numNoseTop));
fprintf('- numNoseBottom: %s\n', mat2str(optimisationOptions.ranges.numNoseBottom));
fprintf('- ribSpacing_in: %s in\n', mat2str(optimisationOptions.ranges.ribSpacing_in));
fprintf('- min spar separation x/c: %.3f\n\n', ...
    optimisationOptions.ranges.minSparSeparation_xc);

if ADAPTIVE_MONTE_CARLO
    results = runAdaptiveFinalMonteCarlo(loads, designOptions, analysisOptions, optimisationOptions);

    totalDesigns = height(results.allDesigns);
    totalPassing = nnz(results.allDesigns.overall_pass == true & ...
        results.allDesigns.error_flag == false);
    fprintf('\nAdaptive final Monte Carlo summary\n');
    fprintf('Stages completed: %d / %d\n', results.completedStages, NUM_STAGES);
    fprintf('Total designs evaluated: %d\n', totalDesigns);
    fprintf('Total passing designs: %d (%.2f%%)\n', totalPassing, ...
        100 * totalPassing / max(totalDesigns, 1));
    fprintf('Seed mode: %s\n', SEED_MODE);
    fprintf('\nPer-stage summary:\n');
    for iStage = 1:height(results.stageSummary)
        row = results.stageSummary(iStage, :);
        fprintf('- Stage %d: N = %d, passing = %d (%.2f%%), best weight = %.4f lb, elite count = %d, seed = %d\n', ...
            row.stage(1), row.n_requested(1), row.n_passing(1), ...
            100 * row.pass_rate(1), ...
            row.best_weight_including_ribs_excluding_misc_lb(1), ...
            row.elite_count(1), row.rngSeed(1));
    end

    if isfield(results, 'stage4BestDesign') && height(results.stage4BestDesign) > 0
        stage4Best = results.stage4BestDesign;
        fprintf('\nFinal-stage lightest passing design\n');
        fprintf('front/rear spar x/c: %.3f / %.3f\n', ...
            stage4Best.frontSpar_xc(1), stage4Best.rearSpar_xc(1));
        fprintf('Stringer: %s (sec_no %.0f)\n', ...
            stage4Best.stringer_section_key(1), stage4Best.stringer_sec_no(1));
        fprintf('Spar cap sec_no: %.0f\n', stage4Best.spar_cap_sec_no(1));
        fprintf('main skin / nose skin: %.4f / %.4f in\n', ...
            stage4Best.main_skin_t(1), stage4Best.noseSkin_t(1));
        fprintf('front/rear spar web: %.4f / %.4f in\n', ...
            stage4Best.t1_frontSpar_in(1), stage4Best.t2_rearSpar_in(1));
        fprintf('top/bottom/noseTop/noseBottom stringers: %.0f / %.0f / %.0f / %.0f\n', ...
            stage4Best.numTop(1), stage4Best.numBottom(1), ...
            stage4Best.numNoseTop(1), stage4Best.numNoseBottom(1));
        fprintf('rib spacing: %.1f in\n', stage4Best.ribSpacing_in(1));
        fprintf('weight including ribs/excluding misc: %.4f lb\n', ...
            stage4Best.weight_including_ribs_excluding_misc_lb(1));
        fprintf('yield / column / divergence: %.6f / %.6f / %.6f\n', ...
            stage4Best.max_yield_ratio(1), stage4Best.max_column_buckling_ratio(1), ...
            stage4Best.divergence_ratio(1));
        fprintf('Rc / Rb / Rs / active combined: %.6f / %.6f / %.6f / %.6f\n', ...
            stage4Best.max_Rc(1), stage4Best.max_Rb(1), stage4Best.max_Rs(1), ...
            stage4Best.max_active_combined_panel_ratio(1));
        fprintf('Controlling combined: bay %.0f, %s, %s, %s\n', ...
            stage4Best.controlling_combined_bay(1), ...
            stage4Best.controlling_combined_case(1), ...
            stage4Best.controlling_combined_panel(1), ...
            stage4Best.controlling_combined_mode(1));
    else
        fprintf('\nNo final-stage passing design found.\n');
    end

    if height(results.bestDesign) > 0
        best = results.bestDesign;
        fprintf('\nGlobal best passing design across all stages, used for active final candidate update\n');
        fprintf('stage/index: %.0f / %.0f\n', best.stage(1), best.index(1));
        fprintf('front/rear spar x/c: %.3f / %.3f\n', ...
            best.frontSpar_xc(1), best.rearSpar_xc(1));
        fprintf('Stringer: %s (sec_no %.0f), spar cap sec_no %.0f\n', ...
            best.stringer_section_key(1), best.stringer_sec_no(1), ...
            best.spar_cap_sec_no(1));
        fprintf('weight including ribs/excluding misc: %.4f lb\n', ...
            best.weight_including_ribs_excluding_misc_lb(1));
        fprintf('yield / column / divergence / active combined: %.6f / %.6f / %.6f / %.6f\n', ...
            best.max_yield_ratio(1), best.max_column_buckling_ratio(1), ...
            best.divergence_ratio(1), best.max_active_combined_panel_ratio(1));
        if isfield(results, 'stage4BestDesign') && height(results.stage4BestDesign) > 0
            stage4BestWeight = results.stage4BestDesign.weight_including_ribs_excluding_misc_lb(1);
            globalBestWeight = best.weight_including_ribs_excluding_misc_lb(1);
            if abs(globalBestWeight - stage4BestWeight) > 1e-6
                fprintf('Global best across all stages differs from final-stage best: %.4f lb vs %.4f lb.\n', ...
                    globalBestWeight, stage4BestWeight);
            else
                fprintf('Global best across all stages matches the final-stage best by weight.\n');
            end
        end
    else
        fprintf('\nNo passing design found in any adaptive stage. Active final candidate artifacts were not overwritten.\n');
    end

    if isfield(results, 'finalCandidateUpdate') && ...
            isfield(results.finalCandidateUpdate, 'adaptiveOutputsFolder')
        fprintf('\nAdaptive Monte Carlo outputs for this run were written to:\n');
        fprintf('- %s\n', results.finalCandidateUpdate.adaptiveOutputsFolder);
    else
        fprintf('\nAdaptive Monte Carlo CSV outputs are in outputs/tables.\n');
    end
else
    results = runFinalMonteCarlo(loads, designOptions, analysisOptions, optimisationOptions);

    nDesigns = height(results.allDesigns);
    nPassing = height(results.passingDesigns);
    if nDesigns > 0
        passPercent = 100.0 * nPassing / nDesigns;
    else
        passPercent = NaN;
    end

    fprintf('\nFinal Monte Carlo summary\n');
    fprintf('Designs evaluated: %d\n', nDesigns);
    fprintf('Passing designs: %d (%.2f%%)\n', nPassing, passPercent);
    fprintf('Errors: %d\n', results.failureCounts.n_errors(1));
    fprintf('Yield failures: %d\n', results.failureCounts.n_yield_fail(1));
    fprintf('Column buckling failures: %d\n', ...
        results.failureCounts.n_column_buckling_fail(1));
    fprintf('Divergence failures: %d\n', results.failureCounts.n_divergence_fail(1));
    fprintf('Panel buckling failures: %d\n', ...
        results.failureCounts.n_panel_buckling_fail(1));
    fprintf('RNG: mode = %s, type = %s, seed = %d\n', ...
        results.rngState.seedMode, results.rngState.Type, results.rngState.Seed);

    if nPassing > 0
        best = results.bestDesign;
        bestDesign = results.bestDesignStruct;
        fprintf('\nLightest passing design found\n');
        fprintf('front/rear spar x/c: %.3f / %.3f\n', ...
            best.frontSpar_xc(1), best.rearSpar_xc(1));
        fprintf('Stringer: %s (sec_no %.0f)\n', ...
            best.stringer_section_key(1), best.stringer_sec_no(1));
        fprintf('Spar cap sec_no: %.0f\n', best.spar_cap_sec_no(1));
        fprintf('main skin / nose skin: %.4f / %.4f in\n', ...
            best.main_skin_t(1), best.noseSkin_t(1));
        fprintf('front/rear spar web: %.4f / %.4f in\n', ...
            best.t1_frontSpar_in(1), best.t2_rearSpar_in(1));
        fprintf('top/bottom/noseTop/noseBottom stringers: %.0f / %.0f / %.0f / %.0f\n', ...
            best.numTop(1), best.numBottom(1), ...
            best.numNoseTop(1), best.numNoseBottom(1));
        fprintf('rib spacing: %.1f in\n', best.ribSpacing_in(1));
        fprintf('weight including ribs/excluding misc: %.4f lb\n', ...
            best.weight_including_ribs_excluding_misc_lb(1));
        fprintf('yield / column / divergence: %.6f / %.6f / %.6f\n', ...
            best.max_yield_ratio(1), best.max_column_buckling_ratio(1), ...
            best.divergence_ratio(1));
        fprintf('Rc / Rb / Rs / active combined: %.6f / %.6f / %.6f / %.6f\n', ...
            best.max_Rc(1), best.max_Rb(1), best.max_Rs(1), ...
            best.max_active_combined_panel_ratio(1));
        fprintf('Controlling combined: bay %.0f, %s, %s, %s\n', ...
            best.controlling_combined_bay(1), ...
            best.controlling_combined_case(1), ...
            best.controlling_combined_panel(1), ...
            best.controlling_combined_mode(1));

        fprintf('\nRe-running detailed spanwise analysis for best passing design...\n');
        bestDetailed = analyzeWingDesignSpanwise(bestDesign, loads, designOptions, analysisOptions);
        fprintf('Detailed rerun overall pass: %s\n', string(bestDetailed.summary.overall_pass));
        fprintf('Detailed rerun weight: %.4f lb\n', ...
            bestDetailed.summary.weight_including_ribs_excluding_misc_lb);

        fprintf('\nRegenerating active final candidate tables, notes, and figures...\n');
        exportFinalCandidateArtifacts(bestDesign, bestDetailed, analysisOptions, results);
    else
        fprintf('\nNo passing design found in this validation run. Increase numDesigns or adjust ranges after review.\n');
        fprintf('Previous final candidate artifacts were not overwritten.\n');
    end

    fprintf('\nFinal Monte Carlo CSV outputs are in outputs/tables:\n');
    fprintf('- final_monte_carlo_results.csv\n');
    fprintf('- final_monte_carlo_summary.csv\n');
    fprintf('- final_monte_carlo_best_design.csv\n');
end
