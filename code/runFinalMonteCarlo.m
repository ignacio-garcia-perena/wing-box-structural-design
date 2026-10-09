function results = runFinalMonteCarlo(loads, designOptions, analysisOptions, optimisationOptions)
% runFinalMonteCarlo Run final-design random wing optimisation.
%
% The final workflow evaluates each random design with
% analyzeWingDesignSpanwise and the final structural checks.

if nargin < 1 || isempty(loads) || ~isstruct(loads)
    error('runFinalMonteCarlo:MissingLoads', ...
        'loads struct from criticalLoads is required.');
end
if nargin < 2 || isempty(designOptions) || ~isstruct(designOptions)
    error('runFinalMonteCarlo:MissingDesignOptions', ...
        'designOptions struct from loadDesignOptions is required.');
end
if nargin < 3 || isempty(analysisOptions)
    analysisOptions = defaultFinalAnalysisOptions();
end
if nargin < 4 || isempty(optimisationOptions)
    optimisationOptions = defaultFinalOptimisationOptions(designOptions);
end

validateInputs(analysisOptions, optimisationOptions);
rngState = initializeRng(optimisationOptions);

N = optimisationOptions.numDesigns;
rows = repmat(makeEmptyResultRow(), N, 1);
designs = cell(N, 1);

useParallel = shouldUseParallel(optimisationOptions);
if useParallel
    [rows, designs] = generateRandomDesignRows(N, designOptions, optimisationOptions);
    parfor iDesign = 1:N
        rows(iDesign) = evaluateGeneratedDesign(rows(iDesign), designs{iDesign}, ...
            loads, designOptions, analysisOptions);
    end
else
    for iDesign = 1:N
        row = makeEmptyResultRow();
        row.index = iDesign;
        try
            design = randomFinalDesign(designOptions, optimisationOptions);
            designs{iDesign} = design;
            row = fillDesignFields(row, design);

            result = analyzeWingDesignSpanwise(design, loads, designOptions, analysisOptions);
            row = fillResultFields(row, result);
            clear result
        catch ME
            row.error_flag = true;
            row.error_message = string(ME.message);
        end
        rows(iDesign) = row;
    end
end

allDesigns = struct2table(rows);
passingDesigns = allDesigns(allDesigns.overall_pass == true & ...
    allDesigns.error_flag == false, :);
if height(passingDesigns) > 0
    passingDesigns = sortrows(passingDesigns, ...
        'weight_including_ribs_excluding_misc_lb', 'ascend');
    bestDesign = passingDesigns(1, :);
    bestIndex = bestDesign.index(1);
else
    bestDesign = passingDesigns;
    bestIndex = NaN;
end

failureCounts = makeFailureCounts(allDesigns);
summary = summarizeFinalMonteCarloResults(allDesigns, passingDesigns, ...
    failureCounts, analysisOptions, optimisationOptions, rngState);

results = struct();
results.allDesigns = allDesigns;
results.passingDesigns = passingDesigns;
results.bestDesign = bestDesign;
results.bestDesignStruct = [];
if isfinite(bestIndex)
    results.bestDesignStruct = designs{bestIndex};
end
results.failureCounts = failureCounts;
results.summary = summary;
results.analysisOptions = analysisOptions;
results.optimisationOptions = optimisationOptions;
results.rngState = rngState;
results.useParallel = useParallel;
results.designs = designs;

exportFinalMonteCarloResults(results);
end

function [rows, designs] = generateRandomDesignRows(N, designOptions, optimisationOptions)
rows = repmat(makeEmptyResultRow(), N, 1);
designs = cell(N, 1);
for iDesign = 1:N
    row = makeEmptyResultRow();
    row.index = iDesign;
    try
        design = randomFinalDesign(designOptions, optimisationOptions);
        designs{iDesign} = design;
        row = fillDesignFields(row, design);
    catch ME
        row.error_flag = true;
        row.error_message = string(ME.message);
    end
    rows(iDesign) = row;
end
end

function row = evaluateGeneratedDesign(row, design, loads, designOptions, analysisOptions)
if row.error_flag || isempty(design)
    return
end
try
    result = analyzeWingDesignSpanwise(design, loads, designOptions, analysisOptions);
    row = fillResultFields(row, result);
catch ME
    row.error_flag = true;
    row.error_message = string(ME.message);
end
end

function useParallel = shouldUseParallel(optimisationOptions)
useParallel = false;
if ~isfield(optimisationOptions, 'useParallel') || ~optimisationOptions.useParallel
    return
end
if ~parallelToolboxAvailable()
    fprintf('Parallel execution requested, but Parallel Computing Toolbox is not available. Using serial for-loop.\n');
    return
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
    useParallel = true;
    fprintf('Parallel execution enabled for Monte Carlo design analysis.\n');
catch ME
    fprintf('Parallel execution requested but could not start a pool: %s\n', ME.message);
    fprintf('Using serial for-loop.\n');
end
end

function tf = parallelToolboxAvailable()
tf = false;
try
    hasLicense = license('test', 'Distrib_Computing_Toolbox');
    hasParpool = exist('parpool', 'file') == 2 || exist('parpool', 'file') == 6;
    hasGcp = exist('gcp', 'file') == 2 || exist('gcp', 'file') == 6;
    tf = hasLicense && hasParpool && hasGcp;
catch
    tf = false;
end
end

function validateInputs(analysisOptions, optimisationOptions)
if ~isfield(analysisOptions, 'loadFactor') || analysisOptions.loadFactor < 1.5
    error('runFinalMonteCarlo:InvalidFinalLoadFactor', ...
        'Final Monte Carlo requires analysisOptions.loadFactor >= 1.5.');
end
if ~isfield(analysisOptions, 'usePanelIdealisation') || ...
        ~analysisOptions.usePanelIdealisation
    error('runFinalMonteCarlo:PanelIdealisationDisabled', ...
        'Final Monte Carlo requires panel idealisation enabled.');
end
if ~isfield(analysisOptions, 'enablePanelBuckling') || ...
        ~analysisOptions.enablePanelBuckling
    error('runFinalMonteCarlo:PanelBucklingDisabled', ...
        'Final Monte Carlo requires panel buckling enabled.');
end
requiredPanelOptions = ["enableCombinedRatios", ...
    "useBayAveragedRecovery", "enableCurvedPanelCorrection"];
missing = setdiff(requiredPanelOptions, string(fieldnames(analysisOptions.panelBuckling)));
if ~isempty(missing)
    error('runFinalMonteCarlo:MissingPanelOptions', ...
        'analysisOptions.panelBuckling is missing fields: %s', strjoin(missing, ', '));
end
if ~analysisOptions.panelBuckling.enableCombinedRatios || ...
        ~analysisOptions.panelBuckling.useBayAveragedRecovery || ...
        ~analysisOptions.panelBuckling.enableCurvedPanelCorrection
    error('runFinalMonteCarlo:IncompleteFinalPanelOptions', ...
        'Final Monte Carlo requires combined, bay-averaged, and curved panel buckling options enabled.');
end
if ~isfield(optimisationOptions, 'numDesigns') || ...
        ~(isnumeric(optimisationOptions.numDesigns) && isscalar(optimisationOptions.numDesigns) && ...
        isfinite(optimisationOptions.numDesigns) && optimisationOptions.numDesigns > 0 && ...
        mod(optimisationOptions.numDesigns, 1) == 0)
    error('runFinalMonteCarlo:InvalidN', ...
        'optimisationOptions.numDesigns must be a positive integer.');
end
end

function rngState = initializeRng(optimisationOptions)
seedMode = string(optimisationOptions.seedMode);
switch seedMode
    case "shuffle"
        rng('shuffle');
    case "fixed"
        if ~isfield(optimisationOptions, 'seed') || ...
                ~(isnumeric(optimisationOptions.seed) && isscalar(optimisationOptions.seed) && ...
                isfinite(optimisationOptions.seed))
            error('runFinalMonteCarlo:InvalidSeed', ...
                'optimisationOptions.seed must be finite when seedMode = fixed.');
        end
        rng(optimisationOptions.seed);
    otherwise
        error('runFinalMonteCarlo:InvalidSeedMode', ...
            'seedMode must be "shuffle" or "fixed".');
end
rngState = rng;
rngState.seedMode = seedMode;
rngState.stateString = mat2str(rngState.State(:)');
end

function row = makeEmptyResultRow()
row = struct();
row.index = NaN;
row.frontSpar_xc = NaN;
row.rearSpar_xc = NaN;
row.stringer_section_key = "";
row.stringer_sec_no = NaN;
row.spar_cap_sec_no = NaN;
row.main_skin_t = NaN;
row.noseSkin_t = NaN;
row.t1_frontSpar_in = NaN;
row.t2_rearSpar_in = NaN;
row.numTop = NaN;
row.numBottom = NaN;
row.numNoseTop = NaN;
row.numNoseBottom = NaN;
row.ribSpacing_in = NaN;
row.K_column = NaN;
row.rib_count_fullspan = NaN;
row.rib_weight_fullspan_lb = NaN;
row.weight_excluding_ribs_misc_lb = NaN;
row.weight_including_ribs_excluding_misc_lb = NaN;
row.max_yield_ratio = NaN;
row.max_column_buckling_ratio = NaN;
row.divergence_ratio = NaN;
row.max_Rc = NaN;
row.max_Rb = NaN;
row.max_Rs = NaN;
row.max_R_SC = NaN;
row.max_R_BS = NaN;
row.max_active_combined_panel_ratio = NaN;
row.controlling_combined_bay = NaN;
row.controlling_combined_bay_start_in = NaN;
row.controlling_combined_bay_end_in = NaN;
row.controlling_combined_case = "";
row.controlling_combined_panel = "";
row.controlling_combined_mode = "";
row.yield_pass = false;
row.column_buckling_pass = false;
row.divergence_pass = false;
row.panel_buckling_pass = false;
row.overall_pass = false;
row.error_flag = false;
row.error_message = "";
end

function row = fillDesignFields(row, design)
row.frontSpar_xc = design.box.frontSpar_xc;
row.rearSpar_xc = design.box.rearSpar_xc;
row.stringer_section_key = string(design.stringers.section_key);
row.stringer_sec_no = design.stringers.sec_no;
row.spar_cap_sec_no = design.sparCaps.sec_no;
row.main_skin_t = design.thickness.t3_upperMainSkin_in;
row.noseSkin_t = design.thickness.t5_noseSkin_in;
row.t1_frontSpar_in = design.thickness.t1_frontSpar_in;
row.t2_rearSpar_in = design.thickness.t2_rearSpar_in;
row.numTop = design.stringers.numTop;
row.numBottom = design.stringers.numBottom;
row.numNoseTop = design.stringers.numNoseTop;
row.numNoseBottom = design.stringers.numNoseBottom;
row.ribSpacing_in = design.ribs.ribSpacing_in;
row.K_column = design.columnBuckling.K_column;
end

function row = fillResultFields(row, result)
summary = result.summary;
row.rib_count_fullspan = result.weight.components.ribs.count_fullspan;
row.rib_weight_fullspan_lb = result.weight.components.ribs.fullspan_weight_lb;
row.weight_excluding_ribs_misc_lb = summary.weight_excluding_ribs_misc_lb;
row.weight_including_ribs_excluding_misc_lb = ...
    summary.weight_including_ribs_excluding_misc_lb;
row.max_yield_ratio = summary.max_yield_ratio;
row.max_column_buckling_ratio = summary.max_column_buckling_ratio;
row.divergence_ratio = summary.divergence_ratio;
row.max_Rc = summary.panelBuckling.max_Rc;
row.max_Rb = summary.panelBuckling.max_Rb;
row.max_Rs = summary.panelBuckling.max_Rs;
row.max_R_SC = summary.panelBuckling.max_R_SC;
row.max_R_BS = summary.panelBuckling.max_R_BS;
row.max_active_combined_panel_ratio = summary.panelBuckling.max_combined_active;
if isfield(summary.panelBuckling, 'maxCombinedActive')
    item = summary.panelBuckling.maxCombinedActive;
    row.controlling_combined_bay = getFieldOrNaN(item, 'bay_index');
    row.controlling_combined_bay_start_in = getFieldOrNaN(item, 'bay_start_y_in');
    row.controlling_combined_bay_end_in = getFieldOrNaN(item, 'bay_end_y_in');
    row.controlling_combined_case = string(getFieldOrBlank(item, 'case_name'));
    row.controlling_combined_panel = string(getFieldOrBlank(item, 'panel_name'));
    row.controlling_combined_mode = string(getFieldOrBlank(item, 'combined_mode'));
end
row.yield_pass = summary.yield_pass;
row.column_buckling_pass = summary.column_buckling_pass;
row.divergence_pass = summary.divergence_pass;
row.panel_buckling_pass = summary.panel_buckling_pass;
row.overall_pass = summary.overall_pass;
end

function value = getFieldOrNaN(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = NaN;
end
end

function value = getFieldOrBlank(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = "";
end
end

function failureCounts = makeFailureCounts(allDesigns)
validRows = allDesigns.error_flag == false;
failureCounts = table( ...
    height(allDesigns), ...
    nnz(allDesigns.error_flag), ...
    nnz(allDesigns.yield_pass == false & validRows), ...
    nnz(allDesigns.column_buckling_pass == false & validRows), ...
    nnz(allDesigns.divergence_pass == false & validRows), ...
    nnz(allDesigns.panel_buckling_pass == false & validRows), ...
    nnz(allDesigns.overall_pass == true & validRows), ...
    'VariableNames', { ...
    'n_designs', ...
    'n_errors', ...
    'n_yield_fail', ...
    'n_column_buckling_fail', ...
    'n_divergence_fail', ...
    'n_panel_buckling_fail', ...
    'n_overall_pass'});
end
