function optimisationOptions = defaultFinalOptimisationOptions(designOptions)
% defaultFinalOptimisationOptions Return final Monte Carlo design ranges.
%
% Inputs:
%   designOptions - Optional struct from loadDesignOptions. When supplied,
%                   valid A3.10 Z-stringer and A3.13 T-spar-cap candidates
%                   are discovered from the active option tables.

if nargin < 1
    designOptions = [];
end

optimisationOptions = struct();
optimisationOptions.numDesigns = 10;
optimisationOptions.seedMode = "shuffle";
optimisationOptions.seed = NaN;
optimisationOptions.useParallel = false;
optimisationOptions.objective = "minimize_weight_including_ribs_excluding_misc";

optimisationOptions.ranges.frontSpar_xc = [0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22];
optimisationOptions.ranges.rearSpar_xc = [0.58, 0.60, 0.62, 0.64, 0.66, 0.68];
optimisationOptions.ranges.minSparSeparation_xc = 0.20;
optimisationOptions.ranges.main_skin_t = [0.038, 0.040, 0.042, 0.045, 0.050];
optimisationOptions.ranges.noseSkin_t = [0.056, 0.060, 0.063, 0.067];
optimisationOptions.ranges.t1_frontSpar_in = [0.125, 0.135, 0.145, 0.150, 0.156];
optimisationOptions.ranges.t2_rearSpar_in = [0.063, 0.070, 0.080, 0.090, 0.100];
optimisationOptions.ranges.numTop = [12, 14, 16];
optimisationOptions.ranges.numBottom = [12, 14, 16];
optimisationOptions.ranges.numNoseTop = [2, 3, 4];
optimisationOptions.ranges.numNoseBottom = [2, 3, 4];
optimisationOptions.ranges.ribSpacing_in = [16, 18, 20];
optimisationOptions.ranges.K_column = 0.7;

if ~isempty(designOptions)
    [stringerKeys, stringerFilter] = validFinalStringerKeys(designOptions);
    [sparCapSecNos, sparCapFilter] = validFinalSparCapSecNos(designOptions);
    optimisationOptions.ranges.stringer_section_key = stringerKeys;
    optimisationOptions.ranges.spar_cap_sec_no = sparCapSecNos;
    optimisationOptions.optionFilters.stringers = stringerFilter;
    optimisationOptions.optionFilters.sparCaps = sparCapFilter;
    optimisationOptions.optionCounts.n_valid_stringers = numel(stringerKeys);
    optimisationOptions.optionCounts.n_valid_spar_caps = numel(sparCapSecNos);
else
    optimisationOptions.ranges.stringer_section_key = strings(0, 1);
    optimisationOptions.ranges.spar_cap_sec_no = [];
    optimisationOptions.optionFilters.stringers = ...
        "section_key starts with A3_10_Z_ and section_type contains Z-section";
    optimisationOptions.optionFilters.sparCaps = ...
        "section_type contains T-section";
    optimisationOptions.optionCounts.n_valid_stringers = 0;
    optimisationOptions.optionCounts.n_valid_spar_caps = 0;
end

optimisationOptions.notes = [
    "Validation default is numDesigns = 10; increase for final optimisation runs."
    "seedMode = shuffle samples a fresh random design set on each run."
    "Final pass/fail is evaluated by analyzeWingDesignSpanwise with defaultFinalAnalysisOptions."
    ];
end

function [keys, filterDescription] = validFinalStringerKeys(designOptions)
if ~isfield(designOptions, 'stringers')
    error('defaultFinalOptimisationOptions:MissingStringers', ...
        'designOptions.stringers is required.');
end
stringers = designOptions.stringers;
required = ["section_key", "section_type"];
missing = setdiff(required, string(stringers.Properties.VariableNames));
if ~isempty(missing)
    error('defaultFinalOptimisationOptions:MissingStringerColumns', ...
        'Stringer table is missing columns: %s', strjoin(missing, ', '));
end
mask = startsWith(string(stringers.section_key), "A3_10_Z_") & ...
    contains(string(stringers.section_type), "Z-section", 'IgnoreCase', true);
keys = string(stringers.section_key(mask));
keys = keys(:);
if isempty(keys)
    error('defaultFinalOptimisationOptions:NoStringerCandidates', ...
        'No valid A3.10 Z-section stringer candidates were found.');
end
filterDescription = ...
    "section_key starts with A3_10_Z_ and section_type contains Z-section";
end

function [secNos, filterDescription] = validFinalSparCapSecNos(designOptions)
if ~isfield(designOptions, 'sparCaps')
    error('defaultFinalOptimisationOptions:MissingSparCaps', ...
        'designOptions.sparCaps is required.');
end
sparCaps = designOptions.sparCaps;
required = ["sec_no", "section_type"];
missing = setdiff(required, string(sparCaps.Properties.VariableNames));
if ~isempty(missing)
    error('defaultFinalOptimisationOptions:MissingSparCapColumns', ...
        'Spar cap table is missing columns: %s', strjoin(missing, ', '));
end
mask = contains(string(sparCaps.section_type), "T-section", 'IgnoreCase', true);
secNos = sparCaps.sec_no(mask);
secNos = secNos(:)';
if isempty(secNos)
    error('defaultFinalOptimisationOptions:NoSparCapCandidates', ...
        'No valid A3.13 T-section spar cap candidates were found.');
end
filterDescription = "section_type contains T-section";
end
