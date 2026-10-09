function spanLoads = interpolateSpanwiseLoads(loads, wing, analysisOptions)
% interpolateSpanwiseLoads Interpolate saved spanwise loads at rib stations.
%
% Inputs:
%   loads           - Struct from criticalLoads with loads.spanwise data.
%   wing            - Struct from setGeometry with rib spacing and semispan.
%   analysisOptions - Struct from defaultAnalysisOptions.
%
% Outputs:
%   spanLoads - Struct containing unscaled and load-factor-scaled spanwise
%               loads at rib stations in structural-analysis units.
%

validateInputs(loads, wing, analysisOptions);

loadFactor = analysisOptions.loadFactor;
caseNames = string(loads.caseNames(:));
nCases = loads.nCases;

ySource_ft = loads.spanwise.y_ft(:);
ySource_in = loads.spanwise.y_in(:);
stations_in = makeRibStations(wing.ribs.spacing_in, wing.planform.semispan_in);
stations_ft = stations_in / 12.0;

unscaled = struct();
unscaled.Vx_lb = zeros(nCases, numel(stations_in));
unscaled.Vz_lb = zeros(nCases, numel(stations_in));
unscaled.Mx_lbin = zeros(nCases, numel(stations_in));
unscaled.Mz_lbin = zeros(nCases, numel(stations_in));
unscaled.My_lbin = zeros(nCases, numel(stations_in));

for iCase = 1:nCases
    zCase = loads.spanwise.distZ{iCase};
    xCase = loads.spanwise.distX{iCase};
    myCase = loads.spanwise.distMy{iCase};

    unscaled.Vz_lb(iCase, :) = interpSpan(ySource_in, zCase.V, stations_in);
    unscaled.Mx_lbin(iCase, :) = 12.0 * interpSpan(ySource_in, zCase.M, stations_in);
    unscaled.Vx_lb(iCase, :) = interpSpan(ySource_in, xCase.Vx, stations_in);
    unscaled.Mz_lbin(iCase, :) = 12.0 * interpSpan(ySource_in, xCase.Mz, stations_in);
    unscaled.My_lbin(iCase, :) = 12.0 * interpSpan(ySource_in, myCase.My_dist, stations_in);
end

spanLoads = struct();
spanLoads.caseNames = caseNames;
spanLoads.nCases = nCases;
spanLoads.stations_in = stations_in;
spanLoads.stations_ft = stations_ft;
spanLoads.nStations = numel(stations_in);
spanLoads.Vx_lb = loadFactor * unscaled.Vx_lb;
spanLoads.Vz_lb = loadFactor * unscaled.Vz_lb;
spanLoads.Mx_lbin = loadFactor * unscaled.Mx_lbin;
spanLoads.Mz_lbin = loadFactor * unscaled.Mz_lbin;
spanLoads.My_lbin = loadFactor * unscaled.My_lbin;
spanLoads.unscaled = unscaled;
spanLoads.loadFactorApplied = loadFactor;
spanLoads.loadFactorPurpose = string(analysisOptions.loadFactorPurpose);
spanLoads.units.y = "in";
spanLoads.units.force = "lb";
spanLoads.units.moment = "lb-in";
spanLoads.units.source_y = "ft";
spanLoads.units.source_moment = "lb-ft";
spanLoads.sourceFields.y = "loads.spanwise.y_ft/y_in";
spanLoads.sourceFields.Vz = "loads.spanwise.distZ{iCase}.V";
spanLoads.sourceFields.Mx = "12 * loads.spanwise.distZ{iCase}.M";
spanLoads.sourceFields.Vx = "loads.spanwise.distX{iCase}.Vx";
spanLoads.sourceFields.Mz = "12 * loads.spanwise.distX{iCase}.Mz";
spanLoads.sourceFields.My = "12 * loads.spanwise.distMy{iCase}.My_dist";
spanLoads.sourceFields.caseNames = "loads.caseNames and loads.spanwise.caseNames";
spanLoads.discoveredFields = discoverSpanwiseFields(loads);
spanLoads.rootCheck = computeRootCheck(loads, unscaled, caseNames);
spanLoads.status = "spanwise_loads_interpolated_at_rib_stations_only";
spanLoads.warnings = makeWarnings(spanLoads.rootCheck);
end

function validateInputs(loads, wing, analysisOptions)
if nargin < 1 || isempty(loads) || ~isstruct(loads)
    error('interpolateSpanwiseLoads:MissingLoads', ...
        'loads struct from criticalLoads is required.');
end
if ~isfield(loads, 'hasSpanwise') || ~loads.hasSpanwise || ...
        ~isfield(loads, 'spanwise')
    error('interpolateSpanwiseLoads:MissingSpanwiseLoads', ...
        'loads must include spanwise data from criticalLoads.');
end
requiredSpanwise = ["y_ft", "y_in", "caseNames", "distZ", "distX", "distMy"];
missingSpanwise = setdiff(requiredSpanwise, string(fieldnames(loads.spanwise)));
if ~isempty(missingSpanwise)
    error('interpolateSpanwiseLoads:MissingSpanwiseFields', ...
        'loads.spanwise is missing required fields: %s', strjoin(missingSpanwise, ', '));
end
if ~isfield(loads, 'root') || ~isfield(loads, 'caseNames') || ~isfield(loads, 'nCases')
    error('interpolateSpanwiseLoads:InvalidLoads', ...
        'loads must include root, caseNames, and nCases fields.');
end
requiredRoot = ["Vx_lb", "Vz_lb", "Mx_lbin", "Mz_lbin", "My_lbin"];
missingRoot = setdiff(requiredRoot, string(fieldnames(loads.root)));
if ~isempty(missingRoot)
    error('interpolateSpanwiseLoads:MissingRootFields', ...
        'loads.root is missing required fields: %s', strjoin(missingRoot, ', '));
end
if numel(loads.spanwise.distZ) ~= loads.nCases || ...
        numel(loads.spanwise.distX) ~= loads.nCases || ...
        numel(loads.spanwise.distMy) ~= loads.nCases
    error('interpolateSpanwiseLoads:SpanwiseCaseCountMismatch', ...
        'distZ, distX, and distMy must contain one cell per load case.');
end
if any(string(loads.spanwise.caseNames(:)) ~= string(loads.caseNames(:)))
    error('interpolateSpanwiseLoads:CaseNameMismatch', ...
        'loads.spanwise.caseNames must match loads.caseNames.');
end

for iCase = 1:loads.nCases
    validateDistFields(loads.spanwise.distZ{iCase}, ["V", "M"], ...
        'loads.spanwise.distZ', iCase);
    validateDistFields(loads.spanwise.distX{iCase}, ["Vx", "Mz"], ...
        'loads.spanwise.distX', iCase);
    validateDistFields(loads.spanwise.distMy{iCase}, "My_dist", ...
        'loads.spanwise.distMy', iCase);
end

if nargin < 2 || isempty(wing) || ~isstruct(wing) || ...
        ~isfield(wing, 'ribs') || ~isfield(wing.ribs, 'spacing_in') || ...
        ~isfield(wing, 'planform') || ~isfield(wing.planform, 'semispan_in')
    error('interpolateSpanwiseLoads:InvalidWing', ...
        'wing.ribs.spacing_in and wing.planform.semispan_in are required.');
end
if ~(isfinite(wing.ribs.spacing_in) && wing.ribs.spacing_in > 0 && ...
        isfinite(wing.planform.semispan_in) && wing.planform.semispan_in > 0)
    error('interpolateSpanwiseLoads:InvalidRibGeometry', ...
        'Rib spacing and semispan must be positive finite values.');
end
if nargin < 3 || isempty(analysisOptions) || ~isstruct(analysisOptions) || ...
        ~isfield(analysisOptions, 'loadFactor') || ...
        ~(isnumeric(analysisOptions.loadFactor) && isscalar(analysisOptions.loadFactor) && ...
        isfinite(analysisOptions.loadFactor) && analysisOptions.loadFactor > 0)
    error('interpolateSpanwiseLoads:InvalidAnalysisOptions', ...
        'analysisOptions.loadFactor must be a positive finite scalar.');
end
if ~isfield(analysisOptions, 'loadFactorPurpose')
    analysisOptions.loadFactorPurpose = "";
end

y_ft = loads.spanwise.y_ft(:);
y_in = loads.spanwise.y_in(:);
if numel(y_ft) < 2 || any(~isfinite(y_ft)) || any(~isfinite(y_in))
    error('interpolateSpanwiseLoads:InvalidYData', ...
        'Spanwise y arrays must contain at least two finite stations.');
end
if any(diff(y_in) <= 0)
    error('interpolateSpanwiseLoads:NonMonotonicY', ...
        'Spanwise y_in data must be strictly increasing.');
end
end

function validateDistFields(s, requiredFields, sourceName, iCase)
if ~isstruct(s)
    error('interpolateSpanwiseLoads:InvalidDistributionCell', ...
        '%s{%d} must be a struct.', sourceName, iCase);
end
requiredFields = string(requiredFields);
missing = setdiff(requiredFields, string(fieldnames(s)));
if ~isempty(missing)
    error('interpolateSpanwiseLoads:MissingDistributionFields', ...
        '%s{%d} is missing fields: %s', sourceName, iCase, strjoin(missing, ', '));
end
for iField = 1:numel(requiredFields)
    values = s.(requiredFields(iField));
    if ~isnumeric(values) || any(~isfinite(values(:)))
        error('interpolateSpanwiseLoads:InvalidDistributionValues', ...
            '%s{%d}.%s must be finite numeric values.', ...
            sourceName, iCase, requiredFields(iField));
    end
end
end

function stations_in = makeRibStations(ribSpacing_in, semispan_in)
stations_in = 0:ribSpacing_in:semispan_in;
if stations_in(end) < semispan_in
    stations_in(end + 1) = semispan_in;
end
stations_in = stations_in(:).';
end

function valuesAtStations = interpSpan(ySource_in, values, stations_in)
values = values(:);
if numel(values) ~= numel(ySource_in)
    error('interpolateSpanwiseLoads:DistributionLengthMismatch', ...
        'Spanwise distribution length does not match y array length.');
end
valuesAtStations = interp1(ySource_in, values, stations_in, 'linear');
if any(~isfinite(valuesAtStations))
    error('interpolateSpanwiseLoads:InterpolationFailed', ...
        'Interpolated spanwise load values must be finite.');
end
end

function rootCheck = computeRootCheck(loads, unscaled, caseNames)
components = ["Vx_lb", "Vz_lb", "Mx_lbin", "Mz_lbin", "My_lbin"];
nCases = loads.nCases;
nComponents = numel(components);
absDiff = zeros(nCases, nComponents);
relDiff = zeros(nCases, nComponents);
interpolatedRoot = zeros(nCases, nComponents);
referenceRoot = zeros(nCases, nComponents);

for iComponent = 1:nComponents
    component = components(iComponent);
    interpolatedRoot(:, iComponent) = unscaled.(component)(:, 1);
    referenceRoot(:, iComponent) = loads.root.(component)(:);
    absDiff(:, iComponent) = abs(interpolatedRoot(:, iComponent) - referenceRoot(:, iComponent));
    denom = max(abs(referenceRoot(:, iComponent)), 1.0);
    relDiff(:, iComponent) = absDiff(:, iComponent) ./ denom;
end

rootCheck = struct();
rootCheck.components = components;
rootCheck.caseNames = caseNames;
rootCheck.interpolatedRoot = interpolatedRoot;
rootCheck.referenceRoot = referenceRoot;
rootCheck.absoluteDifference = absDiff;
rootCheck.relativeDifference = relDiff;
rootCheck.maxAbsoluteDifference = max(absDiff, [], 'all');
rootCheck.maxRelativeDifference = max(relDiff, [], 'all');
[caseIdx, componentIdx] = find(absDiff == rootCheck.maxAbsoluteDifference, 1);
rootCheck.maxAbsoluteDifferenceCase = caseNames(caseIdx);
rootCheck.maxAbsoluteDifferenceComponent = components(componentIdx);
[caseIdxRel, componentIdxRel] = find(relDiff == rootCheck.maxRelativeDifference, 1);
rootCheck.maxRelativeDifferenceCase = caseNames(caseIdxRel);
rootCheck.maxRelativeDifferenceComponent = components(componentIdxRel);
rootCheck.tolerance.absolute = 1e-6;
rootCheck.tolerance.relative = 1e-9;
rootCheck.pass = rootCheck.maxAbsoluteDifference <= rootCheck.tolerance.absolute || ...
    rootCheck.maxRelativeDifference <= rootCheck.tolerance.relative;
end

function discovered = discoverSpanwiseFields(loads)
discovered = struct();
discovered.loads_spanwise_fields = string(fieldnames(loads.spanwise));
discovered.distZ_fields = string(fieldnames(loads.spanwise.distZ{1}));
discovered.distX_fields = string(fieldnames(loads.spanwise.distX{1}));
discovered.distMy_fields = string(fieldnames(loads.spanwise.distMy{1}));
if isfield(loads.spanwise, 'allCases') && ~isempty(loads.spanwise.allCases)
    discovered.allCases_fields = string(fieldnames(loads.spanwise.allCases(1)));
else
    discovered.allCases_fields = strings(0, 1);
end
end

function warnings = makeWarnings(rootCheck)
warnings = strings(0, 1);
if ~rootCheck.pass
    warnings(end + 1, 1) = sprintf( ...
        "Root consistency check exceeds tolerance: max abs %.6g on %s/%s; max rel %.6g on %s/%s.", ...
        rootCheck.maxAbsoluteDifference, rootCheck.maxAbsoluteDifferenceCase, ...
        rootCheck.maxAbsoluteDifferenceComponent, rootCheck.maxRelativeDifference, ...
        rootCheck.maxRelativeDifferenceCase, rootCheck.maxRelativeDifferenceComponent);
end
warnings(end + 1, 1) = ...
    "Interpolates spanwise loads only; no stress, failure, or optimisation logic is run here.";
end
