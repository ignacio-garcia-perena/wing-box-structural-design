function columnBuckling = computeColumnBuckling(wing, section, stresses, ribSpacing_in, K)
% computeColumnBuckling Euler column buckling for booms.
%
% Inputs:
%   wing          - Geometry/material struct from setGeometry.
%   section       - Struct from computeSectionUnitLoads.
%   stresses      - Struct from computeStresses.
%   ribSpacing_in - User-provided unsupported column length between ribs.
%   K             - User-provided effective length factor.
%
% Outputs:
%   columnBuckling - Struct containing Euler column buckling
%                    ratios for spar caps and stringers.
%
% Scope:
%   Euler column buckling only for spar caps and stringers using
%   compressive bending stresses.
%
% Excluded:
%   panel buckling, crippling, combined R-ratios, divergence changes, weight
%   changes, optimisation, and final failure margins.

validateInputs(wing, section, stresses, ribSpacing_in, K);

entries = stresses.bending.entries;
sigma_psi = stresses.bending.sigma_psi;
nCases = stresses.nCases;
nEntries = stresses.nEntries;

E_psi = wing.material.sparCapStringer.E_psi;
stringerRhoMin_in = min(wing.stringers.geometry.rho_xx_in, wing.stringers.geometry.rho_yy_in);
sparCapRhoMin_in = min(wing.sparCaps.section.rho_xx_in, wing.sparCaps.section.rho_yy_in);

entryProperties = repmat(makeEntryProperty(), nEntries, 1);
sigmaCr_psi = zeros(1, nEntries);
for iEntry = 1:nEntries
    entry = entries(iEntry);
    [componentType, rhoMin_in, sectionSource] = classifyEntry(entry, ...
        stringerRhoMin_in, sparCapRhoMin_in, wing);
    sigmaCr_psi(iEntry) = pi^2 * E_psi / (K * ribSpacing_in / rhoMin_in)^2;
    entryProperties(iEntry) = makeEntryProperty(entry, iEntry, ...
        componentType, rhoMin_in, sigmaCr_psi(iEntry), sectionSource);
end

sigmaCompression_psi = max(0, -sigma_psi);
bucklingRatio = sigmaCompression_psi ./ sigmaCr_psi;

if any(~isfinite(bucklingRatio), 'all') || any(bucklingRatio < 0, 'all')
    error('computeColumnBuckling:InvalidBucklingRatio', ...
        'Computed buckling ratios must be finite and nonnegative.');
end

columnBuckling = struct();
columnBuckling.caseNames = stresses.caseNames;
columnBuckling.nCases = nCases;
columnBuckling.nEntries = nEntries;
columnBuckling.entries = entries;
columnBuckling.entryProperties = entryProperties;
columnBuckling.sigma_psi = sigma_psi;
columnBuckling.sigma_compression_psi = sigmaCompression_psi;
columnBuckling.sigma_cr_psi = sigmaCr_psi;
columnBuckling.buckling_ratio = bucklingRatio;
columnBuckling.maxRatioByCase = repmat(makeSummaryEntry(), nCases, 1);
columnBuckling.passByCase = false(nCases, 1);

for iCase = 1:nCases
    [maxRatio, idxEntry] = max(bucklingRatio(iCase, :));
    columnBuckling.maxRatioByCase(iCase) = makeSummaryEntry( ...
        stresses.caseNames(iCase), maxRatio, idxEntry, entries(idxEntry), ...
        entryProperties(idxEntry), sigma_psi(iCase, idxEntry), ...
        sigmaCompression_psi(iCase, idxEntry));
    columnBuckling.passByCase(iCase) = maxRatio <= 1.0;
end

[globalMaxRatio, linearIdx] = max(bucklingRatio(:));
[globalCaseIdx, globalEntryIdx] = ind2sub(size(bucklingRatio), linearIdx);
columnBuckling.globalMaxRatio = makeSummaryEntry( ...
    stresses.caseNames(globalCaseIdx), globalMaxRatio, globalEntryIdx, ...
    entries(globalEntryIdx), entryProperties(globalEntryIdx), ...
    sigma_psi(globalCaseIdx, globalEntryIdx), ...
    sigmaCompression_psi(globalCaseIdx, globalEntryIdx));
columnBuckling.globalMaxRatio.case_index = globalCaseIdx;
columnBuckling.globalPass = globalMaxRatio <= 1.0;

columnBuckling.inputs.ribSpacing_in = ribSpacing_in;
columnBuckling.inputs.K = K;
columnBuckling.material.E_psi = E_psi;
columnBuckling.material.source = wing.material.sparCapStringer.source;
columnBuckling.material.confidence = wing.material.sparCapStringer.confidence;
columnBuckling.properties.stringer.rho_min_in = stringerRhoMin_in;
columnBuckling.properties.stringer.section = wing.stringers.geometry;
columnBuckling.properties.stringer.sigma_cr_psi = ...
    pi^2 * E_psi / (K * ribSpacing_in / stringerRhoMin_in)^2;
columnBuckling.properties.sparCap.rho_min_in = sparCapRhoMin_in;
columnBuckling.properties.sparCap.section = wing.sparCaps.section;
columnBuckling.properties.sparCap.sigma_cr_psi = ...
    pi^2 * E_psi / (K * ribSpacing_in / sparCapRhoMin_in)^2;
columnBuckling.units.stress = "psi";
columnBuckling.units.length = "in";
columnBuckling.units.ratio = "dimensionless";
columnBuckling.method = ...
    "sigma_cr_psi = pi^2 * E_psi / (K * ribSpacing_in / rho_min)^2; compression only.";
columnBuckling.assumptions = [
    "Euler column buckling only."
    "ribSpacing_in and K are user-provided."
    "rho_min is taken from Bruhn section properties."
    "E uses 7075-T6 extrusion properties from the verified project material table."
    "Uses root bending stresses only."
    "Does not use average stress between rib stations yet."
    "Negative sigma is treated as compression; tensile stresses get zero buckling ratio."
    ];
columnBuckling.warnings = [
    "This module checks Euler column buckling only."
    "Combined panel buckling ratios are handled in computePanelBuckling and computeSpanwiseFailureSummary."
    "Crippling interaction and final certification margin stack-up are not included in this module."
    "This is not an optimisation result."
    ];
columnBuckling.status = "euler_column_buckling_screening";
columnBuckling.exclusions = [
    "crippling"
    "crippling interaction"
    "final certification margin stack-up"
    "divergence changes"
    "weight changes"
    "optimisation"
    ];
end

function validateInputs(wing, section, stresses, ribSpacing_in, K)
if nargin < 1 || isempty(wing) || ~isstruct(wing)
    error('computeColumnBuckling:MissingWing', 'wing struct from setGeometry is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computeColumnBuckling:MissingSection', 'section struct from computeSectionUnitLoads is required.');
end
if nargin < 3 || isempty(stresses) || ~isstruct(stresses)
    error('computeColumnBuckling:MissingStresses', 'stresses struct from computeStresses is required.');
end
if nargin < 4 || isempty(ribSpacing_in) || ~(isnumeric(ribSpacing_in) && isscalar(ribSpacing_in)) || ...
        ~(isfinite(ribSpacing_in) && ribSpacing_in > 0)
    error('computeColumnBuckling:InvalidRibSpacing', ...
        'ribSpacing_in must be explicitly provided as a positive finite scalar.');
end
if nargin < 5 || isempty(K) || ~(isnumeric(K) && isscalar(K)) || ...
        ~(isfinite(K) && K > 0)
    error('computeColumnBuckling:InvalidEffectiveLengthFactor', ...
        'K must be explicitly provided as a positive finite scalar effective length factor.');
end

if ~isfield(wing, 'material') || ~isfield(wing.material, 'sparCapStringer') || ...
        ~isfield(wing.material.sparCapStringer, 'E_psi') || ...
        ~(isfinite(wing.material.sparCapStringer.E_psi) && wing.material.sparCapStringer.E_psi > 0)
    error('computeColumnBuckling:InvalidMaterialE', ...
        'wing.material.sparCapStringer.E_psi must be positive and finite.');
end
if ~isfield(wing, 'stringers') || ~isfield(wing.stringers, 'geometry') || ...
        ~isfield(wing.stringers.geometry, 'rho_xx_in') || ...
        ~isfield(wing.stringers.geometry, 'rho_yy_in')
    error('computeColumnBuckling:MissingStringerRho', ...
        'wing.stringers.geometry.rho_xx_in and rho_yy_in are required.');
end
if ~isfield(wing, 'sparCaps') || ~isfield(wing.sparCaps, 'section') || ...
        ~isfield(wing.sparCaps.section, 'rho_xx_in') || ...
        ~isfield(wing.sparCaps.section, 'rho_yy_in')
    error('computeColumnBuckling:MissingSparCapRho', ...
        'wing.sparCaps.section.rho_xx_in and rho_yy_in are required.');
end

rhoValues = [wing.stringers.geometry.rho_xx_in, wing.stringers.geometry.rho_yy_in, ...
    wing.sparCaps.section.rho_xx_in, wing.sparCaps.section.rho_yy_in];
if any(~isfinite(rhoValues)) || any(rhoValues <= 0)
    error('computeColumnBuckling:InvalidRhoValues', ...
        'Stringer and spar-cap rho values must be positive and finite.');
end

if ~isfield(stresses, 'nCases') || ~isfield(stresses, 'nEntries') || ...
        ~isfield(stresses, 'caseNames') || ~isfield(stresses, 'bending') || ...
        ~isfield(stresses.bending, 'entries') || ~isfield(stresses.bending, 'sigma_psi')
    error('computeColumnBuckling:InvalidStresses', ...
        'stresses must include nCases, nEntries, caseNames, bending.entries, and bending.sigma_psi.');
end
if ~isequal(size(stresses.bending.sigma_psi), [stresses.nCases, stresses.nEntries])
    error('computeColumnBuckling:StressMatrixSizeMismatch', ...
        'stresses.bending.sigma_psi must have size nCases x nEntries.');
end
if numel(stresses.bending.entries) ~= stresses.nEntries
    error('computeColumnBuckling:EntryCountMismatch', ...
        'stresses.bending.entries count must match nEntries.');
end
if any(~isfinite(stresses.bending.sigma_psi), 'all')
    error('computeColumnBuckling:InvalidStressMatrix', ...
        'Bending stress matrix must contain finite values.');
end

validGroups = ["sparCaps", "topStringers", "bottomStringers", ...
    "noseTopStringers", "noseBottomStringers"];
entryGroups = [stresses.bending.entries.group];
unknownGroups = setdiff(unique(entryGroups), validGroups);
if ~isempty(unknownGroups)
    error('computeColumnBuckling:UnknownEntryGroups', ...
        'Unsupported bending entry groups for column buckling: %s', ...
        strjoin(unknownGroups, ', '));
end

% The section argument is not otherwise needed yet, but verify that it
% matches the same basic boom count used by stresses.
if ~isfield(section, 'unitBending') || ~isfield(section.unitBending, 'entries') || ...
        numel(section.unitBending.entries) ~= stresses.nEntries
    error('computeColumnBuckling:SectionStressMismatch', ...
        'section.unitBending.entries must exist and match stresses.nEntries.');
end
end

function [componentType, rhoMin_in, sectionSource] = classifyEntry(entry, stringerRhoMin_in, sparCapRhoMin_in, wing)
switch string(entry.group)
    case "sparCaps"
        componentType = "sparCap";
        rhoMin_in = sparCapRhoMin_in;
        sectionSource = wing.sparCaps.section.source;
    case {"topStringers", "bottomStringers", "noseTopStringers", "noseBottomStringers"}
        componentType = "stringer";
        rhoMin_in = stringerRhoMin_in;
        sectionSource = wing.stringers.geometry.source;
    otherwise
        error('computeColumnBuckling:UnsupportedEntryGroup', ...
            'Unsupported entry group: %s', entry.group);
end
end

function property = makeEntryProperty(entry, entryIndex, componentType, rhoMin_in, sigmaCr_psi, sectionSource)
if nargin == 0
    entry = makeEntryPlaceholder();
    entryIndex = NaN;
    componentType = "";
    rhoMin_in = NaN;
    sigmaCr_psi = NaN;
    sectionSource = "";
end

property = struct();
property.entry_index = entryIndex;
property.name = entry.name;
property.group = entry.group;
property.component_type = string(componentType);
property.x_in = entry.x_in;
property.z_in = entry.z_in;
property.rho_min_in = rhoMin_in;
property.sigma_cr_psi = sigmaCr_psi;
property.section_source = string(sectionSource);
end

function summary = makeSummaryEntry(caseName, bucklingRatio, entryIndex, entry, property, sigma_psi, sigmaCompression_psi)
if nargin == 0
    caseName = "";
    bucklingRatio = NaN;
    entryIndex = NaN;
    entry = makeEntryPlaceholder();
    property = makeEntryProperty();
    sigma_psi = NaN;
    sigmaCompression_psi = NaN;
end

summary = struct();
summary.case_name = string(caseName);
summary.buckling_ratio = bucklingRatio;
summary.pass = bucklingRatio <= 1.0;
summary.entry_index = entryIndex;
summary.group = entry.group;
summary.name = entry.name;
summary.component_type = property.component_type;
summary.x_in = entry.x_in;
summary.z_in = entry.z_in;
summary.rho_min_in = property.rho_min_in;
summary.sigma_cr_psi = property.sigma_cr_psi;
summary.sigma_psi = sigma_psi;
summary.sigma_compression_psi = sigmaCompression_psi;
end

function entry = makeEntryPlaceholder()
entry = struct();
entry.name = "";
entry.group = "";
entry.x_in = NaN;
entry.z_in = NaN;
end

