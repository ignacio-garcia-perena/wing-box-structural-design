function stresses = computeStresses(loads, section)
% computeStresses Scale unit bending coefficients to saved root moments.
%
% Inputs:
%   loads   - Struct from criticalLoads.
%   section - Struct from computeSectionUnitLoads.
%
% Outputs:
%   stresses - Struct containing bending normal stress results for every
%              saved load case and every boom/stringer coefficient entry.
%
% Scope:
%   Bending normal stress only.
%
% Excluded:
%   shear stress, shear flow, torsion, von Mises, buckling, failure
%   margins, yield checks, divergence, weight, and optimisation.

validateInputs(loads, section);

entries = section.unitBending.entries;
nCases = loads.nCases;
nEntries = numel(entries);

Mx_lbin = loads.root.Mx_lbin(:);
Mz_lbin = loads.root.Mz_lbin(:);
sigY_Mx = [entries.sigY_Mx_per_lbin];
sigY_Mz = [entries.sigY_Mz_per_lbin];

sigma_psi = Mx_lbin * sigY_Mx + Mz_lbin * sigY_Mz;

if any(~isfinite(sigma_psi), 'all')
    error('computeStresses:InvalidStressMatrix', ...
        'Computed bending stress matrix contains non-finite values.');
end

stresses = struct();
stresses.caseNames = loads.caseNames;
stresses.nCases = nCases;
stresses.nEntries = nEntries;
stresses.bending.entries = entries;
stresses.bending.sigma_psi = sigma_psi;
stresses.bending.maxTensionByCase = repmat(makeSummaryEntry(), nCases, 1);
stresses.bending.maxCompressionByCase = repmat(makeSummaryEntry(), nCases, 1);
stresses.bending.maxAbsByCase = repmat(makeSummaryEntry(), nCases, 1);

for iCase = 1:nCases
    caseSigma = sigma_psi(iCase, :);

    [maxTension, idxTension] = max(caseSigma);
    [maxCompression, idxCompression] = min(caseSigma);
    [maxAbs, idxAbs] = max(abs(caseSigma));

    stresses.bending.maxTensionByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), maxTension, idxTension, entries(idxTension), ...
        Mx_lbin(iCase), Mz_lbin(iCase));
    stresses.bending.maxCompressionByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), maxCompression, idxCompression, entries(idxCompression), ...
        Mx_lbin(iCase), Mz_lbin(iCase));
    stresses.bending.maxAbsByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), caseSigma(idxAbs), idxAbs, entries(idxAbs), ...
        Mx_lbin(iCase), Mz_lbin(iCase), maxAbs);
end

[globalMaxAbs, linearIdx] = max(abs(sigma_psi(:)));
[globalCaseIdx, globalEntryIdx] = ind2sub(size(sigma_psi), linearIdx);
stresses.bending.globalMaxAbs = makeSummaryEntry( ...
    loads.caseNames(globalCaseIdx), sigma_psi(globalCaseIdx, globalEntryIdx), ...
    globalEntryIdx, entries(globalEntryIdx), ...
    Mx_lbin(globalCaseIdx), Mz_lbin(globalCaseIdx), globalMaxAbs);
stresses.bending.globalMaxAbs.case_index = globalCaseIdx;

stresses.bending.units.stress = "psi";
stresses.bending.units.moment = "lb-in";
stresses.bending.method = ...
    "sigma_psi = Mx_lbin*sigY_Mx_per_lbin + Mz_lbin*sigY_Mz_per_lbin";
stresses.bending.signConvention = ...
    "Positive sigma is reported as tension and negative sigma as compression for summary and screening.";

stresses.status = "bending_normal_stress_recovery";
stresses.downstreamChecks = [
    "shear stress"
    "shear flow"
    "torsion"
    "von Mises"
    "buckling"
    "failure margins"
    "yield checks"
    "divergence"
    "weight"
    "optimisation"
    ];
end

function validateInputs(loads, section)
if nargin < 1 || isempty(loads) || ~isstruct(loads)
    error('computeStresses:MissingLoads', 'loads struct from criticalLoads is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computeStresses:MissingSection', 'section struct from computeSectionUnitLoads is required.');
end
if ~isfield(loads, 'nCases') || ~isfield(loads, 'caseNames') || ~isfield(loads, 'root')
    error('computeStresses:InvalidLoads', ...
        'loads must include nCases, caseNames, and root fields.');
end
if ~isfield(loads.root, 'Mx_lbin') || ~isfield(loads.root, 'Mz_lbin')
    error('computeStresses:MissingRootMoments', ...
        'loads.root.Mx_lbin and loads.root.Mz_lbin are required.');
end
if numel(loads.root.Mx_lbin) ~= loads.nCases || numel(loads.root.Mz_lbin) ~= loads.nCases
    error('computeStresses:RootMomentSizeMismatch', ...
        'Root moment arrays must have one value per load case.');
end
if numel(loads.caseNames) ~= loads.nCases
    error('computeStresses:CaseNameSizeMismatch', ...
        'loads.caseNames must have one name per load case.');
end
if any(~isfinite(loads.root.Mx_lbin)) || any(~isfinite(loads.root.Mz_lbin))
    error('computeStresses:InvalidRootMoments', ...
        'Root moments must be finite lb-in values.');
end

if ~isfield(section, 'unitBending') || ~isfield(section.unitBending, 'entries')
    error('computeStresses:MissingUnitBending', ...
        'section.unitBending.entries is required.');
end
entries = section.unitBending.entries;
if isempty(entries)
    error('computeStresses:EmptyUnitBending', ...
        'section.unitBending.entries must not be empty.');
end

requiredEntryFields = ["name", "group", "x_in", "z_in", "dx_in", "dz_in", ...
    "sigY_Mx_per_lbin", "sigY_Mz_per_lbin"];
entryFieldNames = string(fieldnames(entries));
missingEntryFields = setdiff(requiredEntryFields, entryFieldNames);
if ~isempty(missingEntryFields)
    error('computeStresses:MissingEntryFields', ...
        'Unit bending entries are missing required fields: %s', ...
        strjoin(missingEntryFields, ', '));
end

if any(~isfinite([entries.sigY_Mx_per_lbin])) || any(~isfinite([entries.sigY_Mz_per_lbin]))
    error('computeStresses:InvalidUnitCoefficients', ...
        'Unit bending coefficients must be finite.');
end
if isfield(section.unitBending, 'units') && isfield(section.unitBending.units, 'coefficient')
    expectedUnits = "psi per lb-in";
    if string(section.unitBending.units.coefficient) ~= expectedUnits
        error('computeStresses:UnexpectedCoefficientUnits', ...
            'Expected unit bending coefficients in "%s".', expectedUnits);
    end
end
end

function entry = makeSummaryEntry(caseName, sigma_psi, entryIndex, coeffEntry, Mx_lbin, Mz_lbin, absSigma_psi)
if nargin == 0
    caseName = "";
    sigma_psi = NaN;
    entryIndex = NaN;
    coeffEntry = makeCoeffPlaceholder();
    Mx_lbin = NaN;
    Mz_lbin = NaN;
    absSigma_psi = NaN;
elseif nargin < 7
    absSigma_psi = abs(sigma_psi);
end

entry = struct();
entry.case_name = string(caseName);
entry.sigma_psi = sigma_psi;
entry.abs_sigma_psi = absSigma_psi;
entry.entry_index = entryIndex;
entry.group = coeffEntry.group;
entry.name = coeffEntry.name;
entry.x_in = coeffEntry.x_in;
entry.z_in = coeffEntry.z_in;
entry.Mx_lbin = Mx_lbin;
entry.Mz_lbin = Mz_lbin;
end

function coeffEntry = makeCoeffPlaceholder()
coeffEntry = struct();
coeffEntry.name = "";
coeffEntry.group = "";
coeffEntry.x_in = NaN;
coeffEntry.z_in = NaN;
end
