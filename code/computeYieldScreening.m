function yieldScreening = computeYieldScreening(vonMises, wing)
% computeYieldScreening Panel yield screening.
%
% Inputs:
%   vonMises - Struct from computeVonMises containing panel
%              von Mises stresses.
%   wing     - Struct from setGeometry containing material properties.
%
% Outputs:
%   yieldScreening - Struct containing yield ratios and
%                    pass/fail flags for every load case and physical panel.
%
% Scope:
%   Panel yield screening for the current section model.
%
% Method:
%   yield_ratio = sigma_vm_psi / yield_allowable_psi
%   yield_allowable_psi = min(wing.material.skinWeb.Fty_psi,
%                             wing.material.skinWeb.Fcy_psi)
%
% Excluded:
%   buckling, divergence, weight, optimisation, and final failure margins.

validateInputs(vonMises, wing);

sigmaVm_psi = vonMises.panel.sigma_vm_psi;
skinWeb = wing.material.skinWeb;
yieldAllowable_psi = min(skinWeb.Fty_psi, skinWeb.Fcy_psi);

yieldRatio = sigmaVm_psi ./ yieldAllowable_psi;

if any(~isfinite(yieldRatio), 'all')
    error('computeYieldScreening:InvalidYieldRatio', ...
        'Computed yield ratio matrix contains non-finite values.');
end
if any(yieldRatio < 0, 'all')
    error('computeYieldScreening:NegativeYieldRatio', ...
        'Yield ratios should be nonnegative because von Mises stresses are nonnegative.');
end

nCases = vonMises.nCases;
nPanels = vonMises.nPanels;
panels = vonMises.panels;

yieldScreening = struct();
yieldScreening.caseNames = vonMises.caseNames;
yieldScreening.nCases = nCases;
yieldScreening.nPanels = nPanels;
yieldScreening.panels = panels;
yieldScreening.sigma_vm_psi = sigmaVm_psi;
yieldScreening.yield_allowable_psi = yieldAllowable_psi;
yieldScreening.yield_ratio = yieldRatio;
yieldScreening.maxRatioByCase = repmat(makeSummaryEntry(), nCases, 1);
yieldScreening.passByCase = false(nCases, 1);

for iCase = 1:nCases
    [maxRatio, idxPanel] = max(yieldRatio(iCase, :));
    yieldScreening.maxRatioByCase(iCase) = makeSummaryEntry( ...
        vonMises.caseNames(iCase), maxRatio, idxPanel, panels(idxPanel), ...
        sigmaVm_psi(iCase, idxPanel), yieldAllowable_psi);
    yieldScreening.passByCase(iCase) = maxRatio <= 1.0;
end

[globalMaxRatio, linearIdx] = max(yieldRatio(:));
[globalCaseIdx, globalPanelIdx] = ind2sub(size(yieldRatio), linearIdx);
yieldScreening.globalMaxRatio = makeSummaryEntry( ...
    vonMises.caseNames(globalCaseIdx), globalMaxRatio, globalPanelIdx, ...
    panels(globalPanelIdx), sigmaVm_psi(globalCaseIdx, globalPanelIdx), ...
    yieldAllowable_psi);
yieldScreening.globalMaxRatio.case_index = globalCaseIdx;
yieldScreening.globalPass = globalMaxRatio <= 1.0;

yieldScreening.materialUsed.alloy = skinWeb.name;
yieldScreening.materialUsed.temper = skinWeb.temper;
yieldScreening.materialUsed.product_form = skinWeb.product_form;
yieldScreening.materialUsed.Fty_psi = skinWeb.Fty_psi;
yieldScreening.materialUsed.Fcy_psi = skinWeb.Fcy_psi;
yieldScreening.materialUsed.Ftu_psi = skinWeb.Ftu_psi;
yieldScreening.materialUsed.source = skinWeb.source;
yieldScreening.materialUsed.confidence = skinWeb.confidence;
yieldScreening.materialUsed.allowable_rule = "min(Fty_psi, Fcy_psi)";

yieldScreening.units.stress = "psi";
yieldScreening.units.yieldRatio = "dimensionless";
yieldScreening.method = ...
    "yield_ratio = sigma_vm_psi / min(wing.material.skinWeb.Fty_psi, wing.material.skinWeb.Fcy_psi)";
yieldScreening.warnings = [
    "Panel yield screening for the current section model."
    "Material properties are read from the verified project material table."
    "Panel sigma is approximated from adjacent/nearest boom stress."
    "Optional panel idealisation into boom areas is handled upstream in computeSectionUnitLoads when enabled."
    "Detailed plate stress recovery within each panel is not implemented."
    "Spar cap/stringer material-point yield is not separately evaluated."
    "No factor of safety has been applied."
    ];
yieldScreening.status = "panel_yield_screening";
yieldScreening.exclusions = [
    "buckling"
    "divergence"
    "weight"
    "optimisation"
    "final failure margins"
    ];
end

function validateInputs(vonMises, wing)
if nargin < 1 || isempty(vonMises) || ~isstruct(vonMises)
    error('computeYieldScreening:MissingVonMises', ...
        'vonMises struct from computeVonMises is required.');
end
if nargin < 2 || isempty(wing) || ~isstruct(wing)
    error('computeYieldScreening:MissingWing', ...
        'wing struct from setGeometry is required.');
end

if ~isfield(vonMises, 'nCases') || ~isfield(vonMises, 'nPanels') || ...
        ~isfield(vonMises, 'caseNames') || ~isfield(vonMises, 'panels') || ...
        ~isfield(vonMises, 'panel') || ~isfield(vonMises.panel, 'sigma_vm_psi')
    error('computeYieldScreening:InvalidVonMises', ...
        'vonMises must include nCases, nPanels, caseNames, panels, and panel.sigma_vm_psi.');
end

sigmaVm_psi = vonMises.panel.sigma_vm_psi;
if ~isequal(size(sigmaVm_psi), [vonMises.nCases, vonMises.nPanels])
    error('computeYieldScreening:VonMisesMatrixSizeMismatch', ...
        'vonMises.panel.sigma_vm_psi must have size nCases x nPanels.');
end
if any(~isfinite(sigmaVm_psi), 'all') || any(sigmaVm_psi < 0, 'all')
    error('computeYieldScreening:InvalidVonMisesMatrix', ...
        'von Mises stresses must be finite and nonnegative.');
end
if numel(vonMises.caseNames) ~= vonMises.nCases
    error('computeYieldScreening:CaseNameSizeMismatch', ...
        'vonMises.caseNames must have one name per load case.');
end
if numel(vonMises.panels) ~= vonMises.nPanels
    error('computeYieldScreening:PanelCountMismatch', ...
        'vonMises.panels must contain one panel per von Mises matrix column.');
end

if ~isfield(wing, 'material') || ~isfield(wing.material, 'skinWeb')
    error('computeYieldScreening:MissingSkinWebMaterial', ...
        'wing.material.skinWeb is required for panel yield screening.');
end

skinWeb = wing.material.skinWeb;
requiredFields = ["name", "temper", "product_form", "Fty_psi", ...
    "Fcy_psi", "Ftu_psi", "source", "confidence"];
materialFields = string(fieldnames(skinWeb));
missingFields = setdiff(requiredFields, materialFields);
if ~isempty(missingFields)
    error('computeYieldScreening:MissingMaterialFields', ...
        'wing.material.skinWeb is missing required fields: %s', ...
        strjoin(missingFields, ', '));
end

if ~(isfinite(skinWeb.Fty_psi) && skinWeb.Fty_psi > 0 && ...
        isfinite(skinWeb.Fcy_psi) && skinWeb.Fcy_psi > 0)
    error('computeYieldScreening:InvalidYieldAllowables', ...
        'wing.material.skinWeb.Fty_psi and Fcy_psi must be positive finite values.');
end
if ~(isfinite(skinWeb.Ftu_psi) && skinWeb.Ftu_psi > 0)
    error('computeYieldScreening:InvalidUltimateAllowable', ...
        'wing.material.skinWeb.Ftu_psi must be positive and finite for metadata.');
end
end

function entry = makeSummaryEntry(caseName, yieldRatio, panelIndex, panel, sigmaVm_psi, yieldAllowable_psi)
if nargin == 0
    caseName = "";
    yieldRatio = NaN;
    panelIndex = NaN;
    panel = makePanelPlaceholder();
    sigmaVm_psi = NaN;
    yieldAllowable_psi = NaN;
end

entry = struct();
entry.case_name = string(caseName);
entry.yield_ratio = yieldRatio;
entry.pass = yieldRatio <= 1.0;
entry.panel_index = panelIndex;
entry.panel_group = panel.group;
entry.panel_name = panel.name;
entry.xStart_in = panel.xStart_in;
entry.zStart_in = panel.zStart_in;
entry.xEnd_in = panel.xEnd_in;
entry.zEnd_in = panel.zEnd_in;
entry.thickness_in = panel.thickness_in;
entry.sigma_vm_psi = sigmaVm_psi;
entry.yield_allowable_psi = yieldAllowable_psi;
end

function panel = makePanelPlaceholder()
panel = struct();
panel.name = "";
panel.group = "";
panel.xStart_in = NaN;
panel.zStart_in = NaN;
panel.xEnd_in = NaN;
panel.zEnd_in = NaN;
panel.thickness_in = NaN;
end

