function shearStress = computePanelShearStress(shearFlows)
% computePanelShearStress Convert panel shear flow to panel shear stress.
%
% Inputs:
%   shearFlows - Struct from computeShearFlows.
%
% Outputs:
%   shearStress - Struct containing panel shear stress results for every
%                 saved load case and every physical wing-box panel.
%
% Scope:
%   Panel shear stress only.
%
% Excluded:
%   von Mises, buckling, yield, divergence, failure margins, weight, and
%   optimisation.

validateInputs(shearFlows);

panels = shearFlows.panels;
nCases = shearFlows.nCases;
nPanels = shearFlows.nPanels;

q_lb_per_in = shearFlows.shear.q_lb_per_in;
thickness_in = [panels.thickness_in];

tau_psi = q_lb_per_in ./ thickness_in;

if any(~isfinite(tau_psi), 'all')
    error('computePanelShearStress:InvalidShearStressMatrix', ...
        'Computed panel shear stress matrix contains non-finite values.');
end

shearStress = struct();
shearStress.caseNames = shearFlows.caseNames;
shearStress.nCases = nCases;
shearStress.nPanels = nPanels;
shearStress.panels = panels;
shearStress.panelShear.tau_psi = tau_psi;
shearStress.panelShear.maxPositiveByCase = repmat(makeSummaryEntry(), nCases, 1);
shearStress.panelShear.maxNegativeByCase = repmat(makeSummaryEntry(), nCases, 1);
shearStress.panelShear.maxAbsByCase = repmat(makeSummaryEntry(), nCases, 1);

for iCase = 1:nCases
    caseTau = tau_psi(iCase, :);

    [maxPositive, idxPositive] = max(caseTau);
    [maxNegative, idxNegative] = min(caseTau);
    [maxAbs, idxAbs] = max(abs(caseTau));

    shearStress.panelShear.maxPositiveByCase(iCase) = makeSummaryEntry( ...
        shearFlows.caseNames(iCase), maxPositive, idxPositive, panels(idxPositive), ...
        shearFlows.shear.q_lb_per_in(iCase, idxPositive));
    shearStress.panelShear.maxNegativeByCase(iCase) = makeSummaryEntry( ...
        shearFlows.caseNames(iCase), maxNegative, idxNegative, panels(idxNegative), ...
        shearFlows.shear.q_lb_per_in(iCase, idxNegative));
    shearStress.panelShear.maxAbsByCase(iCase) = makeSummaryEntry( ...
        shearFlows.caseNames(iCase), caseTau(idxAbs), idxAbs, panels(idxAbs), ...
        shearFlows.shear.q_lb_per_in(iCase, idxAbs), maxAbs);
end

[globalMaxAbs, linearIdx] = max(abs(tau_psi(:)));
[globalCaseIdx, globalPanelIdx] = ind2sub(size(tau_psi), linearIdx);
shearStress.panelShear.globalMaxAbs = makeSummaryEntry( ...
    shearFlows.caseNames(globalCaseIdx), tau_psi(globalCaseIdx, globalPanelIdx), ...
    globalPanelIdx, panels(globalPanelIdx), ...
    shearFlows.shear.q_lb_per_in(globalCaseIdx, globalPanelIdx), globalMaxAbs);
shearStress.panelShear.globalMaxAbs.case_index = globalCaseIdx;

shearStress.units.stress = "psi";
shearStress.units.shearFlow = "lb/in";
shearStress.units.thickness = "in";
shearStress.method = "tau_psi = q_lb_per_in / thickness_in";
shearStress.signConvention = ...
    "Sign follows panel traversal and saved load signs; absolute values are used for sizing checks where appropriate.";
shearStress.status = "panel_shear_stress_recovery";
shearStress.downstreamChecks = [
    "von Mises"
    "buckling"
    "yield"
    "divergence"
    "failure margins"
    "weight"
    "optimisation"
    ];
end

function validateInputs(shearFlows)
if nargin < 1 || isempty(shearFlows) || ~isstruct(shearFlows)
    error('computePanelShearStress:MissingShearFlows', ...
        'shearFlows struct from computeShearFlows is required.');
end
if ~isfield(shearFlows, 'nCases') || ~isfield(shearFlows, 'nPanels') || ...
        ~isfield(shearFlows, 'caseNames') || ~isfield(shearFlows, 'panels') || ...
        ~isfield(shearFlows, 'shear')
    error('computePanelShearStress:InvalidShearFlows', ...
        'shearFlows must include nCases, nPanels, caseNames, panels, and shear fields.');
end
if ~isfield(shearFlows.shear, 'q_lb_per_in')
    error('computePanelShearStress:MissingShearFlowMatrix', ...
        'shearFlows.shear.q_lb_per_in is required.');
end

q_lb_per_in = shearFlows.shear.q_lb_per_in;
if ~isequal(size(q_lb_per_in), [shearFlows.nCases, shearFlows.nPanels])
    error('computePanelShearStress:ShearFlowMatrixSizeMismatch', ...
        'q_lb_per_in must have size nCases x nPanels.');
end
if any(~isfinite(q_lb_per_in), 'all')
    error('computePanelShearStress:InvalidShearFlowMatrix', ...
        'q_lb_per_in must contain finite values.');
end
if numel(shearFlows.caseNames) ~= shearFlows.nCases
    error('computePanelShearStress:CaseNameSizeMismatch', ...
        'shearFlows.caseNames must have one name per load case.');
end

panels = shearFlows.panels;
if isempty(panels) || numel(panels) ~= shearFlows.nPanels
    error('computePanelShearStress:PanelCountMismatch', ...
        'shearFlows.panels must contain one panel per shear-flow matrix column.');
end

requiredPanelFields = ["name", "group", "xStart_in", "zStart_in", ...
    "xEnd_in", "zEnd_in", "thickness_in"];
panelFieldNames = string(fieldnames(panels));
missingPanelFields = setdiff(requiredPanelFields, panelFieldNames);
if ~isempty(missingPanelFields)
    error('computePanelShearStress:MissingPanelFields', ...
        'Panel metadata is missing required fields: %s', ...
        strjoin(missingPanelFields, ', '));
end

thickness_in = [panels.thickness_in];
if any(~isfinite(thickness_in)) || any(thickness_in <= 0)
    error('computePanelShearStress:InvalidPanelThickness', ...
        'All panel thickness_in values must be positive and finite.');
end

if ~isfield(shearFlows, 'units') || ~isfield(shearFlows.units, 'shearFlow')
    error('computePanelShearStress:MissingShearFlowUnits', ...
        'shearFlows.units.shearFlow is required.');
end
if string(shearFlows.units.shearFlow) ~= "lb/in"
    error('computePanelShearStress:UnexpectedShearFlowUnits', ...
        'Expected shearFlows.units.shearFlow to be "lb/in".');
end
end

function entry = makeSummaryEntry(caseName, tau_psi, panelIndex, panel, q_lb_per_in, absTau_psi)
if nargin == 0
    caseName = "";
    tau_psi = NaN;
    panelIndex = NaN;
    panel = makePanelPlaceholder();
    q_lb_per_in = NaN;
    absTau_psi = NaN;
elseif nargin < 6
    absTau_psi = abs(tau_psi);
end

entry = struct();
entry.case_name = string(caseName);
entry.tau_psi = tau_psi;
entry.abs_tau_psi = absTau_psi;
entry.panel_index = panelIndex;
entry.group = panel.group;
entry.name = panel.name;
entry.xStart_in = panel.xStart_in;
entry.zStart_in = panel.zStart_in;
entry.xEnd_in = panel.xEnd_in;
entry.zEnd_in = panel.zEnd_in;
entry.thickness_in = panel.thickness_in;
entry.q_lb_per_in = q_lb_per_in;
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
