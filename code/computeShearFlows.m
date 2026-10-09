function shearFlows = computeShearFlows(loads, section)
% computeShearFlows Scale unit shear-flow coefficients to saved loads.
%
% Inputs:
%   loads   - Struct from criticalLoads.
%   section - Struct from computeSectionUnitLoads.
%
% Outputs:
%   shearFlows - Struct containing shear-flow results for every saved
%                load case and every physical wing-box panel.
%
% Scope:
%   Shear flow only.
%
% Excluded:
%   shear stress, von Mises, buckling, yield, divergence, weight, and
%   optimisation.

validateInputs(loads, section);

panels = section.unitShear.stage3C.unitTotals.physicalPanels;
nCases = loads.nCases;
nPanels = numel(panels);

Vx_lb = loads.root.Vx_lb(:);
Vz_lb = loads.root.Vz_lb(:);
My_lbin = loads.root.My_lbin(:);

qtot_X = [panels.qtot_X_per_lb];
qtot_Z = [panels.qtot_Z_per_lb];
qtot_Y = [panels.qtot_Y_per_lbin];

q_lb_per_in = Vx_lb * qtot_X + Vz_lb * qtot_Z + My_lbin * qtot_Y;

if any(~isfinite(q_lb_per_in), 'all')
    error('computeShearFlows:InvalidShearFlowMatrix', ...
        'Computed shear-flow matrix contains non-finite values.');
end

shearFlows = struct();
shearFlows.caseNames = loads.caseNames;
shearFlows.nCases = nCases;
shearFlows.nPanels = nPanels;
shearFlows.panels = panels;
shearFlows.shear.q_lb_per_in = q_lb_per_in;
shearFlows.shear.maxPositiveByCase = repmat(makeSummaryEntry(), nCases, 1);
shearFlows.shear.maxNegativeByCase = repmat(makeSummaryEntry(), nCases, 1);
shearFlows.shear.maxAbsByCase = repmat(makeSummaryEntry(), nCases, 1);

for iCase = 1:nCases
    caseQ = q_lb_per_in(iCase, :);

    [maxPositive, idxPositive] = max(caseQ);
    [maxNegative, idxNegative] = min(caseQ);
    [maxAbs, idxAbs] = max(abs(caseQ));

    shearFlows.shear.maxPositiveByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), maxPositive, idxPositive, panels(idxPositive), ...
        Vx_lb(iCase), Vz_lb(iCase), My_lbin(iCase));
    shearFlows.shear.maxNegativeByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), maxNegative, idxNegative, panels(idxNegative), ...
        Vx_lb(iCase), Vz_lb(iCase), My_lbin(iCase));
    shearFlows.shear.maxAbsByCase(iCase) = makeSummaryEntry( ...
        loads.caseNames(iCase), caseQ(idxAbs), idxAbs, panels(idxAbs), ...
        Vx_lb(iCase), Vz_lb(iCase), My_lbin(iCase), maxAbs);
end

[globalMaxAbs, linearIdx] = max(abs(q_lb_per_in(:)));
[globalCaseIdx, globalPanelIdx] = ind2sub(size(q_lb_per_in), linearIdx);
shearFlows.shear.globalMaxAbs = makeSummaryEntry( ...
    loads.caseNames(globalCaseIdx), q_lb_per_in(globalCaseIdx, globalPanelIdx), ...
    globalPanelIdx, panels(globalPanelIdx), ...
    Vx_lb(globalCaseIdx), Vz_lb(globalCaseIdx), My_lbin(globalCaseIdx), globalMaxAbs);
shearFlows.shear.globalMaxAbs.case_index = globalCaseIdx;

shearFlows.units.force = "lb";
shearFlows.units.moment = "lb-in";
shearFlows.units.shearFlow = "lb/in";
shearFlows.method = ...
    "q_lb_per_in = Vx_lb*qtot_X_per_lb + Vz_lb*qtot_Z_per_lb + My_lbin*qtot_Y_per_lbin";
shearFlows.signConvention = ...
    "Sign follows panel traversal and saved load signs; use absolute panel shear for sizing margins.";
shearFlows.status = "shear_flow_recovery";
shearFlows.downstreamChecks = [
    "shear stress"
    "von Mises"
    "buckling"
    "yield"
    "divergence"
    "weight"
    "optimisation"
    ];
end

function validateInputs(loads, section)
if nargin < 1 || isempty(loads) || ~isstruct(loads)
    error('computeShearFlows:MissingLoads', 'loads struct from criticalLoads is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computeShearFlows:MissingSection', 'section struct from computeSectionUnitLoads is required.');
end
if ~isfield(loads, 'nCases') || ~isfield(loads, 'caseNames') || ~isfield(loads, 'root')
    error('computeShearFlows:InvalidLoads', ...
        'loads must include nCases, caseNames, and root fields.');
end

requiredRootFields = ["Vx_lb", "Vz_lb", "My_lbin"];
rootFieldNames = string(fieldnames(loads.root));
missingRootFields = setdiff(requiredRootFields, rootFieldNames);
if ~isempty(missingRootFields)
    error('computeShearFlows:MissingRootLoads', ...
        'loads.root is missing required fields: %s', strjoin(missingRootFields, ', '));
end

if numel(loads.root.Vx_lb) ~= loads.nCases || ...
        numel(loads.root.Vz_lb) ~= loads.nCases || ...
        numel(loads.root.My_lbin) ~= loads.nCases
    error('computeShearFlows:RootLoadSizeMismatch', ...
        'Vx, Vz, and My root load arrays must have one value per load case.');
end
if numel(loads.caseNames) ~= loads.nCases
    error('computeShearFlows:CaseNameSizeMismatch', ...
        'loads.caseNames must have one name per load case.');
end
if any(~isfinite(loads.root.Vx_lb)) || any(~isfinite(loads.root.Vz_lb)) || ...
        any(~isfinite(loads.root.My_lbin))
    error('computeShearFlows:InvalidRootLoads', ...
        'Root Vx, Vz, and My values must be finite.');
end

if ~isfield(loads, 'units') || ~isfield(loads.units, 'force') || ...
        ~isfield(loads.units, 'moment')
    error('computeShearFlows:MissingLoadUnits', ...
        'loads.units.force and loads.units.moment are required.');
end
if string(loads.units.force) ~= "lb"
    error('computeShearFlows:UnexpectedForceUnits', ...
        'Expected loads.units.force to be "lb".');
end
if string(loads.units.moment) ~= "lb-in"
    error('computeShearFlows:UnexpectedMomentUnits', ...
        'Expected loads.units.moment to be "lb-in".');
end

if ~isfield(section, 'unitShear') || ~isfield(section.unitShear, 'stage3C') || ...
        ~isfield(section.unitShear.stage3C, 'unitTotals') || ...
        ~isfield(section.unitShear.stage3C.unitTotals, 'physicalPanels')
    error('computeShearFlows:MissingUnitTotals', ...
        'section.unitShear.stage3C.unitTotals.physicalPanels is required.');
end

panels = section.unitShear.stage3C.unitTotals.physicalPanels;
if isempty(panels)
    error('computeShearFlows:EmptyPhysicalPanels', ...
        'section.unitShear.stage3C.unitTotals.physicalPanels must not be empty.');
end

requiredPanelFields = ["name", "group", "xStart_in", "zStart_in", ...
    "xEnd_in", "zEnd_in", "thickness_in", ...
    "qtot_X_per_lb", "qtot_Z_per_lb", "qtot_Y_per_lbin"];
panelFieldNames = string(fieldnames(panels));
missingPanelFields = setdiff(requiredPanelFields, panelFieldNames);
if ~isempty(missingPanelFields)
    error('computeShearFlows:MissingPanelFields', ...
        'Unit shear-flow panels are missing required fields: %s', ...
        strjoin(missingPanelFields, ', '));
end

if any(~isfinite([panels.qtot_X_per_lb])) || ...
        any(~isfinite([panels.qtot_Z_per_lb])) || ...
        any(~isfinite([panels.qtot_Y_per_lbin]))
    error('computeShearFlows:InvalidUnitCoefficients', ...
        'All unit shear-flow coefficients must be finite.');
end
if any(~isfinite([panels.xStart_in])) || any(~isfinite([panels.zStart_in])) || ...
        any(~isfinite([panels.xEnd_in])) || any(~isfinite([panels.zEnd_in]))
    error('computeShearFlows:InvalidPanelCoordinates', ...
        'All panel coordinates must be finite.');
end

unitTotals = section.unitShear.stage3C.unitTotals;
if isfield(unitTotals, 'units')
    if ~isfield(unitTotals.units, 'qtot_X') || string(unitTotals.units.qtot_X) ~= "lb/in per lb"
        error('computeShearFlows:UnexpectedQtotXUnits', ...
            'Expected qtot_X units to be "lb/in per lb".');
    end
    if ~isfield(unitTotals.units, 'qtot_Z') || string(unitTotals.units.qtot_Z) ~= "lb/in per lb"
        error('computeShearFlows:UnexpectedQtotZUnits', ...
            'Expected qtot_Z units to be "lb/in per lb".');
    end
    if ~isfield(unitTotals.units, 'qtot_Y') || string(unitTotals.units.qtot_Y) ~= "lb/in per lb-in"
        error('computeShearFlows:UnexpectedQtotYUnits', ...
            'Expected qtot_Y units to be "lb/in per lb-in".');
    end
end
end

function entry = makeSummaryEntry(caseName, q_lb_per_in, panelIndex, panel, Vx_lb, Vz_lb, My_lbin, absQ_lb_per_in)
if nargin == 0
    caseName = "";
    q_lb_per_in = NaN;
    panelIndex = NaN;
    panel = makePanelPlaceholder();
    Vx_lb = NaN;
    Vz_lb = NaN;
    My_lbin = NaN;
    absQ_lb_per_in = NaN;
elseif nargin < 8
    absQ_lb_per_in = abs(q_lb_per_in);
end

entry = struct();
entry.case_name = string(caseName);
entry.q_lb_per_in = q_lb_per_in;
entry.abs_q_lb_per_in = absQ_lb_per_in;
entry.panel_index = panelIndex;
entry.group = panel.group;
entry.name = panel.name;
entry.xStart_in = panel.xStart_in;
entry.zStart_in = panel.zStart_in;
entry.xEnd_in = panel.xEnd_in;
entry.zEnd_in = panel.zEnd_in;
entry.thickness_in = panel.thickness_in;
entry.Vx_lb = Vx_lb;
entry.Vz_lb = Vz_lb;
entry.My_lbin = My_lbin;
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
