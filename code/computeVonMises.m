function vonMises = computeVonMises(stresses, shearStress)
% computeVonMises Compute panel von Mises stresses.
%
% Inputs:
%   stresses    - Struct from computeStresses containing boom bending
%                 normal stresses.
%   shearStress - Struct from computePanelShearStress containing panel
%                 shear stresses.
%
% Outputs:
%   vonMises - Struct containing panel von Mises stresses for
%              every saved load case and every physical panel.
%
% Scope:
%   Panel von Mises stress recovery for the current section model.
%
% Method:
%   sigma_vm_psi = sqrt(sigma_panel_psi.^2 + 3*tau_psi.^2)
%
%   sigma_panel_psi is approximated from the nearby/adjacent boom bending
%   stresses at the panel endpoints. For each panel and load case, the
%   endpoint boom stress with the largest absolute value is selected.
%
% Excluded:
%   yield margins, buckling, divergence, failure margins, weight, and
%   optimisation.

validateInputs(stresses, shearStress);

entries = stresses.bending.entries;
panels = shearStress.panels;
nCases = stresses.nCases;
nPanels = shearStress.nPanels;

sigmaBoom_psi = stresses.bending.sigma_psi;
tau_psi = shearStress.panelShear.tau_psi;

panelMatches = repmat(makePanelMatch(), nPanels, 1);
selectedBooms = repmat(makeSelectedBoom(), nCases, nPanels);
sigmaPanel_psi = zeros(nCases, nPanels);

for iPanel = 1:nPanels
    panel = panels(iPanel);
    [startIdx, startDistance] = findNearestBoom(panel.xStart_in, panel.zStart_in, entries);
    [endIdx, endDistance] = findNearestBoom(panel.xEnd_in, panel.zEnd_in, entries);

    panelMatches(iPanel) = makePanelMatch(panel, iPanel, ...
        startIdx, startDistance, entries(startIdx), ...
        endIdx, endDistance, entries(endIdx));

    if ~isfinite(startDistance) || ~isfinite(endDistance)
        error('computeVonMises:InvalidMatchDistance', ...
            'Panel %s has a non-finite boom matching distance.', panel.name);
    end

    candidateIndices = [startIdx, endIdx];
    for iCase = 1:nCases
        candidateSigma = sigmaBoom_psi(iCase, candidateIndices);
        [~, localIdx] = max(abs(candidateSigma));
        selectedIdx = candidateIndices(localIdx);

        sigmaPanel_psi(iCase, iPanel) = sigmaBoom_psi(iCase, selectedIdx);
        selectedBooms(iCase, iPanel) = makeSelectedBoom( ...
            stresses.caseNames(iCase), iPanel, panel, selectedIdx, ...
            entries(selectedIdx), sigmaPanel_psi(iCase, iPanel));
    end
end

sigmaVm_psi = sqrt(sigmaPanel_psi.^2 + 3.0 * tau_psi.^2);

if any(~isfinite(sigmaPanel_psi), 'all') || any(~isfinite(tau_psi), 'all') || ...
        any(~isfinite(sigmaVm_psi), 'all')
    error('computeVonMises:InvalidStressMatrix', ...
        'sigma_panel_psi, tau_psi, and sigma_vm_psi must all contain finite values.');
end
if any(sigmaVm_psi < 0, 'all')
    error('computeVonMises:NegativeVonMises', ...
        'Von Mises stress values must be nonnegative.');
end

vonMises = struct();
vonMises.caseNames = stresses.caseNames;
vonMises.nCases = nCases;
vonMises.nPanels = nPanels;
vonMises.panels = panels;
vonMises.panelMatches = panelMatches;
vonMises.selectedBooms = selectedBooms;
vonMises.panel.sigma_panel_psi = sigmaPanel_psi;
vonMises.panel.tau_psi = tau_psi;
vonMises.panel.sigma_vm_psi = sigmaVm_psi;
vonMises.panel.maxByCase = repmat(makeSummaryEntry(), nCases, 1);

for iCase = 1:nCases
    [maxVm, idxPanel] = max(sigmaVm_psi(iCase, :));
    vonMises.panel.maxByCase(iCase) = makeSummaryEntry( ...
        stresses.caseNames(iCase), maxVm, idxPanel, panels(idxPanel), ...
        sigmaPanel_psi(iCase, idxPanel), tau_psi(iCase, idxPanel), ...
        selectedBooms(iCase, idxPanel));
end

[globalMaxVm, linearIdx] = max(sigmaVm_psi(:));
[globalCaseIdx, globalPanelIdx] = ind2sub(size(sigmaVm_psi), linearIdx);
vonMises.panel.globalMax = makeSummaryEntry( ...
    stresses.caseNames(globalCaseIdx), globalMaxVm, globalPanelIdx, panels(globalPanelIdx), ...
    sigmaPanel_psi(globalCaseIdx, globalPanelIdx), ...
    tau_psi(globalCaseIdx, globalPanelIdx), ...
    selectedBooms(globalCaseIdx, globalPanelIdx));
vonMises.panel.globalMax.case_index = globalCaseIdx;

vonMises.units.normalStress = "psi";
vonMises.units.shearStress = "psi";
vonMises.units.vonMisesStress = "psi";
vonMises.method = ...
    "sigma_vm_psi = sqrt(sigma_panel_psi.^2 + 3*tau_psi.^2); sigma_panel selected from max-absolute adjacent/nearest endpoint boom stress.";
vonMises.assumptions = [
    "Panel von Mises stress recovery for the current section model."
    "Bending normal stress is approximated from adjacent/nearest boom stresses at panel endpoints."
    "For each panel and case, the endpoint boom stress with largest absolute value is selected."
    "Optional panel idealisation into boom areas is handled upstream in computeSectionUnitLoads when enabled."
    "Detailed plate stress recovery within each panel is not implemented."
    "Yield ratios are computed downstream by computeYieldScreening."
    ];
vonMises.warning = ...
    "Panel sigma values are approximated from representative boom stresses; review panel stress recovery before using results for final certification margins.";
vonMises.status = "panel_von_mises_screening";
vonMises.downstreamChecks = [
    "detailed plate stress recovery"
    "buckling"
    "divergence"
    "failure margins"
    "weight"
    "optimisation"
    ];
end

function validateInputs(stresses, shearStress)
if nargin < 1 || isempty(stresses) || ~isstruct(stresses)
    error('computeVonMises:MissingStresses', ...
        'stresses struct from computeStresses is required.');
end
if nargin < 2 || isempty(shearStress) || ~isstruct(shearStress)
    error('computeVonMises:MissingShearStress', ...
        'shearStress struct from computePanelShearStress is required.');
end

if ~isfield(stresses, 'nCases') || ~isfield(stresses, 'caseNames') || ...
        ~isfield(stresses, 'bending')
    error('computeVonMises:InvalidStresses', ...
        'stresses must include nCases, caseNames, and bending fields.');
end
if ~isfield(stresses.bending, 'entries') || ~isfield(stresses.bending, 'sigma_psi')
    error('computeVonMises:MissingBendingData', ...
        'stresses.bending.entries and stresses.bending.sigma_psi are required.');
end

if ~isfield(shearStress, 'nCases') || ~isfield(shearStress, 'nPanels') || ...
        ~isfield(shearStress, 'caseNames') || ~isfield(shearStress, 'panels') || ...
        ~isfield(shearStress, 'panelShear')
    error('computeVonMises:InvalidShearStress', ...
        'shearStress must include nCases, nPanels, caseNames, panels, and panelShear fields.');
end
if ~isfield(shearStress.panelShear, 'tau_psi')
    error('computeVonMises:MissingTauMatrix', ...
        'shearStress.panelShear.tau_psi is required.');
end

if stresses.nCases ~= shearStress.nCases
    error('computeVonMises:CaseCountMismatch', ...
        'stresses and shearStress must have the same number of load cases.');
end
if numel(stresses.caseNames) ~= stresses.nCases || numel(shearStress.caseNames) ~= shearStress.nCases
    error('computeVonMises:CaseNameSizeMismatch', ...
        'Case name arrays must have one name per load case.');
end
if any(string(stresses.caseNames(:)) ~= string(shearStress.caseNames(:)))
    error('computeVonMises:CaseNameMismatch', ...
        'Case names must match between stresses and shearStress.');
end

entries = stresses.bending.entries;
panels = shearStress.panels;
sigmaBoom_psi = stresses.bending.sigma_psi;
tau_psi = shearStress.panelShear.tau_psi;

if isempty(entries)
    error('computeVonMises:EmptyBendingEntries', ...
        'stresses.bending.entries must not be empty.');
end
if isempty(panels) || numel(panels) ~= shearStress.nPanels
    error('computeVonMises:PanelCountMismatch', ...
        'shearStress.panels must contain one panel per tau_psi column.');
end
if ~isequal(size(sigmaBoom_psi), [stresses.nCases, numel(entries)])
    error('computeVonMises:BendingMatrixSizeMismatch', ...
        'stresses.bending.sigma_psi must have size nCases x nBendingEntries.');
end
if ~isequal(size(tau_psi), [shearStress.nCases, shearStress.nPanels])
    error('computeVonMises:TauMatrixSizeMismatch', ...
        'shearStress.panelShear.tau_psi must have size nCases x nPanels.');
end
if any(~isfinite(sigmaBoom_psi), 'all') || any(~isfinite(tau_psi), 'all')
    error('computeVonMises:NonfiniteInputStress', ...
        'Input bending and shear stress matrices must contain finite values.');
end

requiredEntryFields = ["name", "group", "x_in", "z_in"];
entryFieldNames = string(fieldnames(entries));
missingEntryFields = setdiff(requiredEntryFields, entryFieldNames);
if ~isempty(missingEntryFields)
    error('computeVonMises:MissingEntryFields', ...
        'Bending entries are missing required fields: %s', ...
        strjoin(missingEntryFields, ', '));
end

requiredPanelFields = ["name", "group", "xStart_in", "zStart_in", ...
    "xEnd_in", "zEnd_in", "thickness_in"];
panelFieldNames = string(fieldnames(panels));
missingPanelFields = setdiff(requiredPanelFields, panelFieldNames);
if ~isempty(missingPanelFields)
    error('computeVonMises:MissingPanelFields', ...
        'Panel metadata is missing required fields: %s', ...
        strjoin(missingPanelFields, ', '));
end

if any(~isfinite([entries.x_in])) || any(~isfinite([entries.z_in]))
    error('computeVonMises:InvalidBendingEntryCoordinates', ...
        'All bending entry coordinates must be finite.');
end
if any(~isfinite([panels.xStart_in])) || any(~isfinite([panels.zStart_in])) || ...
        any(~isfinite([panels.xEnd_in])) || any(~isfinite([panels.zEnd_in]))
    error('computeVonMises:InvalidPanelCoordinates', ...
        'All panel endpoint coordinates must be finite.');
end
if any(~isfinite([panels.thickness_in])) || any([panels.thickness_in] <= 0)
    error('computeVonMises:InvalidPanelThickness', ...
        'All panel thicknesses must be positive and finite.');
end

if isfield(stresses.bending, 'units') && isfield(stresses.bending.units, 'stress') && ...
        string(stresses.bending.units.stress) ~= "psi"
    error('computeVonMises:UnexpectedBendingStressUnits', ...
        'Expected bending stress units to be psi.');
end
if isfield(shearStress, 'units') && isfield(shearStress.units, 'stress') && ...
        string(shearStress.units.stress) ~= "psi"
    error('computeVonMises:UnexpectedShearStressUnits', ...
        'Expected shear stress units to be psi.');
end
end

function [idx, distance] = findNearestBoom(x_in, z_in, entries)
dx = [entries.x_in] - x_in;
dz = [entries.z_in] - z_in;
distances = hypot(dx, dz);
[distance, idx] = min(distances);
if isempty(idx) || ~isfinite(distance)
    error('computeVonMises:NoUsableBoomMatch', ...
        'Could not find a usable boom match for panel endpoint at x=%.6g in, z=%.6g in.', ...
        x_in, z_in);
end
end

function match = makePanelMatch(panel, panelIndex, startIdx, startDistance, startEntry, endIdx, endDistance, endEntry)
if nargin == 0
    panel = makePanelPlaceholder();
    panelIndex = NaN;
    startIdx = NaN;
    startDistance = NaN;
    startEntry = makeEntryPlaceholder();
    endIdx = NaN;
    endDistance = NaN;
    endEntry = makeEntryPlaceholder();
end

match = struct();
match.panel_index = panelIndex;
match.panel_name = panel.name;
match.panel_group = panel.group;
match.xStart_in = panel.xStart_in;
match.zStart_in = panel.zStart_in;
match.xEnd_in = panel.xEnd_in;
match.zEnd_in = panel.zEnd_in;
match.start_boom_index = startIdx;
match.start_boom_name = startEntry.name;
match.start_boom_group = startEntry.group;
match.start_boom_x_in = startEntry.x_in;
match.start_boom_z_in = startEntry.z_in;
match.start_match_distance_in = startDistance;
match.end_boom_index = endIdx;
match.end_boom_name = endEntry.name;
match.end_boom_group = endEntry.group;
match.end_boom_x_in = endEntry.x_in;
match.end_boom_z_in = endEntry.z_in;
match.end_match_distance_in = endDistance;
match.method = "nearest endpoint boom match; exact endpoint matches have zero distance";
end

function selected = makeSelectedBoom(caseName, panelIndex, panel, boomIndex, boomEntry, sigmaPanel_psi)
if nargin == 0
    caseName = "";
    panelIndex = NaN;
    panel = makePanelPlaceholder();
    boomIndex = NaN;
    boomEntry = makeEntryPlaceholder();
    sigmaPanel_psi = NaN;
end

selected = struct();
selected.case_name = string(caseName);
selected.panel_index = panelIndex;
selected.panel_name = panel.name;
selected.panel_group = panel.group;
selected.boom_index = boomIndex;
selected.boom_name = boomEntry.name;
selected.boom_group = boomEntry.group;
selected.boom_x_in = boomEntry.x_in;
selected.boom_z_in = boomEntry.z_in;
selected.sigma_panel_psi = sigmaPanel_psi;
selected.selection_rule = "largest absolute bending stress among nearest start/end boom matches";
end

function entry = makeSummaryEntry(caseName, sigmaVm_psi, panelIndex, panel, sigmaPanel_psi, tau_psi, selectedBoom)
if nargin == 0
    caseName = "";
    sigmaVm_psi = NaN;
    panelIndex = NaN;
    panel = makePanelPlaceholder();
    sigmaPanel_psi = NaN;
    tau_psi = NaN;
    selectedBoom = makeSelectedBoom();
end

entry = struct();
entry.case_name = string(caseName);
entry.sigma_vm_psi = sigmaVm_psi;
entry.panel_index = panelIndex;
entry.panel_group = panel.group;
entry.panel_name = panel.name;
entry.xStart_in = panel.xStart_in;
entry.zStart_in = panel.zStart_in;
entry.xEnd_in = panel.xEnd_in;
entry.zEnd_in = panel.zEnd_in;
entry.thickness_in = panel.thickness_in;
entry.sigma_panel_psi = sigmaPanel_psi;
entry.tau_psi = tau_psi;
entry.selected_boom_index = selectedBoom.boom_index;
entry.selected_boom_group = selectedBoom.boom_group;
entry.selected_boom_name = selectedBoom.boom_name;
entry.selected_boom_x_in = selectedBoom.boom_x_in;
entry.selected_boom_z_in = selectedBoom.boom_z_in;
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

function entry = makeEntryPlaceholder()
entry = struct();
entry.name = "";
entry.group = "";
entry.x_in = NaN;
entry.z_in = NaN;
end
