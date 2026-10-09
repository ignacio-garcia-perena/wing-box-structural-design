function panelBuckling = computePanelBuckling(wing, shearStress, vonMises, analysisOptions)
% computePanelBuckling Station-level flat-panel buckling diagnostic.
%
% Scope:
%   Separate panel compression, bending, and shear buckling screens using
%   project-approved constant high-aspect-ratio coefficients. When enabled,
%   Combined load ratios are also computed. Crippling, final
%   margins, and coefficient interpolation are intentionally excluded.

if nargin < 4 || isempty(analysisOptions)
    analysisOptions = defaultAnalysisOptions();
end

validateInputs(wing, shearStress, vonMises, analysisOptions);

Kc = analysisOptions.panelBuckling.Kc;
Kb = analysisOptions.panelBuckling.Kb;
Ks = analysisOptions.panelBuckling.Ks;

material = wing.material.skinWeb;
E_psi = material.E_psi;
if isfield(material, 'Ec_psi') && isfinite(material.Ec_psi) && material.Ec_psi > 0
    Ec_psi = material.Ec_psi;
else
    Ec_psi = material.E_psi;
end
nu = material.nu;

panels = shearStress.panels;
nCases = shearStress.nCases;
nPanels = shearStress.nPanels;
a_in = wing.ribs.spacing_in;

b_in = zeros(1, nPanels);
t_in = [panels.thickness_in];
for iPanel = 1:nPanels
    b_in(iPanel) = panelWidth(panels(iPanel));
end

denominator = 12.0 * (1.0 - nu^2);
sigma_cr_c_psi = Kc * pi^2 * Ec_psi ./ denominator .* (t_in ./ b_in).^2;
sigma_cr_b_psi = Kb * pi^2 * Ec_psi ./ denominator .* (t_in ./ b_in).^2;
tau_cr_s_psi = Ks * pi^2 * E_psi ./ denominator .* (t_in ./ b_in).^2;

sigma_panel_psi = vonMises.panel.sigma_panel_psi;
tau_psi = shearStress.panelShear.tau_psi;
sigma_compression_psi = max(0.0, -sigma_panel_psi);
tau_applied_psi = abs(tau_psi);

Rc = sigma_compression_psi ./ sigma_cr_c_psi;
Rb = sigma_compression_psi ./ sigma_cr_b_psi;
Rs = tau_applied_psi ./ tau_cr_s_psi;
R_SC = Rs.^1.5 + Rc;
R_BS = Rb.^2 + Rs.^2;
[R_combined_active, combinedMode, combinedWarnings] = assignCombinedRatios( ...
    R_SC, R_BS, panels, shearStress.caseNames, analysisOptions);

if any(~isfinite(Rc), 'all') || any(~isfinite(Rb), 'all') || ...
        any(~isfinite(Rs), 'all') || any(~isfinite(R_SC), 'all') || ...
        any(~isfinite(R_BS), 'all') || any(~isfinite(R_combined_active), 'all')
    error('computePanelBuckling:InvalidRatios', ...
        'Panel buckling ratios must be finite.');
end

panelBuckling = struct();
panelBuckling.caseNames = shearStress.caseNames;
panelBuckling.nCases = nCases;
panelBuckling.nPanels = nPanels;
panelBuckling.panels = panels;
panelBuckling.enabled = true;
panelBuckling.coefficients.Kc = Kc;
panelBuckling.coefficients.Kb = Kb;
panelBuckling.coefficients.Ks = Ks;
panelBuckling.coefficients.method = analysisOptions.panelBuckling.method;
panelBuckling.coefficients.assumption = analysisOptions.panelBuckling.assumption;
panelBuckling.material.name = material.name;
panelBuckling.material.temper = material.temper;
panelBuckling.material.product_form = material.product_form;
panelBuckling.material.E_psi = E_psi;
panelBuckling.material.Ec_psi = Ec_psi;
panelBuckling.material.nu = nu;

panelBuckling.geometry.a_in = a_in;
panelBuckling.geometry.b_in = b_in;
panelBuckling.geometry.t_in = t_in;
panelBuckling.geometry.aspect_ratio_a_over_b = a_in ./ b_in;

panelBuckling.critical.sigma_cr_c_psi = sigma_cr_c_psi;
panelBuckling.critical.sigma_cr_b_psi = sigma_cr_b_psi;
panelBuckling.critical.tau_cr_s_psi = tau_cr_s_psi;
panelBuckling.applied.sigma_panel_psi = sigma_panel_psi;
panelBuckling.applied.sigma_compression_psi = sigma_compression_psi;
panelBuckling.applied.tau_psi = tau_psi;
panelBuckling.applied.tau_abs_psi = tau_applied_psi;
panelBuckling.ratios.Rc = Rc;
panelBuckling.ratios.Rb = Rb;
panelBuckling.ratios.Rs = Rs;
panelBuckling.R_SC = R_SC;
panelBuckling.R_BS = R_BS;
panelBuckling.R_combined_active = R_combined_active;
panelBuckling.combinedMode = combinedMode;
panelBuckling.combinedWarnings = combinedWarnings;
panelBuckling.combined.enabled = getCombinedEnabled(analysisOptions);
panelBuckling.combined.method = string(analysisOptions.panelBuckling.combinedMethod);
panelBuckling.combined.notes = string(analysisOptions.panelBuckling.combinedNotes);

panelBuckling.maxRc = makeGlobalSummary(Rc, "Rc", shearStress.caseNames, panels, ...
    sigma_compression_psi, sigma_cr_c_psi, b_in, t_in, a_in);
panelBuckling.maxRb = makeGlobalSummary(Rb, "Rb", shearStress.caseNames, panels, ...
    sigma_compression_psi, sigma_cr_b_psi, b_in, t_in, a_in);
panelBuckling.maxRs = makeGlobalSummary(Rs, "Rs", shearStress.caseNames, panels, ...
    tau_applied_psi, tau_cr_s_psi, b_in, t_in, a_in);
panelBuckling.maxR_SC = makeGlobalSummary(R_SC, "R_SC", shearStress.caseNames, panels, ...
    R_SC, ones(1, nPanels), b_in, t_in, a_in);
panelBuckling.maxR_BS = makeGlobalSummary(R_BS, "R_BS", shearStress.caseNames, panels, ...
    R_BS, ones(1, nPanels), b_in, t_in, a_in);
panelBuckling.maxCombinedActive = makeGlobalSummary(R_combined_active, ...
    "R_combined_active", shearStress.caseNames, panels, R_combined_active, ...
    ones(1, nPanels), b_in, t_in, a_in);
panelBuckling.maxCombinedActive.combined_mode = ...
    combinedMode(panelBuckling.maxCombinedActive.case_index, ...
    panelBuckling.maxCombinedActive.panel_index);

panelBuckling.maxByCase = repmat(makeCaseSummary(), nCases, 1);
for iCase = 1:nCases
    panelBuckling.maxByCase(iCase) = makeCaseSummary( ...
        shearStress.caseNames(iCase), iCase, Rc(iCase, :), Rb(iCase, :), ...
        Rs(iCase, :), panels);
end

panelBuckling.pass = panelBuckling.maxRc.ratio <= 1.0 && ...
    panelBuckling.maxRb.ratio <= 1.0 && panelBuckling.maxRs.ratio <= 1.0;
panelBuckling.combinedPass = panelBuckling.maxCombinedActive.ratio <= 1.0;
if panelBuckling.combined.enabled
    panelBuckling.pass = panelBuckling.pass && panelBuckling.combinedPass;
end
panelBuckling.units.stress = "psi";
panelBuckling.units.length = "in";
panelBuckling.method = ...
    "Flat-panel buckling with constant high-aspect-ratio Kc, Kb, Ks coefficients; Rc/Rb/Rs are separate screens only.";
panelBuckling.assumptions = [
    "All physical skin/web panels use wing.material.skinWeb properties."
    "Panel bay length a_in is wing.ribs.spacing_in."
    "Panel width b_in is panel ds_in when available, otherwise endpoint distance."
    "Compression and bending screens use max(0, -sigma_panel_psi) from current panel stress recovery."
    "Shear screen uses abs(tau_psi)."
    "Combined S+C and B+S ratios are computed when enabled."
    "No crippling, postbuckling, or final margins are included."
    ];
end

function validateInputs(wing, shearStress, vonMises, analysisOptions)
if nargin < 1 || isempty(wing) || ~isstruct(wing)
    error('computePanelBuckling:MissingWing', 'wing struct is required.');
end
if ~isfield(wing, 'ribs') || ~isfield(wing.ribs, 'spacing_in') || ...
        ~(isnumeric(wing.ribs.spacing_in) && isscalar(wing.ribs.spacing_in) && ...
        isfinite(wing.ribs.spacing_in) && wing.ribs.spacing_in > 0)
    error('computePanelBuckling:InvalidRibSpacing', ...
        'wing.ribs.spacing_in must be a positive finite scalar.');
end
if ~isfield(wing, 'material') || ~isfield(wing.material, 'skinWeb')
    error('computePanelBuckling:MissingSkinWebMaterial', ...
        'wing.material.skinWeb is required.');
end
material = wing.material.skinWeb;
requiredMaterialFields = ["E_psi", "nu"];
missingMaterialFields = setdiff(requiredMaterialFields, string(fieldnames(material)));
if ~isempty(missingMaterialFields)
    error('computePanelBuckling:MissingMaterialFields', ...
        'skinWeb material is missing fields: %s', strjoin(missingMaterialFields, ', '));
end
if ~(isnumeric(material.E_psi) && isscalar(material.E_psi) && ...
        isfinite(material.E_psi) && material.E_psi > 0)
    error('computePanelBuckling:InvalidE', 'skinWeb E_psi must be positive and finite.');
end
if ~(isnumeric(material.nu) && isscalar(material.nu) && ...
        isfinite(material.nu) && material.nu > -1 && material.nu < 0.5)
    error('computePanelBuckling:InvalidNu', 'skinWeb nu must be finite and between -1 and 0.5.');
end

if nargin < 2 || isempty(shearStress) || ~isstruct(shearStress)
    error('computePanelBuckling:MissingShearStress', ...
        'shearStress struct from computePanelShearStress is required.');
end
if nargin < 3 || isempty(vonMises) || ~isstruct(vonMises)
    error('computePanelBuckling:MissingVonMises', ...
        'vonMises struct from computeVonMises is required.');
end
if ~isfield(shearStress, 'nCases') || ~isfield(shearStress, 'nPanels') || ...
        ~isfield(shearStress, 'caseNames') || ~isfield(shearStress, 'panels') || ...
        ~isfield(shearStress, 'panelShear') || ~isfield(shearStress.panelShear, 'tau_psi')
    error('computePanelBuckling:InvalidShearStress', ...
        'shearStress must include nCases, nPanels, caseNames, panels, and panelShear.tau_psi.');
end
if ~isfield(vonMises, 'panel') || ~isfield(vonMises.panel, 'sigma_panel_psi')
    error('computePanelBuckling:MissingSigmaPanel', ...
        'vonMises.panel.sigma_panel_psi is required.');
end
if ~isequal(size(shearStress.panelShear.tau_psi), [shearStress.nCases, shearStress.nPanels])
    error('computePanelBuckling:TauSizeMismatch', ...
        'tau_psi must be nCases x nPanels.');
end
if ~isequal(size(vonMises.panel.sigma_panel_psi), [shearStress.nCases, shearStress.nPanels])
    error('computePanelBuckling:SigmaSizeMismatch', ...
        'sigma_panel_psi must match shearStress nCases x nPanels.');
end
if any(string(shearStress.caseNames(:)) ~= string(vonMises.caseNames(:)))
    error('computePanelBuckling:CaseNameMismatch', ...
        'Case names must match between shearStress and vonMises.');
end
if nargin < 4 || isempty(analysisOptions) || ~isstruct(analysisOptions) || ...
        ~isfield(analysisOptions, 'panelBuckling')
    error('computePanelBuckling:MissingPanelBucklingOptions', ...
        'analysisOptions.panelBuckling is required.');
end
if ~isfield(analysisOptions.panelBuckling, 'enableCombinedRatios')
    analysisOptions.panelBuckling.enableCombinedRatios = true;
end
if ~isfield(analysisOptions.panelBuckling, 'combinedMethod')
    analysisOptions.panelBuckling.combinedMethod = "combined_loads_interaction";
end
if ~isfield(analysisOptions.panelBuckling, 'combinedNotes')
    analysisOptions.panelBuckling.combinedNotes = ...
        "Combined ratios: S+C uses Rs^1.5 + Rc < 1; B+S uses Rb^2 + Rs^2 < 1.";
end
requiredCoeffFields = ["Kc", "Kb", "Ks"];
missingCoeffFields = setdiff(requiredCoeffFields, string(fieldnames(analysisOptions.panelBuckling)));
if ~isempty(missingCoeffFields)
    error('computePanelBuckling:MissingCoefficientFields', ...
        'panelBuckling options are missing fields: %s', strjoin(missingCoeffFields, ', '));
end
coeffs = [analysisOptions.panelBuckling.Kc, analysisOptions.panelBuckling.Kb, ...
    analysisOptions.panelBuckling.Ks];
if any(~isfinite(coeffs)) || any(coeffs <= 0)
    error('computePanelBuckling:InvalidCoefficients', ...
        'Kc, Kb, and Ks must be positive finite scalars.');
end

requiredPanelFields = ["name", "group", "xStart_in", "zStart_in", ...
    "xEnd_in", "zEnd_in", "thickness_in"];
missingPanelFields = setdiff(requiredPanelFields, string(fieldnames(shearStress.panels)));
if ~isempty(missingPanelFields)
    error('computePanelBuckling:MissingPanelFields', ...
        'Panel metadata is missing fields: %s', strjoin(missingPanelFields, ', '));
end
if any(~isfinite([shearStress.panels.thickness_in])) || any([shearStress.panels.thickness_in] <= 0)
    error('computePanelBuckling:InvalidPanelThickness', ...
        'Panel thickness values must be positive and finite.');
end
end

function b_in = panelWidth(panel)
if isfield(panel, 'ds_in') && isfinite(panel.ds_in) && panel.ds_in > 0
    b_in = panel.ds_in;
    return;
end
dx = panel.xEnd_in - panel.xStart_in;
dz = panel.zEnd_in - panel.zStart_in;
b_in = hypot(dx, dz);
if ~isfinite(b_in) || b_in <= 0
    error('computePanelBuckling:InvalidPanelWidth', ...
        'Panel %s has invalid width.', panel.name);
end
end

function entry = makeGlobalSummary(ratioMatrix, ratioName, caseNames, panels, ...
    appliedMatrix, criticalByPanel, b_in, t_in, a_in)
[ratio, linearIdx] = max(ratioMatrix(:));
[caseIdx, panelIdx] = ind2sub(size(ratioMatrix), linearIdx);
panel = panels(panelIdx);

entry = struct();
entry.ratio_name = string(ratioName);
entry.ratio = ratio;
entry.case_index = caseIdx;
entry.case_name = string(caseNames(caseIdx));
entry.panel_index = panelIdx;
entry.panel_name = panel.name;
entry.panel_group = panel.group;
entry.applied_stress_psi = appliedMatrix(caseIdx, panelIdx);
entry.critical_stress_psi = criticalByPanel(panelIdx);
entry.a_in = a_in;
entry.b_in = b_in(panelIdx);
entry.t_in = t_in(panelIdx);
entry.aspect_ratio_a_over_b = a_in / b_in(panelIdx);
end

function [RActive, modes, warnings] = assignCombinedRatios(R_SC, R_BS, panels, caseNames, analysisOptions)
RActive = zeros(size(R_SC));
modes = strings(size(R_SC));
warnings = strings(0, 1);

combinedEnabled = getCombinedEnabled(analysisOptions);
if ~combinedEnabled
    RActive(:) = 0.0;
    modes(:) = "disabled";
    return;
end

for iCase = 1:numel(caseNames)
    caseName = upper(string(caseNames(iCase)));
    for iPanel = 1:numel(panels)
        panel = panels(iPanel);
        [mode, warningText] = classifyCombinedMode(panel, caseName);
        if strlength(warningText) > 0
            warnings(end + 1, 1) = warningText; %#ok<AGROW>
        end

        switch mode
            case "S+C"
                RActive(iCase, iPanel) = R_SC(iCase, iPanel);
            case "B+S"
                RActive(iCase, iPanel) = R_BS(iCase, iPanel);
            otherwise
                RActive(iCase, iPanel) = max(R_SC(iCase, iPanel), R_BS(iCase, iPanel));
        end
        modes(iCase, iPanel) = mode;
    end
end
warnings = unique(warnings);
end

function enabled = getCombinedEnabled(analysisOptions)
enabled = isfield(analysisOptions, 'panelBuckling') && ...
    isfield(analysisOptions.panelBuckling, 'enableCombinedRatios') && ...
    logical(analysisOptions.panelBuckling.enableCombinedRatios);
end

function [mode, warningText] = classifyCombinedMode(panel, caseName)
panelName = lower(string(panel.name));
panelGroup = lower(string(panel.group));
warningText = "";

isUpperSkin = contains(panelGroup, "top") || contains(panelName, "top");
isLowerSkin = contains(panelGroup, "bottom") || contains(panelName, "bottom");
isNose = contains(panelGroup, "nose") || contains(panelName, "nose") || ...
    contains(panelName, "leadingedge") || contains(panelGroup, "leadingedge");
isSpar = contains(panelGroup, "spar") || contains(panelName, "spar");
isTrailing = contains(panelGroup, "trailing") || contains(panelName, "trailing");

if isUpperSkin && (caseName == "PHAA" || caseName == "PLAA")
    mode = "S+C";
elseif isLowerSkin && (startsWith(caseName, "NHAA") || startsWith(caseName, "NLAA"))
    mode = "S+C";
elseif isNose || isSpar || isTrailing
    mode = "B+S";
else
    mode = "max_due_to_ambiguous_mapping";
    warningText = "Ambiguous combined panel mapping for panel " + ...
        string(panel.name) + " group " + string(panel.group) + ...
        " case " + caseName + "; used max(R_SC, R_BS).";
end
end

function entry = makeCaseSummary(caseName, caseIdx, RcRow, RbRow, RsRow, panels)
entry = struct();
if nargin == 0
    entry.case_name = "";
    entry.case_index = NaN;
    entry.max_Rc = NaN;
    entry.max_Rc_panel = "";
    entry.max_Rb = NaN;
    entry.max_Rb_panel = "";
    entry.max_Rs = NaN;
    entry.max_Rs_panel = "";
    entry.pass = false;
    return;
end

[maxRc, idxRc] = max(RcRow);
[maxRb, idxRb] = max(RbRow);
[maxRs, idxRs] = max(RsRow);
entry.case_name = string(caseName);
entry.case_index = caseIdx;
entry.max_Rc = maxRc;
entry.max_Rc_panel = panels(idxRc).name;
entry.max_Rb = maxRb;
entry.max_Rb_panel = panels(idxRb).name;
entry.max_Rs = maxRs;
entry.max_Rs_panel = panels(idxRs).name;
entry.pass = maxRc <= 1.0 && maxRb <= 1.0 && maxRs <= 1.0;
end
