function panelBucklingBay = computePanelBucklingBayAveraged(wing, section, stationResultInboard, stationResultOutboard, analysisOptions)
% computePanelBucklingBayAveraged Bay-level panel buckling recovery.
%
% Uses bay-averaged recovery for panel buckling:
%   - average compressive stress from the four surrounding boom/stringer
%     stresses at the inboard and outboard rib stations,
%   - average total shear flow between the same two rib stations.
%
% Existing station endpoint-envelope panel buckling remains useful as a
% diagnostic, but this function is intended for the final analysis
% panel buckling pass/fail when enabled.

if nargin < 5 || isempty(analysisOptions)
    analysisOptions = defaultAnalysisOptions();
end

validateInputs(wing, section, stationResultInboard, stationResultOutboard, analysisOptions);

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

panelsIn = stationResultInboard.shearFlows.panels;
panelsOut = stationResultOutboard.shearFlows.panels;
validatePanelConsistency(panelsIn, panelsOut);

caseNames = string(stationResultInboard.stresses.caseNames(:));
nCases = numel(caseNames);
nPanels = numel(panelsIn);
bayStart_in = stationResultInboard.y_in;
bayEnd_in = stationResultOutboard.y_in;
bayLength_in = bayEnd_in - bayStart_in;

entriesIn = stationResultInboard.stresses.bending.entries;
entriesOut = stationResultOutboard.stresses.bending.entries;
sigmaBoomIn = stationResultInboard.stresses.bending.sigma_psi;
sigmaBoomOut = stationResultOutboard.stresses.bending.sigma_psi;
qIn = stationResultInboard.shearFlows.shear.q_lb_per_in;
qOut = stationResultOutboard.shearFlows.shear.q_lb_per_in;

b_in = zeros(1, nPanels);
t_in = [panelsIn.thickness_in];
sigmaCompressionAvg_psi = zeros(nCases, nPanels);
qAvg_lb_per_in = zeros(nCases, nPanels);
tauAvg_psi = zeros(nCases, nPanels);
panelMatches = repmat(makePanelMatch(), nPanels, 1);

for iPanel = 1:nPanels
    panel = panelsIn(iPanel);
    b_in(iPanel) = panelWidth(panel);
    [startIdxIn, startDistIn] = findNearestBoom(panel.xStart_in, panel.zStart_in, entriesIn);
    [endIdxIn, endDistIn] = findNearestBoom(panel.xEnd_in, panel.zEnd_in, entriesIn);
    [startIdxOut, startDistOut] = findNearestBoom(panel.xStart_in, panel.zStart_in, entriesOut);
    [endIdxOut, endDistOut] = findNearestBoom(panel.xEnd_in, panel.zEnd_in, entriesOut);

    panelMatches(iPanel) = makePanelMatch(panel, iPanel, ...
        startIdxIn, endIdxIn, startIdxOut, endIdxOut, ...
        startDistIn, endDistIn, startDistOut, endDistOut, ...
        entriesIn(startIdxIn), entriesIn(endIdxIn), ...
        entriesOut(startIdxOut), entriesOut(endIdxOut));

    for iCase = 1:nCases
        sigmaFour = [
            sigmaBoomIn(iCase, startIdxIn)
            sigmaBoomIn(iCase, endIdxIn)
            sigmaBoomOut(iCase, startIdxOut)
            sigmaBoomOut(iCase, endIdxOut)
            ];
        compressionFour = max(0.0, -sigmaFour);
        sigmaCompressionAvg_psi(iCase, iPanel) = mean(compressionFour);
        qAvg_lb_per_in(iCase, iPanel) = 0.5 * (qIn(iCase, iPanel) + qOut(iCase, iPanel));
        tauAvg_psi(iCase, iPanel) = abs(qAvg_lb_per_in(iCase, iPanel)) / t_in(iPanel);
    end
end

denominator = 12.0 * (1.0 - nu^2);
sigma_cr_c_psi = Kc * pi^2 * Ec_psi ./ denominator .* (t_in ./ b_in).^2;
sigma_cr_b_psi = Kb * pi^2 * Ec_psi ./ denominator .* (t_in ./ b_in).^2;
tau_cr_s_psi = Ks * pi^2 * E_psi ./ denominator .* (t_in ./ b_in).^2;
[sigma_cr_c_effective_psi, curvedCorrection] = applyCurvedPanelCorrection( ...
    panelsIn, sigma_cr_c_psi, Ec_psi, bayLength_in, wing, analysisOptions);

Rc = sigmaCompressionAvg_psi ./ sigma_cr_c_effective_psi;
Rb = sigmaCompressionAvg_psi ./ sigma_cr_b_psi;
Rs = tauAvg_psi ./ tau_cr_s_psi;
R_SC = Rs.^1.5 + Rc;
R_BS = Rb.^2 + Rs.^2;
[R_combined_active, combinedMode, combinedWarnings] = assignCombinedRatios( ...
    R_SC, R_BS, panelsIn, caseNames, analysisOptions);

if any(~isfinite(Rc), 'all') || any(~isfinite(Rb), 'all') || ...
        any(~isfinite(Rs), 'all') || any(~isfinite(R_SC), 'all') || ...
        any(~isfinite(R_BS), 'all') || any(~isfinite(R_combined_active), 'all')
    error('computePanelBucklingBayAveraged:InvalidRatios', ...
        'Bay-averaged panel buckling ratios must be finite.');
end

panelBucklingBay = struct();
panelBucklingBay.enabled = true;
panelBucklingBay.recoveryMethod = string(analysisOptions.panelBuckling.recoveryMethod);
panelBucklingBay.method = ...
    "Bay-averaged panel buckling recovery with constant Kc/Kb/Ks coefficients.";
panelBucklingBay.caseNames = caseNames;
panelBucklingBay.nCases = nCases;
panelBucklingBay.nPanels = nPanels;
panelBucklingBay.panels = panelsIn;
panelBucklingBay.panelMatches = panelMatches;
panelBucklingBay.bay.start_y_in = bayStart_in;
panelBucklingBay.bay.end_y_in = bayEnd_in;
panelBucklingBay.bay.mid_y_in = 0.5 * (bayStart_in + bayEnd_in);
panelBucklingBay.bay.length_in = bayLength_in;
panelBucklingBay.coefficients.Kc = Kc;
panelBucklingBay.coefficients.Kb = Kb;
panelBucklingBay.coefficients.Ks = Ks;
panelBucklingBay.coefficients.method = analysisOptions.panelBuckling.method;
panelBucklingBay.coefficients.assumption = analysisOptions.panelBuckling.assumption;
panelBucklingBay.material.name = material.name;
panelBucklingBay.material.temper = material.temper;
panelBucklingBay.material.product_form = material.product_form;
panelBucklingBay.material.E_psi = E_psi;
panelBucklingBay.material.Ec_psi = Ec_psi;
panelBucklingBay.material.nu = nu;
panelBucklingBay.geometry.a_in = bayLength_in;
panelBucklingBay.geometry.b_in = b_in;
panelBucklingBay.geometry.t_in = t_in;
panelBucklingBay.geometry.aspect_ratio_a_over_b = bayLength_in ./ b_in;
panelBucklingBay.critical.sigma_cr_c_psi = sigma_cr_c_psi;
panelBucklingBay.critical.sigma_cr_c_flat_psi = sigma_cr_c_psi;
panelBucklingBay.critical.sigma_cr_c_effective_psi = sigma_cr_c_effective_psi;
panelBucklingBay.critical.sigma_cr_cylinder_psi = curvedCorrection.sigma_cr_cylinder_psi;
panelBucklingBay.critical.sigma_cr_b_psi = sigma_cr_b_psi;
panelBucklingBay.critical.tau_cr_s_psi = tau_cr_s_psi;
panelBucklingBay.curvedCorrection = curvedCorrection;
panelBucklingBay.applied.sigma_compression_avg_psi = sigmaCompressionAvg_psi;
panelBucklingBay.applied.q_avg_lb_per_in = qAvg_lb_per_in;
panelBucklingBay.applied.tau_avg_psi = tauAvg_psi;
panelBucklingBay.ratios.Rc = Rc;
panelBucklingBay.ratios.Rb = Rb;
panelBucklingBay.ratios.Rs = Rs;
panelBucklingBay.R_SC = R_SC;
panelBucklingBay.R_BS = R_BS;
panelBucklingBay.R_combined_active = R_combined_active;
panelBucklingBay.combinedMode = combinedMode;
panelBucklingBay.combinedWarnings = combinedWarnings;
panelBucklingBay.combined.enabled = getCombinedEnabled(analysisOptions);
panelBucklingBay.combined.method = string(analysisOptions.panelBuckling.combinedMethod);
panelBucklingBay.combined.notes = string(analysisOptions.panelBuckling.combinedNotes);
panelBucklingBay.maxRc = makeGlobalSummary(Rc, "Rc", caseNames, panelsIn, ...
    sigmaCompressionAvg_psi, sigma_cr_c_effective_psi, b_in, t_in, bayLength_in);
panelBucklingBay.maxRb = makeGlobalSummary(Rb, "Rb", caseNames, panelsIn, ...
    sigmaCompressionAvg_psi, sigma_cr_b_psi, b_in, t_in, bayLength_in);
panelBucklingBay.maxRs = makeGlobalSummary(Rs, "Rs", caseNames, panelsIn, ...
    tauAvg_psi, tau_cr_s_psi, b_in, t_in, bayLength_in);
panelBucklingBay.maxR_SC = makeGlobalSummary(R_SC, "R_SC", caseNames, panelsIn, ...
    R_SC, ones(1, nPanels), b_in, t_in, bayLength_in);
panelBucklingBay.maxR_BS = makeGlobalSummary(R_BS, "R_BS", caseNames, panelsIn, ...
    R_BS, ones(1, nPanels), b_in, t_in, bayLength_in);
panelBucklingBay.maxCombinedActive = makeGlobalSummary(R_combined_active, ...
    "R_combined_active", caseNames, panelsIn, R_combined_active, ...
    ones(1, nPanels), b_in, t_in, bayLength_in);
panelBucklingBay.maxCombinedActive.combined_mode = ...
    combinedMode(panelBucklingBay.maxCombinedActive.case_index, ...
    panelBucklingBay.maxCombinedActive.panel_index);
panelBucklingBay.pass = panelBucklingBay.maxRc.ratio <= 1.0 && ...
    panelBucklingBay.maxRb.ratio <= 1.0 && panelBucklingBay.maxRs.ratio <= 1.0;
panelBucklingBay.combinedPass = panelBucklingBay.maxCombinedActive.ratio <= 1.0;
if panelBucklingBay.combined.enabled
    panelBucklingBay.pass = panelBucklingBay.pass && panelBucklingBay.combinedPass;
end
panelBucklingBay.units.stress = "psi";
panelBucklingBay.units.length = "in";
panelBucklingBay.units.shearFlow = "lb/in";
panelBucklingBay.assumptions = [
    "Compression uses the mean of max(0,-sigma) from the four panel endpoint boom stresses at the two bay stations."
    "Shear uses abs of the bay-averaged total shear flow divided by panel thickness."
    "Curved-panel correction augments compression critical stress only for curved airfoil skin/nose panels when enabled."
    "No postbuckling, crippling, or final margins are included."
    ];
end

function validateInputs(wing, section, stationIn, stationOut, analysisOptions)
if nargin < 1 || isempty(wing) || ~isstruct(wing) || ...
        ~isfield(wing, 'material') || ~isfield(wing.material, 'skinWeb')
    error('computePanelBucklingBayAveraged:InvalidWing', ...
        'wing.material.skinWeb is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computePanelBucklingBayAveraged:InvalidSection', ...
        'section struct is required.');
end
requiredStationFields = ["y_in", "stresses", "shearFlows"];
missingIn = setdiff(requiredStationFields, string(fieldnames(stationIn)));
missingOut = setdiff(requiredStationFields, string(fieldnames(stationOut)));
if ~isempty(missingIn) || ~isempty(missingOut)
    error('computePanelBucklingBayAveraged:InvalidStationResults', ...
        'Both station results must include y_in, stresses, and shearFlows.');
end
if ~(stationOut.y_in > stationIn.y_in)
    error('computePanelBucklingBayAveraged:InvalidBayLength', ...
        'Outboard station must be farther outboard than inboard station.');
end
if nargin < 5 || isempty(analysisOptions) || ~isstruct(analysisOptions) || ...
        ~isfield(analysisOptions, 'panelBuckling')
    error('computePanelBucklingBayAveraged:MissingOptions', ...
        'analysisOptions.panelBuckling is required.');
end
requiredCoeffFields = ["Kc", "Kb", "Ks"];
missingCoeffFields = setdiff(requiredCoeffFields, string(fieldnames(analysisOptions.panelBuckling)));
if ~isempty(missingCoeffFields)
    error('computePanelBucklingBayAveraged:MissingCoefficients', ...
        'panelBuckling options are missing fields: %s', strjoin(missingCoeffFields, ', '));
end
coeffs = [analysisOptions.panelBuckling.Kc, analysisOptions.panelBuckling.Kb, ...
    analysisOptions.panelBuckling.Ks];
if any(~isfinite(coeffs)) || any(coeffs <= 0)
    error('computePanelBucklingBayAveraged:InvalidCoefficients', ...
        'Kc, Kb, and Ks must be positive finite scalars.');
end
material = wing.material.skinWeb;
requiredMaterialFields = ["E_psi", "nu"];
missingMaterialFields = setdiff(requiredMaterialFields, string(fieldnames(material)));
if ~isempty(missingMaterialFields)
    error('computePanelBucklingBayAveraged:MissingMaterialFields', ...
        'skinWeb material is missing fields: %s', strjoin(missingMaterialFields, ', '));
end
if ~(isfinite(material.E_psi) && material.E_psi > 0 && ...
        isfinite(material.nu) && material.nu > -1 && material.nu < 0.5)
    error('computePanelBucklingBayAveraged:InvalidMaterial', ...
        'skinWeb E_psi and nu must be valid.');
end
end

function validatePanelConsistency(panelsIn, panelsOut)
if numel(panelsIn) ~= numel(panelsOut)
    error('computePanelBucklingBayAveraged:PanelCountMismatch', ...
        'Inboard and outboard stations must have the same panel count.');
end
for iPanel = 1:numel(panelsIn)
    if string(panelsIn(iPanel).name) ~= string(panelsOut(iPanel).name) || ...
            string(panelsIn(iPanel).group) ~= string(panelsOut(iPanel).group)
        error('computePanelBucklingBayAveraged:PanelNameMismatch', ...
            'Panel order/name mismatch between bay stations at index %d.', iPanel);
    end
end
end

function b_in = panelWidth(panel)
if isfield(panel, 'ds_in') && isfinite(panel.ds_in) && panel.ds_in > 0
    b_in = panel.ds_in;
    return;
end
b_in = hypot(panel.xEnd_in - panel.xStart_in, panel.zEnd_in - panel.zStart_in);
if ~isfinite(b_in) || b_in <= 0
    error('computePanelBucklingBayAveraged:InvalidPanelWidth', ...
        'Panel %s has invalid width.', panel.name);
end
end

function [sigma_cr_c_effective_psi, correction] = applyCurvedPanelCorrection( ...
    panels, sigma_cr_c_flat_psi, Ec_psi, bayLength_in, wing, analysisOptions)
nPanels = numel(panels);
enabled = isfield(analysisOptions, 'panelBuckling') && ...
    isfield(analysisOptions.panelBuckling, 'enableCurvedPanelCorrection') && ...
    logical(analysisOptions.panelBuckling.enableCurvedPanelCorrection);

sigma_cr_c_effective_psi = sigma_cr_c_flat_psi;
correction = struct();
correction.enabled = enabled;
correction.method = "disabled";
correction.applyTo = "none";
if enabled
    correction.method = string(analysisOptions.panelBuckling.curvedPanelCorrectionMethod);
    correction.applyTo = string(analysisOptions.panelBuckling.curvedPanelApplyTo);
end
correction.isCurvedPanel = false(1, nPanels);
correction.applied = false(1, nPanels);
correction.regime = strings(1, nPanels);
correction.R_in = NaN(1, nPanels);
correction.R_over_t = NaN(1, nPanels);
correction.L_over_R = NaN(1, nPanels);
correction.deltaTheta_rad = NaN(1, nPanels);
correction.sigma_cr_c_flat_psi = sigma_cr_c_flat_psi;
correction.sigma_cr_cylinder_psi = zeros(1, nPanels);
correction.sigma_cr_c_effective_psi = sigma_cr_c_effective_psi;
correction.warning = strings(1, nPanels);
correction.counts.curvature_valid_range = 0;
correction.counts.curvature_high_R_over_t = 0;
correction.counts.flat_only_outside_validity = 0;
correction.counts.not_curved_panel = 0;
correction.counts.fallback_invalid_geometry = 0;

if ~enabled
    correction.regime(:) = "disabled";
    correction.counts.not_curved_panel = nPanels;
    return;
end

if ~(isfield(wing, 'planform') && isfield(wing.planform, 'chord_in') && ...
        isfinite(wing.planform.chord_in) && wing.planform.chord_in > 0)
    error('computePanelBucklingBayAveraged:MissingChord', ...
        'wing.planform.chord_in is required for curved-panel correction.');
end

chord_in = wing.planform.chord_in;
for iPanel = 1:nPanels
    panel = panels(iPanel);
    correction.isCurvedPanel(iPanel) = isCurvedAirfoilPanel(panel);
    if ~correction.isCurvedPanel(iPanel)
        correction.regime(iPanel) = "not_curved_panel";
        correction.counts.not_curved_panel = correction.counts.not_curved_panel + 1;
        continue;
    end

    [R_in, deltaTheta_rad, warningText] = estimatePanelRadius(panel, chord_in);
    correction.R_in(iPanel) = R_in;
    correction.deltaTheta_rad(iPanel) = deltaTheta_rad;
    correction.warning(iPanel) = warningText;

    t_in = panel.thickness_in;
    if ~(isfinite(R_in) && R_in > 0 && isfinite(t_in) && t_in > 0 && ...
            isfinite(bayLength_in) && bayLength_in > 0)
        correction.regime(iPanel) = "flat_only_invalid_geometry";
        correction.counts.fallback_invalid_geometry = ...
            correction.counts.fallback_invalid_geometry + 1;
        continue;
    end

    R_over_t = R_in / t_in;
    L_over_R = bayLength_in / R_in;
    correction.R_over_t(iPanel) = R_over_t;
    correction.L_over_R(iPanel) = L_over_R;

    sigmaCylinder_psi = 0.0;
    if R_over_t > 500 && R_over_t < 3000 && L_over_R > 0.1 && L_over_R < 2.5
        sigmaCylinder_psi = Ec_psi * (9.0 * (t_in / R_in)^1.6 + ...
            0.16 * (t_in / bayLength_in)^1.2);
        correction.regime(iPanel) = "curvature_valid_range";
        correction.applied(iPanel) = true;
        correction.counts.curvature_valid_range = ...
            correction.counts.curvature_valid_range + 1;
    elseif R_over_t >= 3000
        sigmaCylinder_psi = 0.3 * Ec_psi * (t_in / R_in);
        correction.regime(iPanel) = "curvature_high_R_over_t";
        correction.applied(iPanel) = true;
        correction.counts.curvature_high_R_over_t = ...
            correction.counts.curvature_high_R_over_t + 1;
    else
        correction.regime(iPanel) = "flat_only_outside_validity";
        correction.counts.flat_only_outside_validity = ...
            correction.counts.flat_only_outside_validity + 1;
    end

    correction.sigma_cr_cylinder_psi(iPanel) = sigmaCylinder_psi;
    sigma_cr_c_effective_psi(iPanel) = sigma_cr_c_flat_psi(iPanel) + sigmaCylinder_psi;
end

correction.sigma_cr_c_effective_psi = sigma_cr_c_effective_psi;
end

function tf = isCurvedAirfoilPanel(panel)
panelName = lower(string(panel.name));
panelGroup = lower(string(panel.group));
tf = contains(panelGroup, "top") || contains(panelGroup, "bottom") || ...
    contains(panelGroup, "uppernose") || contains(panelGroup, "lowernose") || ...
    startsWith(panelName, "top_") || startsWith(panelName, "bottom_") || ...
    startsWith(panelName, "uppernose_") || startsWith(panelName, "lowernose_");
end

function [R_in, deltaTheta_rad, warningText] = estimatePanelRadius(panel, chord_in)
warningText = "";
R_in = NaN;
deltaTheta_rad = NaN;

if ~(isfield(panel, 'surfaceFlag') && isfinite(panel.surfaceFlag))
    warningText = "Curved panel has no finite surfaceFlag; flat-panel compression critical stress retained.";
    return;
end
if ~(isfield(panel, 'ds_in') && isfinite(panel.ds_in) && panel.ds_in > 0)
    warningText = "Curved panel has no finite ds_in; flat-panel compression critical stress retained.";
    return;
end

xStart_c = panel.xStart_in / chord_in;
xEnd_c = panel.xEnd_in / chord_in;
upperFlag = panel.surfaceFlag >= 0.5;
theta1 = naca2412SurfaceSlopeAngle(xStart_c, upperFlag);
theta2 = naca2412SurfaceSlopeAngle(xEnd_c, upperFlag);
deltaTheta_rad = abs(theta2 - theta1);
if ~(isfinite(deltaTheta_rad) && deltaTheta_rad > 1.0e-8)
    warningText = "Panel slope change is near zero; flat-panel compression critical stress retained.";
    return;
end

R_in = panel.ds_in / deltaTheta_rad;
if ~(isfinite(R_in) && R_in > 0)
    warningText = "Estimated panel radius is nonfinite; flat-panel compression critical stress retained.";
end
end

function theta = naca2412SurfaceSlopeAngle(x_c, upperFlag)
h = 1.0e-5;
x_c = min(max(x_c, 1.0e-6), 1.0);
x1 = max(1.0e-6, x_c - h);
x2 = min(1.0, x_c + h);
if x2 <= x1
    x1 = max(1.0e-6, x_c - 2.0 * h);
    x2 = min(1.0, x_c);
end
z1 = naca2412Z(x1, upperFlag);
z2 = naca2412Z(x2, upperFlag);
theta = atan((z2 - z1) / (x2 - x1));
end

function z = naca2412Z(x, upperFlag)
M = 0.02;
P = 0.4;
T = 0.12;
a0 = 0.2969;
a1 = -0.126;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

if x < P
    z_camber = M / P^2 * (2 * P * x - x^2);
else
    z_camber = (M / (1 - P)^2) * (1 - 2 * P + 2 * P * x - x^2);
end
z_thickness = (T / 0.2) * (a0 * sqrt(x) + a1 * x + a2 * x^2 + ...
    a3 * x^3 + a4 * x^4);

if upperFlag
    z = z_camber + z_thickness;
else
    z = z_camber - z_thickness;
end
end

function [idx, distance] = findNearestBoom(x_in, z_in, entries)
dx = [entries.x_in] - x_in;
dz = [entries.z_in] - z_in;
distances = hypot(dx, dz);
[distance, idx] = min(distances);
if isempty(idx) || ~isfinite(distance)
    error('computePanelBucklingBayAveraged:NoUsableBoomMatch', ...
        'Could not find a usable boom match for panel endpoint.');
end
end

function [RActive, modes, warnings] = assignCombinedRatios(R_SC, R_BS, panels, caseNames, analysisOptions)
RActive = zeros(size(R_SC));
modes = strings(size(R_SC));
warnings = strings(0, 1);

if ~getCombinedEnabled(analysisOptions)
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

function match = makePanelMatch(panel, panelIndex, startIdxIn, endIdxIn, startIdxOut, endIdxOut, ...
    startDistIn, endDistIn, startDistOut, endDistOut, startEntryIn, endEntryIn, startEntryOut, endEntryOut)
if nargin == 0
    panel = makePanelPlaceholder();
    panelIndex = NaN;
    startIdxIn = NaN; endIdxIn = NaN; startIdxOut = NaN; endIdxOut = NaN;
    startDistIn = NaN; endDistIn = NaN; startDistOut = NaN; endDistOut = NaN;
    startEntryIn = makeEntryPlaceholder();
    endEntryIn = makeEntryPlaceholder();
    startEntryOut = makeEntryPlaceholder();
    endEntryOut = makeEntryPlaceholder();
end
match = struct();
match.panel_index = panelIndex;
match.panel_name = panel.name;
match.panel_group = panel.group;
match.inboard_start_boom_index = startIdxIn;
match.inboard_end_boom_index = endIdxIn;
match.outboard_start_boom_index = startIdxOut;
match.outboard_end_boom_index = endIdxOut;
match.inboard_start_match_distance_in = startDistIn;
match.inboard_end_match_distance_in = endDistIn;
match.outboard_start_match_distance_in = startDistOut;
match.outboard_end_match_distance_in = endDistOut;
match.inboard_start_boom_name = startEntryIn.name;
match.inboard_end_boom_name = endEntryIn.name;
match.outboard_start_boom_name = startEntryOut.name;
match.outboard_end_boom_name = endEntryOut.name;
match.method = "nearest panel endpoint boom match at inboard and outboard bay stations";
end

function panel = makePanelPlaceholder()
panel = struct();
panel.name = "";
panel.group = "";
end

function entry = makeEntryPlaceholder()
entry = struct();
entry.name = "";
end
