function divergence = computeDivergence(wing, section, loads)
% computeDivergence Wing divergence screening check.
%
% Inputs:
%   wing    - Geometry/material struct from setGeometry.
%   section - Struct from computeSectionUnitLoads.
%   loads   - Optional loads struct, reserved for future VD metadata.
%
% Outputs:
%   divergence - Struct containing torsional stiffness and
%                divergence speed diagnostics.
%
% Scope:
%   Divergence screening check for the current wing-box model.
%
% Formula:
%   V_divergence = sqrt((pi^2*G*J) / ...
%       (2*rho*CL_alpha*e*c^2*semiSpan^2))
%
% Units:
%   G in lb/in^2, J in in^4, rho in lb*s^2/in^4, c and semiSpan in inches,
%   and V_divergence in in/s.
%
% Excluded:
%   buckling, optimisation, final failure margins, yield-screening changes,
%   and weight changes.

if nargin < 3
    loads = [];
end

validateInputs(wing, section);

G_used_psi = wing.material.skinWeb.G_psi;
chord_in = wing.planform.chord_in;
semispan_in = wing.planform.semispan_in;
x_sc_in = section.shearCenter.x_sc_in;
x_ac_in = 0.25 * chord_in;
e_nondim = abs((x_sc_in - x_ac_in) / chord_in);

CL_alpha_per_rad = 4.5942;
rho_slug_ft3 = 0.00238;
rho_lbs2_in4 = rho_slug_ft3 / 12^4;
VD_mph = 270;

torsionEstimate = estimateTorsionConstant(section);
J_prelim_in4 = torsionEstimate.J_prelim_in4;

V_divergence_in_s = sqrt((pi^2 * G_used_psi * J_prelim_in4) / ...
    (2 * rho_lbs2_in4 * CL_alpha_per_rad * e_nondim * chord_in^2 * semispan_in^2));
V_divergence_ft_s = V_divergence_in_s / 12.0;
V_divergence_mph = V_divergence_ft_s * 3600.0 / 5280.0;
V_divergence_over_VD = V_divergence_mph / VD_mph;
divergencePass = V_divergence_mph > VD_mph;

validateOutputs(G_used_psi, J_prelim_in4, e_nondim, rho_lbs2_in4, ...
    CL_alpha_per_rad, V_divergence_in_s, torsionEstimate);

divergence = struct();
divergence.J_prelim_in4 = J_prelim_in4;
divergence.G_used_psi = G_used_psi;
divergence.CL_alpha_per_rad = CL_alpha_per_rad;
divergence.rho_slug_ft3 = rho_slug_ft3;
divergence.rho_lbs2_in4 = rho_lbs2_in4;
divergence.e_nondim = e_nondim;
divergence.x_sc_in = x_sc_in;
divergence.x_ac_in = x_ac_in;
divergence.chord_in = chord_in;
divergence.semispan_in = semispan_in;
divergence.VD_mph = VD_mph;
divergence.V_divergence_in_s = V_divergence_in_s;
divergence.V_divergence_ft_s = V_divergence_ft_s;
divergence.V_divergence_mph = V_divergence_mph;
divergence.V_divergence_over_VD = V_divergence_over_VD;
divergence.pass = divergencePass;
divergence.torsion = torsionEstimate;
divergence.inputs.loadsProvided = ~isempty(loads);
divergence.units.G = "lb/in^2";
divergence.units.J = "in^4";
divergence.units.rho = "lb*s^2/in^4";
divergence.units.velocity = ["in/s"; "ft/s"; "mph"];
divergence.method = ...
    "Divergence formula with J estimated from the current two-cell thin-walled torsion model.";
divergence.assumptions = [
    "Divergence screening check for the current wing-box model."
    "Formula follows the project aeroelastic screening convention."
    "J is estimated from the current two-cell thin-walled torsion model."
    "G uses wing.material.skinWeb.G_psi for the thin-walled skin/web torsion model."
    "CL_alpha_per_rad = 4.5942 and rho_slug_ft3 = 0.00238 are carried from the saved aerodynamic load assumptions."
    "e = abs((x_sc - 0.25*c)/c)."
    "VD_mph = 270."
    ];
divergence.warnings = [
    "Divergence uses the current two-cell thin-walled torsion model and remains a screening-level aeroelastic check."
    "Material properties are read from the verified project material table."
    "No final aeroelastic margin is claimed."
    "This check does not include buckling, optimisation, or final failure margins."
    ];
divergence.status = "divergence_screening_check";
divergence.exclusions = [
    "buckling"
    "optimisation"
    "final failure margins"
    "yield-screening changes"
    "weight changes"
    ];
end

function validateInputs(wing, section)
if nargin < 1 || isempty(wing) || ~isstruct(wing)
    error('computeDivergence:MissingWing', 'wing struct from setGeometry is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computeDivergence:MissingSection', 'section struct from computeSectionUnitLoads is required.');
end

if ~isfield(wing, 'material') || ~isfield(wing.material, 'skinWeb') || ...
        ~isfield(wing.material.skinWeb, 'G_psi')
    error('computeDivergence:MissingShearModulus', ...
        'wing.material.skinWeb.G_psi is required.');
end
if ~(isfinite(wing.material.skinWeb.G_psi) && wing.material.skinWeb.G_psi > 0)
    error('computeDivergence:InvalidShearModulus', ...
        'wing.material.skinWeb.G_psi must be positive and finite.');
end
if ~isfield(wing, 'planform') || ~isfield(wing.planform, 'chord_in') || ...
        ~isfield(wing.planform, 'semispan_in')
    error('computeDivergence:MissingPlanform', ...
        'wing.planform.chord_in and semispan_in are required.');
end
if ~(isfinite(wing.planform.chord_in) && wing.planform.chord_in > 0 && ...
        isfinite(wing.planform.semispan_in) && wing.planform.semispan_in > 0)
    error('computeDivergence:InvalidPlanform', ...
        'Chord and semispan must be positive and finite.');
end
if ~isfield(section, 'shearCenter') || ~isfield(section.shearCenter, 'x_sc_in') || ...
        ~isfinite(section.shearCenter.x_sc_in)
    error('computeDivergence:MissingShearCenter', ...
        'section.shearCenter.x_sc_in is required and must be finite.');
end
if ~isfield(section, 'unitShear') || ~isfield(section.unitShear, 'stage3B') || ...
        ~isfield(section.unitShear.stage3B, 'compatibility') || ...
        ~isfield(section.unitShear.stage3B, 'preTorsion') || ...
        ~isfield(section.unitShear, 'stage3C') || ...
        ~isfield(section.unitShear.stage3C, 'torsion')
    error('computeDivergence:MissingTorsionData', ...
        'Shear-flow compatibility/preTorsion and torsion data are required.');
end
end

function torsionEstimate = estimateTorsionConstant(section)
compatibility = section.unitShear.stage3B.compatibility;
preTorsion = section.unitShear.stage3B.preTorsion;
torsion = section.unitShear.stage3C.torsion;

Area1 = preTorsion.area.Area1_mainCell_in2;
Area2 = preTorsion.area.Area2_noseCell_in2;
A11 = compatibility.A11;
A22 = compatibility.A22;
shared_ds_over_t = compatibility.sharedWall.mainBox_ds_over_t;
q1 = torsion.qt_Y_mainCell_per_lbin;
q2 = torsion.qt_Y_noseCell_per_lbin;

required = [Area1, Area2, A11, A22, shared_ds_over_t, q1, q2];
if any(~isfinite(required)) || Area1 <= 0 || Area2 <= 0 || A11 <= 0 || A22 <= 0 || shared_ds_over_t <= 0
    error('computeDivergence:InvalidTorsionInputs', ...
        'Two-cell torsion inputs must be positive and finite, except signed q values which must be finite.');
end

% This follows the same shared-wall sign convention used in
% computeSectionUnitLoads pure-torsion unit solution. For unit torque, these no-G
% compliances are equal to 1/J for each cell if compatibility is satisfied.
complianceCell1_noG_per_in4 = (q1 * A11 - q2 * shared_ds_over_t) / (2 * Area1);
complianceCell2_noG_per_in4 = (q2 * A22 - q1 * shared_ds_over_t) / (2 * Area2);
compatibilityAbsDiff = abs(complianceCell1_noG_per_in4 - complianceCell2_noG_per_in4);
compatibilityRelDiff = compatibilityAbsDiff / max(1e-12, ...
    abs((complianceCell1_noG_per_in4 + complianceCell2_noG_per_in4) / 2));
averageCompliance_noG_per_in4 = mean([complianceCell1_noG_per_in4, complianceCell2_noG_per_in4]);
J_prelim_in4 = 1 / averageCompliance_noG_per_in4;

torsionEstimate = struct();
torsionEstimate.J_prelim_in4 = J_prelim_in4;
torsionEstimate.complianceCell1_noG_per_in4 = complianceCell1_noG_per_in4;
torsionEstimate.complianceCell2_noG_per_in4 = complianceCell2_noG_per_in4;
torsionEstimate.averageCompliance_noG_per_in4 = averageCompliance_noG_per_in4;
torsionEstimate.compatibilityAbsDiff = compatibilityAbsDiff;
torsionEstimate.compatibilityRelDiff = compatibilityRelDiff;
torsionEstimate.compatibilityTolerance = 1e-6;
torsionEstimate.compatibilityPass = compatibilityRelDiff <= torsionEstimate.compatibilityTolerance;
torsionEstimate.Area1_mainCell_in2 = Area1;
torsionEstimate.Area2_noseCell_in2 = Area2;
torsionEstimate.A11 = A11;
torsionEstimate.A22 = A22;
torsionEstimate.shared_ds_over_t = shared_ds_over_t;
torsionEstimate.qt_Y_mainCell_per_lbin = q1;
torsionEstimate.qt_Y_noseCell_per_lbin = q2;
torsionEstimate.method = ...
    "J = 1 / average no-G twist compliance from the pure My=1 two-cell torsion solution.";
end

function validateOutputs(G_used_psi, J_prelim_in4, e_nondim, rho_lbs2_in4, ...
    CL_alpha_per_rad, V_divergence_in_s, torsionEstimate)
if ~(isfinite(G_used_psi) && G_used_psi > 0)
    error('computeDivergence:InvalidG', 'G must be positive and finite.');
end
if ~(isfinite(J_prelim_in4) && J_prelim_in4 > 0)
    error('computeDivergence:InvalidJ', 'J must be positive and finite.');
end
if ~(isfinite(e_nondim) && e_nondim > 0)
    error('computeDivergence:InvalidE', ...
        'e must be positive and finite; shear center cannot coincide with quarter chord for this formula.');
end
if ~(isfinite(rho_lbs2_in4) && rho_lbs2_in4 > 0)
    error('computeDivergence:InvalidDensity', 'rho must be positive and finite.');
end
if ~(isfinite(CL_alpha_per_rad) && CL_alpha_per_rad > 0)
    error('computeDivergence:InvalidCLAlpha', 'CL_alpha must be positive and finite.');
end
if ~(isreal(V_divergence_in_s) && isfinite(V_divergence_in_s) && V_divergence_in_s > 0)
    error('computeDivergence:InvalidDivergenceSpeed', ...
        'V_divergence must be finite, real, and positive.');
end
if ~torsionEstimate.compatibilityPass
    error('computeDivergence:TorsionCompatibilityMismatch', ...
        'Two-cell twist compatibility relative difference %.6g exceeds tolerance %.6g.', ...
        torsionEstimate.compatibilityRelDiff, torsionEstimate.compatibilityTolerance);
end
end

