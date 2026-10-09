clear;clc; close all;

%% ============================================================
%  Flight-envelope and spanwise wing loads (V2: tip-inward sign fix)
%  Utility-category aircraft, NACA 2412 rectangular wing
%
%  This script:
%  1)  Builds the complete V-n diagram with gust lines
%  2)  Evaluates all maneuver-envelope corner cases:
%        PHAA, PLAA, NHAA, NLAA
%  3)  Evaluates gust cases at VC and VD
%  4)  Identifies governing positive and negative vertical cases
%  5)  Computes spanwise load distributions for ALL nine cases:
%        qz(y), Vz(y), Mx(y)   -- vertical bending
%        qx(y), Vx(y), Mz(y)   -- in-plane (drag) bending
%        qm(y), My(y)           -- torsion
%  6)  Saves two figures per case to outputs/figures/load_cases/:
%        qvm_<case>.png         -- vertical: qz, Vz, Mx
%        drag_torsion_<case>.png -- in-plane: qx, Vx, Mz, qm, My
%  7)  Prints a consolidated summary table for all cases
%  8)  Saves the V-n diagram as outputs/figures/load_cases/vn_diagram.png
%
%  Governing cases for structural sizing remain PHAA (+) and
%  Negative gust at VC (-).  All other cases are computed and
%  plotted for completeness.
%
%  NOTE: My(y) here represents only the aerodynamic pitching-moment
%  contribution.  The torsional contribution from Fz acting at a
%  moment arm from the shear center is added in the structural analysis once the
%  wing-box geometry (and hence shear-center location) is known.
%% ============================================================

%% Create output folder for figures
figFolder = fullfile(fileparts(mfilename('fullpath')), '..', 'outputs', 'figures', 'load_cases');
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

%% 1) GIVEN DATA
rho = 0.00238;          % air density [slug/ft^3]
g   = 32.174;           % gravity [ft/s^2]

% Geometry
c = 5.0;                % chord [ft]
cbar = c;               % mean aerodynamic chord [ft]
semiSpan = 17.0;        % semispan [ft]
b = 2*semiSpan;         % full span [ft]
S = b*c;                % wing area [ft^2]
AR = b^2/S;             % aspect ratio

% Weights
Wmax = 3200.0;          % [lb]
Wmin = 1900.0;          % [lb]

% Speeds
Vc_mph = 230.0;         % [mph]
Vd_mph = 270.0;         % [mph]
mph_to_fts = 1.46667;
Vc = Vc_mph * mph_to_fts;   % [ft/s]
Vd = Vd_mph * mph_to_fts;   % [ft/s]

% FAR 23 utility category limit load factors
n_max = +4.4;
n_min = -1.76;

% Finite-wing efficiency factor
e = 0.8;

%% 2) XFOIL-BASED DATA USED IN THE REPORT
% High-alpha positive peak
Cl_max_2D_raw   = 1.8523;
alpha_pos_raw_d = 19.5;      % [deg]
Cd_pos_raw      = 0.04610;
% CM for PHAA: XFOIL CM interpolated at alpha_PHAA = 1.25*19.5 = 24.375 deg
% between the 24.0 deg (-0.0452) and 24.5 deg (-0.0553) polar rows (-0.052775),
% rounded to 4 decimals. Not the CM at the 19.5 deg peak (-0.0071).
Cm_pos_raw      = -0.0528;

% High-alpha negative peak
Cl_min_2D       = -1.6245;
alpha_neg_d     = -19.0;     % [deg]
Cd_neg          = 0.02987;
Cm_neg          = -0.0716;

% Additional positive high-alpha CD data used for PHAA interpolation
alpha_cd_hi = [23.5, 24.0, 24.5, 25.0];                % [deg]
cd_hi_table = [0.13068, 0.14669, 0.16283, 0.18161];

%% 3) FINITE-WING CORRECTION
CL_max_3D_raw = Cl_max_2D_raw / (1 + Cl_max_2D_raw/(pi*AR*e));
CL_min_3D     = Cl_min_2D     / (1 + abs(Cl_min_2D)/(pi*AR*e));

% 3D lift-curve slope for low-alpha cases and gusts
Cl_alpha_2D = 2*pi;   % [1/rad]
CL_alpha_3D = Cl_alpha_2D / (1 + Cl_alpha_2D/(pi*AR*e));   % [1/rad]

%% 4) PHAA SPECIAL +25% MARGIN
alpha_PHAA_d = 1.25 * alpha_pos_raw_d;                 % [deg]
alpha_PHAA   = deg2rad(alpha_PHAA_d);                  % [rad]
CL_PHAA      = 1.25 * CL_max_3D_raw;
CD_PHAA      = interp1(alpha_cd_hi, cd_hi_table, alpha_PHAA_d, 'linear', 'extrap');
CM_PHAA      = Cm_pos_raw;

CN_PHAA = CL_PHAA*cos(alpha_PHAA) + CD_PHAA*sin(alpha_PHAA);
CX_PHAA = CD_PHAA*cos(alpha_PHAA) - CL_PHAA*sin(alpha_PHAA);

%% 5) NHAA SPECIAL HIGH-ALPHA CASE
alpha_NHAA = deg2rad(alpha_neg_d);                     % [rad]
CL_NHAA    = CL_min_3D;
CD_NHAA    = Cd_neg;
CM_NHAA    = Cm_neg;

CN_NHAA = CL_NHAA*cos(alpha_NHAA) + CD_NHAA*sin(alpha_NHAA);
CX_NHAA = CD_NHAA*cos(alpha_NHAA) - CL_NHAA*sin(alpha_NHAA);

%% 6) LOW-ALPHA MODEL FOR PLAA, NLAA AND GUST POINTS
alpha_L0_d = -2.13;                          % [deg]
alpha_L0   = deg2rad(alpha_L0_d);            % [rad]

cd_fit_alpha = [-3.95, -2.88, 3.40];
cd_fit_vals  = [0.00606, 0.00581, 0.00538];
cd_poly = polyfit(cd_fit_alpha, cd_fit_vals, 2);

CM_low = -0.0520;

%% 7) MANEUVERING SPEEDS
VA_pos = sqrt(2*Wmax*n_max/(rho*S*CN_PHAA));           % [ft/s]
VA_neg = sqrt(2*Wmin*abs(n_min)/(rho*S*abs(CN_NHAA))); % [ft/s]

VA_pos_mph = VA_pos / mph_to_fts;
VA_neg_mph = VA_neg / mph_to_fts;

%% 8) V-n CURVES
Vplot_mph = linspace(0, 300, 2000);
Vplot = Vplot_mph * mph_to_fts;                       % [ft/s]

n_pos_curve = 0.5*rho*Vplot.^2*S*CN_PHAA/Wmax;
n_neg_curve = 0.5*rho*Vplot.^2*S*CN_NHAA/Wmin;

mask_pos = Vplot <= VA_pos;
mask_neg = Vplot <= VA_neg;

V_pos_curve_mph = Vplot_mph(mask_pos);
n_pos_env       = n_pos_curve(mask_pos);

V_neg_curve_mph = Vplot_mph(mask_neg);
n_neg_env       = n_neg_curve(mask_neg);

% Negative branch prescribed by FAR 23
V_neg_line1_mph = [VA_neg_mph, Vc_mph];
n_neg_line1     = [n_min, n_min];

V_neg_line2_mph = [Vc_mph, Vd_mph];
n_neg_line2     = [n_min, -1.0];

%% 9) GUST LINES
Ude_Vc = 50.0;   % [ft/s]
Ude_Vd = 25.0;   % [ft/s]

W_gust = Wmin;

mu = 2*(W_gust/S) / (rho*g*cbar*CL_alpha_3D);
k  = 0.88*mu / (5.3 + mu);

delta_n_Vc = rho*S*CL_alpha_3D*Vc*Ude_Vc*k/(2*W_gust);
delta_n_Vd = rho*S*CL_alpha_3D*Vd*Ude_Vd*k/(2*W_gust);

n_gust_pos_Vc = 1 + delta_n_Vc;
n_gust_neg_Vc = 1 - delta_n_Vc;

n_gust_pos_Vd = 1 + delta_n_Vd;
n_gust_neg_Vd = 1 - delta_n_Vd;

Vg1_mph = linspace(0, Vc_mph, 500);
Vg1 = Vg1_mph * mph_to_fts;
delta_n_g1 = rho*S*CL_alpha_3D.*Vg1*Ude_Vc*k/(2*W_gust);
n_gust_pos_1 = 1 + delta_n_g1;
n_gust_neg_1 = 1 - delta_n_g1;

Vg2_mph = linspace(Vc_mph, Vd_mph, 300);
n_gust_pos_2 = linspace(n_gust_pos_Vc, n_gust_pos_Vd, numel(Vg2_mph));
n_gust_neg_2 = linspace(n_gust_neg_Vc, n_gust_neg_Vd, numel(Vg2_mph));

%% 10) CHARACTERISTIC MANEUVER POINTS
PHAA_pt = [VA_pos_mph, n_max];
NHAA_pt = [VA_neg_mph, n_min];
PLAA_pt = [Vd_mph,     n_max];
NLAA_pt = [Vd_mph,    -1.0];

PG_VC_pt = [Vc_mph, n_gust_pos_Vc];
NG_VC_pt = [Vc_mph, n_gust_neg_Vc];
PG_VD_pt = [Vd_mph, n_gust_pos_Vd];
NG_VD_pt = [Vd_mph, n_gust_neg_Vd];

%% 11) EVALUATE ALL LOAD CASES

% ---- PHAA ----
PHAA.name     = 'PHAA';
PHAA.V_mph    = VA_pos_mph;
PHAA.V        = VA_pos;
PHAA.n        = n_max;
PHAA.W        = Wmax;
PHAA.alpha_d  = alpha_PHAA_d;
PHAA.alpha    = alpha_PHAA;
PHAA.CL       = CL_PHAA;
PHAA.CD       = CD_PHAA;
PHAA.CM       = CM_PHAA;
PHAA.CN       = CN_PHAA;
PHAA.CX       = CX_PHAA;
PHAA.L        = 0.5*rho*PHAA.V^2*S*PHAA.CL;
PHAA.D        = 0.5*rho*PHAA.V^2*S*PHAA.CD;
PHAA.Fz       = 0.5*rho*PHAA.V^2*S*PHAA.CN;
PHAA.Fx       = 0.5*rho*PHAA.V^2*S*PHAA.CX;
PHAA.My       = 0.5*rho*PHAA.V^2*S*cbar*PHAA.CM;

% ---- NHAA ----
NHAA.name     = 'NHAA';
NHAA.V_mph    = VA_neg_mph;
NHAA.V        = VA_neg;
NHAA.n        = n_min;
NHAA.W        = Wmin;
NHAA.alpha_d  = alpha_neg_d;
NHAA.alpha    = alpha_NHAA;
NHAA.CL       = CL_NHAA;
NHAA.CD       = CD_NHAA;
NHAA.CM       = CM_NHAA;
NHAA.CN       = CN_NHAA;
NHAA.CX       = CX_NHAA;
NHAA.L        = 0.5*rho*NHAA.V^2*S*NHAA.CL;
NHAA.D        = 0.5*rho*NHAA.V^2*S*NHAA.CD;
NHAA.Fz       = 0.5*rho*NHAA.V^2*S*NHAA.CN;
NHAA.Fx       = 0.5*rho*NHAA.V^2*S*NHAA.CX;
NHAA.My       = 0.5*rho*NHAA.V^2*S*cbar*NHAA.CM;

% ---- PLAA ----
PLAA = solve_low_alpha_case('PLAA', Vd_mph, n_max, Wmax, rho, S, cbar, ...
                            CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

% ---- NLAA at VC ----
NLAA_VC = solve_low_alpha_case('NLAA_VC', Vc_mph, n_min, Wmin, rho, S, cbar, ...
                               CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

% ---- NLAA at VD ----
NLAA_VD = solve_low_alpha_case('NLAA_VD', Vd_mph, -1.0, Wmin, rho, S, cbar, ...
                               CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

% ---- Gust cases ----
PG_VC = solve_low_alpha_case('PG_VC', Vc_mph, n_gust_pos_Vc, Wmin, rho, S, cbar, ...
                             CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

NG_VC = solve_low_alpha_case('NG_VC', Vc_mph, n_gust_neg_Vc, Wmin, rho, S, cbar, ...
                             CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

PG_VD = solve_low_alpha_case('PG_VD', Vd_mph, n_gust_pos_Vd, Wmin, rho, S, cbar, ...
                             CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

NG_VD = solve_low_alpha_case('NG_VD', Vd_mph, n_gust_neg_Vd, Wmin, rho, S, cbar, ...
                             CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts);

%% 12) GOVERNING VERTICAL CASES
allCases = [PHAA, NHAA, PLAA, NLAA_VC, NLAA_VD, PG_VC, NG_VC, PG_VD, NG_VD];
Fz_all = [allCases.Fz];

[~, idx_pos] = max(Fz_all);
[~, idx_neg] = min(Fz_all);

GOV_POS = allCases(idx_pos);
GOV_NEG = allCases(idx_neg);

%% 13) SPANWISE GRID
y = linspace(0, semiSpan, 800);     % root (y=0) to tip (y=17 ft)

%% ============================================================
%  Compute load distributions for ALL nine cases
%  and save two figures per case.
%% ============================================================

nCases = numel(allCases);

% Pre-allocate cell arrays to hold distribution structs and
% file-safe case name strings for figure naming.
distZ  = cell(nCases, 1);   % vertical distributions
distX  = cell(nCases, 1);   % in-plane (drag) distributions
distMy = cell(nCases, 1);   % torsion distributions

for iCase = 1:nCases
    cas = allCases(iCase);

    % --- Vertical (Fz) distributions ---
    distZ{iCase} = build_vertical_distribution(cas.Fz/2, semiSpan, y);

    % --- In-plane (Fx / drag) distributions ---
    distX{iCase} = build_drag_distribution(cas.Fx/2, semiSpan, y);

    % --- Torsion (pitching moment) distributions ---
    distMy{iCase} = build_uniform_distribution(cas.My/2, semiSpan, y);

    % Build a file-safe name (spaces replaced with underscores)
    safeName = strrep(cas.name, ' ', '_');
    safeName = strrep(safeName, '@', 'at');

    % --- Figure A: vertical qz, Vz, Mx ---
    figA = plot_qVM(y, distZ{iCase}, cas.name);
    fNameA = fullfile(figFolder, ['qvm_' safeName '.png']);
    exportgraphics(figA, fNameA, 'Resolution', 300);

    % --- Figure B: in-plane and torsion qx, Vx, Mz, qm, My ---
    figB = plot_drag_torsion(y, distX{iCase}, distMy{iCase}, cas.name);
    fNameB = fullfile(figFolder, ['drag_torsion_' safeName '.png']);
    exportgraphics(figB, fNameB, 'Resolution', 300);
end

% Identify governing-case distribution structs by the same indices used above
dist_pos  = distZ{idx_pos};
dist_neg  = distZ{idx_neg};
drag_pos  = distX{idx_pos};
drag_neg  = distX{idx_neg};
my_pos    = distMy{idx_pos};
my_neg    = distMy{idx_neg};

%% 14) V-n DIAGRAM (unchanged)
fig1 = figure('Name','Complete V-n Diagram','Color','w');
hold on; grid on; box on;

plot(V_pos_curve_mph, n_pos_env, 'b', 'LineWidth', 2);
plot([VA_pos_mph, Vd_mph], [n_max, n_max], 'b', 'LineWidth', 2);

plot(V_neg_curve_mph, n_neg_env, 'b', 'LineWidth', 2);
plot(V_neg_line1_mph, n_neg_line1, 'b', 'LineWidth', 2);
plot(V_neg_line2_mph, n_neg_line2, 'b', 'LineWidth', 2);

plot([Vd_mph, Vd_mph], [-1.0, n_max], 'b', 'LineWidth', 2);

plot(Vg1_mph, n_gust_pos_1, '--m', 'LineWidth', 1.8);
plot(Vg2_mph, n_gust_pos_2, '--m', 'LineWidth', 1.8);
plot(Vg1_mph, n_gust_neg_1, '--m', 'LineWidth', 1.8);
plot(Vg2_mph, n_gust_neg_2, '--m', 'LineWidth', 1.8);

xline(VA_pos_mph, '--r', 'V_A^+', 'LabelVerticalAlignment','bottom');
xline(VA_neg_mph, '--r', 'V_A^-', 'LabelVerticalAlignment','bottom');
xline(Vc_mph, ':k', 'V_C', 'LabelVerticalAlignment','bottom');
xline(Vd_mph, ':k', 'V_D', 'LabelVerticalAlignment','bottom');

yline(n_max, '--k');
yline(n_min, '--k');
yline(-1.0, '--k');
yline(1.0, ':k');

plot(PHAA_pt(1), PHAA_pt(2), 'ko', 'MarkerFaceColor', 'k');
plot(NHAA_pt(1), NHAA_pt(2), 'ko', 'MarkerFaceColor', 'k');
plot(PLAA_pt(1), PLAA_pt(2), 'ko', 'MarkerFaceColor', 'k');
plot(NLAA_pt(1), NLAA_pt(2), 'ko', 'MarkerFaceColor', 'k');

plot(PG_VC_pt(1), PG_VC_pt(2), 'ms', 'MarkerFaceColor', 'm');
plot(NG_VC_pt(1), NG_VC_pt(2), 'ms', 'MarkerFaceColor', 'm');
plot(PG_VD_pt(1), PG_VD_pt(2), 'md', 'MarkerFaceColor', 'm');
plot(NG_VD_pt(1), NG_VD_pt(2), 'md', 'MarkerFaceColor', 'm');

text(PHAA_pt(1)+3, PHAA_pt(2)+0.15, 'PHAA');
text(NHAA_pt(1)+3, NHAA_pt(2)-0.25, 'NHAA');
text(PLAA_pt(1)-22, PLAA_pt(2)+0.15, 'PLAA');
text(NLAA_pt(1)-22, NLAA_pt(2)-0.25, 'NLAA');

text(PG_VC_pt(1)+3, PG_VC_pt(2)+0.08, 'Gust + @ V_C', 'Color', 'm');
text(NG_VC_pt(1)+3, NG_VC_pt(2)-0.18, 'Gust - @ V_C', 'Color', 'm');
text(PG_VD_pt(1)-40, PG_VD_pt(2)+0.08, 'Gust + @ V_D', 'Color', 'm');
text(NG_VD_pt(1)-40, NG_VD_pt(2)-0.18, 'Gust - @ V_D', 'Color', 'm');

xlabel('Velocity V [mph]');
ylabel('Load factor n');
title('Complete V-n Diagram with Gust Lines - Utility Aircraft');
xlim([0 300]);
ylim([-5 8]);

exportgraphics(fig1, fullfile(figFolder, 'vn_diagram.png'), 'Resolution', 300);

%% ============================================================
%  Consolidated summary table
%% ============================================================
fprintf('\n============================================================\n');
fprintf('COMPLETE LOADS ANALYSIS\n');
fprintf('============================================================\n');

fprintf('\n--- Geometry ---\n');
fprintf('Chord c                  = %.2f ft\n', c);
fprintf('Semi-span                = %.2f ft\n', semiSpan);
fprintf('Wing area S              = %.2f ft^2\n', S);
fprintf('Aspect ratio AR          = %.3f\n', AR);

fprintf('\n--- Key aerodynamic values ---\n');
fprintf('CL_alpha_3D              = %.4f 1/rad\n', CL_alpha_3D);
fprintf('alpha_L0 (low-alpha fit) = %.3f deg\n', alpha_L0_d);
fprintf('CN_PHAA                  = %.4f\n', CN_PHAA);
fprintf('CN_NHAA                  = %.4f\n', CN_NHAA);

fprintf('\n--- Maneuvering speeds ---\n');
fprintf('V_A+                     = %.2f mph\n', VA_pos_mph);
fprintf('V_A-                     = %.2f mph\n', VA_neg_mph);

fprintf('\n--- Gust lines ---\n');
fprintf('mu                       = %.4f\n', mu);
fprintf('k                        = %.4f\n', k);
fprintf('Delta n at Vc            = %.4f\n', delta_n_Vc);
fprintf('Delta n at Vd            = %.4f\n', delta_n_Vd);
fprintf('Positive gust at Vc      = %.4f\n', n_gust_pos_Vc);
fprintf('Negative gust at Vc      = %.4f\n', n_gust_neg_Vc);
fprintf('Positive gust at Vd      = %.4f\n', n_gust_pos_Vd);
fprintf('Negative gust at Vd      = %.4f\n', n_gust_neg_Vd);

fprintf('\n--- Per-case aerodynamic loads ---\n');
for iCase = 1:nCases
    print_case(allCases(iCase));
end

fprintf('\n--- Governing vertical cases ---\n');
fprintf('Positive governing case  = %s, Fz = %.1f lb\n', GOV_POS.name, GOV_POS.Fz);
fprintf('Negative governing case  = %s, Fz = %.1f lb\n', GOV_NEG.name, GOV_NEG.Fz);

% Updated semi-wing results block: now includes Vx_root, Mz_root, My_root
fprintf('\n--- Semi-wing root loads (governing cases) ---\n');
fprintf('%s:\n', GOV_POS.name);
fprintf('  Fz_wing   = %10.1f lb\n',    dist_pos.F_wing);
fprintf('  Vz_root   = %10.1f lb\n',    dist_pos.V_root);
fprintf('  Mx_root   = %10.1f lb-ft\n', dist_pos.M_root);
fprintf('  Vx_root   = %10.1f lb\n',    drag_pos.Vx_root);
fprintf('  Mz_root   = %10.1f lb-ft\n', drag_pos.Mz_root);
fprintf('  My_root   = %10.1f lb-ft\n', my_pos.My_root);

fprintf('%s:\n', GOV_NEG.name);
fprintf('  Fz_wing   = %10.1f lb\n',    dist_neg.F_wing);
fprintf('  Vz_root   = %10.1f lb\n',    dist_neg.V_root);
fprintf('  Mx_root   = %10.1f lb-ft\n', dist_neg.M_root);
fprintf('  Vx_root   = %10.1f lb\n',    drag_neg.Vx_root);
fprintf('  Mz_root   = %10.1f lb-ft\n', drag_neg.Mz_root);
fprintf('  My_root   = %10.1f lb-ft\n', my_neg.My_root);

% Consolidated summary table: one row per case
fprintf('\n--- Consolidated summary table (semi-wing root loads) ---\n');
hdr = sprintf('%-22s  %7s  %6s  %7s  %9s  %9s  %9s  %9s  %10s  %9s  %10s  %10s', ...
    'Case', 'V(mph)', 'n', 'a(deg)', ...
    'Fz(lb)', 'Fx(lb)', 'My(lb-ft)', ...
    'Vz_r(lb)', 'Mx_r(lb-ft)', ...
    'Vx_r(lb)', 'Mz_r(lb-ft)', 'My_r(lb-ft)');
fprintf('\n%s\n', hdr);
fprintf('%s\n', repmat('-', 1, length(hdr)));

for iCase = 1:nCases
    cas = allCases(iCase);
    dz  = distZ{iCase};
    dx  = distX{iCase};
    dm  = distMy{iCase};
    fprintf('%-22s  %7.2f  %6.3f  %7.3f  %9.1f  %9.1f  %9.1f  %9.1f  %10.1f  %9.1f  %10.1f  %10.1f\n', ...
        cas.name, cas.V_mph, cas.n, cas.alpha_d, ...
        cas.Fz, cas.Fx, cas.My, ...
        dz.V_root, dz.M_root, ...
        dx.Vx_root, dx.Mz_root, dm.My_root);
end

run(fullfile(fileparts(mfilename('fullpath')), 'exportLoadCases.m'));

fprintf('\n--- Saved figures ---\n');
fprintf('  %s\n', fullfile(figFolder, 'vn_diagram.png'));
for iCase = 1:nCases
    cas = allCases(iCase);
    safeName = strrep(cas.name, ' ', '_');
    safeName = strrep(safeName, '@', 'at');
    fprintf('  %s\n', fullfile(figFolder, ['qvm_' safeName '.png']));
    fprintf('  %s\n', fullfile(figFolder, ['drag_torsion_' safeName '.png']));
end

fprintf('============================================================\n');

%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

function caseOut = solve_low_alpha_case(name, V_mph, n, W, rho, S, cbar, ...
                                        CL_alpha_3D, alpha_L0, cd_poly, CM_low, mph_to_fts)
    V = V_mph * mph_to_fts;
    CN_req = 2*n*W/(rho*V^2*S);

    fun = @(a) low_alpha_CN(a, CL_alpha_3D, alpha_L0, cd_poly) - CN_req;

    a0 = alpha_L0 + CN_req/CL_alpha_3D;

    try
        alpha = fzero(fun, a0);
    catch
        alpha = fzero(fun, [deg2rad(-10), deg2rad(10)]);
    end

    alpha_d = rad2deg(alpha);
    CL = CL_alpha_3D*(alpha - alpha_L0);
    CD = polyval(cd_poly, alpha_d);
    CM = CM_low;

    CN = CL*cos(alpha) + CD*sin(alpha);
    CX = CD*cos(alpha) - CL*sin(alpha);

    L  = 0.5*rho*V^2*S*CL;
    D  = 0.5*rho*V^2*S*CD;
    Fz = 0.5*rho*V^2*S*CN;
    Fx = 0.5*rho*V^2*S*CX;
    My = 0.5*rho*V^2*S*cbar*CM;

    caseOut.name    = name;
    caseOut.V_mph   = V_mph;
    caseOut.V       = V;
    caseOut.n       = n;
    caseOut.W       = W;
    caseOut.alpha_d = alpha_d;
    caseOut.alpha   = alpha;
    caseOut.CL      = CL;
    caseOut.CD      = CD;
    caseOut.CM      = CM;
    caseOut.CN      = CN;
    caseOut.CX      = CX;
    caseOut.L       = L;
    caseOut.D       = D;
    caseOut.Fz      = Fz;
    caseOut.Fx      = Fx;
    caseOut.My      = My;
end

function CN = low_alpha_CN(alpha, CL_alpha_3D, alpha_L0, cd_poly)
    CL = CL_alpha_3D*(alpha - alpha_L0);
    CD = polyval(cd_poly, rad2deg(alpha));
    CN = CL*cos(alpha) + CD*sin(alpha);
end

function dist = build_vertical_distribution(F_wing, L, y)
    % Averages uniform and elliptical load distributions.
    % Integrates tip-inward to obtain Vz(y) and Mx(y).
    %
    % Convention: V(y) = integral_y^L q(xi) dxi
    %   -> V(0) = F_wing  (same sign as applied load)
    %   -> V(L) = 0       (free tip boundary condition)
    %   M(y) = integral_y^L V(xi) dxi
    %   -> M(0) = root bending moment, M(L) = 0
    q_uniform = (F_wing/L) * ones(size(y));

    q0 = 4*F_wing/(pi*L);
    q_elliptic = q0 * sqrt(max(0, 1 - (y/L).^2));

    q = 0.5*(q_uniform + q_elliptic);

    % Tip-inward shear: V(y) = integral_y^L q dxi
    I_q = cumtrapz(y, q);           % I_q(y) = integral_0^y q dxi
    V   = I_q(end) - I_q;           % V(y)   = I_q(L) - I_q(y)

    % Tip-inward moment: M(y) = integral_y^L V dxi
    I_V = cumtrapz(y, V);
    M   = I_V(end) - I_V;

    dist.F_wing     = F_wing;
    dist.q_uniform  = q_uniform;
    dist.q_elliptic = q_elliptic;
    dist.q          = q;
    dist.V          = V;
    dist.M          = M;
    dist.V_root     = V(1);
    dist.M_root     = M(1);
end

function dist = build_drag_distribution(Fx_wing, L, y)
    % Builds the spanwise drag load distribution qx(y).
    % Shape: uniform along the inboard 80% of semispan, with a
    % 25% increase over the outer 20% to model tip effects.
    %
    % Convention: Vx(y) = integral_y^L qx(xi) dxi
    %   -> Vx(0) = Fx_wing  (same sign as applied load)
    %   -> Vx(L) = 0
    %   Mz(y) = integral_y^L Vx(xi) dxi
    %   -> Mz(0) = root in-plane bending moment, Mz(L) = 0
    shape = ones(size(y));
    shape(y >= 0.8*L) = 1.25;

    scale = Fx_wing / trapz(y, shape);
    qx    = scale * shape;

    I_qx = cumtrapz(y, qx);
    Vx   = I_qx(end) - I_qx;

    I_Vx = cumtrapz(y, Vx);
    Mz   = I_Vx(end) - I_Vx;

    dist.qx      = qx;
    dist.Fx_wing = trapz(y, qx);
    dist.Vx      = Vx;
    dist.Mz      = Mz;
    dist.Vx_root = Vx(1);
    dist.Mz_root = Mz(1);
end

function dist = build_uniform_distribution(Mwing, L, y)
    % Distributes the semi-wing pitching moment uniformly along the span.
    %
    % Convention: My(y) = integral_y^L qm(xi) dxi
    %   -> My(0) = Mwing  (same sign as applied pitching moment)
    %   -> My(L) = 0
    % NOTE: this represents only the aerodynamic CM contribution.
    % The torsional contribution from Fz acting at a moment arm from
    % the shear center is added in the structural analysis once the wing-box geometry
    % (and hence shear-center location) is known.
    qm = (Mwing/L) * ones(size(y));

    I_qm    = cumtrapz(y, qm);
    My_dist = I_qm(end) - I_qm;

    dist.qm      = qm;
    dist.Mwing   = trapz(y, qm);
    dist.My_dist = My_dist;
    dist.My_root = My_dist(1);
end

function fig = plot_qVM(y, dist, caseName)
    % Three-panel figure: spanwise vertical load, shear force, bending moment.
    fig = figure('Name',['Vertical q-V-M  ' caseName],'Color','w','NumberTitle','off');
    tiledlayout(3,1);

    nexttile;
    plot(y, dist.q, 'b', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('q_z(y)  [lb/ft]');
    title(['Vertical load distribution - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, dist.V, 'b', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('V_z(y)  [lb]');
    title(['Vertical shear force - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, dist.M, 'b', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('M_x(y)  [lb-ft]');
    title(['Bending moment (vertical) - ' caseName], 'Interpreter', 'none');
end

function fig = plot_drag_torsion(y, distX, distMy, caseName)
    % Five-panel figure: spanwise drag load, in-plane shear, in-plane
    % bending moment, pitching-moment distribution, and cumulative torsion.
    fig = figure('Name',['Drag/Torsion  ' caseName],'Color','w','NumberTitle','off');
    tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(y, distX.qx, 'r', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('q_x(y)  [lb/ft]');
    title(['Drag load distribution - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, distX.Vx, 'r', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('V_x(y)  [lb]');
    title(['In-plane shear force (drag) - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, distX.Mz, 'r', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('M_z(y)  [lb-ft]');
    title(['In-plane bending moment (drag) - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, distMy.qm, 'g', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('q_m(y)  [lb-ft/ft]');
    title(['Pitching-moment distribution - ' caseName], 'Interpreter', 'none');

    nexttile;
    plot(y, distMy.My_dist, 'g', 'LineWidth', 2);
    grid on; box on;
    xlabel('y from root [ft]');
    ylabel('M_y(y)  [lb-ft]');
    title(['Cumulative torsion (aero CM only) - ' caseName], 'Interpreter', 'none');
end

function print_case(c)
    fprintf('\n%s\n', c.name);
    fprintf('  V      = %8.2f mph\n',   c.V_mph);
    fprintf('  n      = %8.4f\n',       c.n);
    fprintf('  alpha  = %8.3f deg\n',   c.alpha_d);
    fprintf('  CL     = %8.4f\n',       c.CL);
    fprintf('  CD     = %8.5f\n',       c.CD);
    fprintf('  CM     = %8.4f\n',       c.CM);
    fprintf('  CN     = %8.4f\n',       c.CN);
    fprintf('  CX     = %8.4f\n',       c.CX);
    fprintf('  L      = %8.1f lb\n',    c.L);
    fprintf('  D      = %8.1f lb\n',    c.D);
    fprintf('  Fz     = %8.1f lb\n',    c.Fz);
    fprintf('  Fx     = %8.1f lb\n',    c.Fx);
    fprintf('  My     = %8.1f lb-ft\n', c.My);
end
