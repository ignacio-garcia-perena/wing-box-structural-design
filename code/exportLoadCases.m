% exportLoadCases
% Export the critical flight load cases for the structural analysis.
%
% This script assumes initialLoadGeneration.m has already computed the
% load cases, spanwise distributions, governing cases, and basic geometry.

requiredVars = {'allCases', 'distZ', 'distX', 'distMy', 'nCases', ...
    'y', 'idx_pos', 'idx_neg', 'GOV_POS', 'GOV_NEG', ...
    'c', 'cbar', 'semiSpan', 'b', 'S', 'AR'};
for iVar = 1:numel(requiredVars)
    if ~exist(requiredVars{iVar}, 'var')
        error('exportLoadCases:MissingVariable', ...
            'Required variable "%s" was not found. Run this after the load cases are computed.', ...
            requiredVars{iVar});
    end
end

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
dataFolder = fullfile(projectRoot, 'data');
if ~exist(dataFolder, 'dir')
    mkdir(dataFolder);
end

case_name = strings(nCases, 1);
V_mph = zeros(nCases, 1);
n_col = zeros(nCases, 1);
W_lb = zeros(nCases, 1);
alpha_deg = zeros(nCases, 1);
CL = zeros(nCases, 1);
CD = zeros(nCases, 1);
CM = zeros(nCases, 1);
CN = zeros(nCases, 1);
CX = zeros(nCases, 1);
Fz_total_lb = zeros(nCases, 1);
Fx_total_lb = zeros(nCases, 1);
My_total_lbft = zeros(nCases, 1);
Fz_semispan_lb = zeros(nCases, 1);
Fx_semispan_lb = zeros(nCases, 1);
My_semispan_lbft = zeros(nCases, 1);
Vz_root_lb = zeros(nCases, 1);
Vx_root_lb = zeros(nCases, 1);
Mx_root_lbft = zeros(nCases, 1);
Mz_root_lbft = zeros(nCases, 1);
My_root_lbft_or_semispan_value = zeros(nCases, 1);
notes = strings(nCases, 1);

for iCase = 1:nCases
    cas = allCases(iCase);
    dz = distZ{iCase};
    dx = distX{iCase};
    dm = distMy{iCase};

    case_name(iCase) = string(cas.name);
    V_mph(iCase) = cas.V_mph;
    n_col(iCase) = cas.n;
    W_lb(iCase) = cas.W;
    alpha_deg(iCase) = cas.alpha_d;
    CL(iCase) = cas.CL;
    CD(iCase) = cas.CD;
    CM(iCase) = cas.CM;
    CN(iCase) = cas.CN;
    CX(iCase) = cas.CX;
    Fz_total_lb(iCase) = cas.Fz;
    Fx_total_lb(iCase) = cas.Fx;
    My_total_lbft(iCase) = cas.My;
    Fz_semispan_lb(iCase) = dz.F_wing;
    Fx_semispan_lb(iCase) = dx.Fx_wing;
    My_semispan_lbft(iCase) = dm.Mwing;
    Vz_root_lb(iCase) = dz.V_root;
    Vx_root_lb(iCase) = dx.Vx_root;
    Mx_root_lbft(iCase) = dz.M_root;
    Mz_root_lbft(iCase) = dx.Mz_root;
    My_root_lbft_or_semispan_value(iCase) = dm.My_root;
    notes(iCase) = "Aero CM torsion only; structural solver adds vertical-load shear-center moment arm.";
end

criticalLoadCases = table(case_name, V_mph, n_col, W_lb, alpha_deg, ...
    CL, CD, CM, CN, CX, ...
    Fz_total_lb, Fx_total_lb, My_total_lbft, ...
    Fz_semispan_lb, Fx_semispan_lb, My_semispan_lbft, ...
    Vz_root_lb, Vx_root_lb, Mx_root_lbft, Mz_root_lbft, ...
    My_root_lbft_or_semispan_value, notes, ...
    'VariableNames', {'case_name', 'V_mph', 'n', 'W_lb', 'alpha_deg', ...
    'CL', 'CD', 'CM', 'CN', 'CX', ...
    'Fz_total_lb', 'Fx_total_lb', 'My_total_lbft', ...
    'Fz_semispan_lb', 'Fx_semispan_lb', 'My_semispan_lbft', ...
    'Vz_root_lb', 'Vx_root_lb', 'Mx_root_lbft', 'Mz_root_lbft', ...
    'My_root_lbft_or_semispan_value', 'notes'});

outFile = fullfile(dataFolder, 'critical_load_cases.csv');
writetable(criticalLoadCases, outFile);

spanwiseLoads = struct();
spanwiseLoads.y = y;
spanwiseLoads.caseNames = case_name;
spanwiseLoads.allCases = allCases;
spanwiseLoads.distZ = distZ;
spanwiseLoads.distX = distX;
spanwiseLoads.distMy = distMy;
spanwiseLoads.idx_pos = idx_pos;
spanwiseLoads.idx_neg = idx_neg;
spanwiseLoads.GOV_POS = GOV_POS;
spanwiseLoads.GOV_NEG = GOV_NEG;
spanwiseLoads.geometry.c = c;
spanwiseLoads.geometry.cbar = cbar;
spanwiseLoads.geometry.semiSpan = semiSpan;
spanwiseLoads.geometry.b = b;
spanwiseLoads.geometry.S = S;
spanwiseLoads.geometry.AR = AR;

spanwiseOutFile = fullfile(dataFolder, 'spanwise_loads.mat');
save(spanwiseOutFile, 'spanwiseLoads');

fprintf('\n--- Load case export ---\n');
fprintf('  %s\n', outFile);
fprintf('  %s\n', spanwiseOutFile);
