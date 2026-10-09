function loads = criticalLoads(loadCaseTable, spanwiseLoads)
% criticalLoads Convert saved aerodynamic load tables into structural load inputs.
%
% Inputs:
%   loadCaseTable    - Table loaded from data/critical_load_cases.csv.
%   spanwiseLoads - Optional struct loaded from data/spanwise_loads.mat.
%
% Outputs:
%   loads - Struct containing original loads and unit-converted fields for
%           cross-section structural analysis.
%
% Notes:
%   This function does not change the saved aerodynamic load methodology. It only
%   packages existing load values and converts ft-based moments to in-lb.

if nargin < 1 || isempty(loadCaseTable)
error('criticalLoads:MissingInput', 'A non-empty aerodynamic load table is required.');
end
if nargin < 2
    spanwiseLoads = [];
end

requiredColumns = ["case_name", "V_mph", "n", "W_lb", "alpha_deg", ...
    "CL", "CD", "CM", "CN", "CX", ...
    "Fz_total_lb", "Fx_total_lb", "My_total_lbft", ...
    "Fz_semispan_lb", "Fx_semispan_lb", "My_semispan_lbft", ...
    "Vz_root_lb", "Vx_root_lb", "Mx_root_lbft", "Mz_root_lbft", ...
    "My_root_lbft_or_semispan_value", "notes"];

missingColumns = setdiff(requiredColumns, string(loadCaseTable.Properties.VariableNames));
if ~isempty(missingColumns)
    error('criticalLoads:MissingColumns', ...
        'Aerodynamic load table is missing required columns: %s', strjoin(missingColumns, ', '));
end

loads.sourceTable = loadCaseTable;
loads.caseNames = loadCaseTable.case_name;
loads.nCases = height(loadCaseTable);

loads.aero.V_mph = loadCaseTable.V_mph;
loads.aero.n = loadCaseTable.n;
loads.aero.W_lb = loadCaseTable.W_lb;
loads.aero.alpha_deg = loadCaseTable.alpha_deg;
loads.aero.CL = loadCaseTable.CL;
loads.aero.CD = loadCaseTable.CD;
loads.aero.CM = loadCaseTable.CM;
loads.aero.CN = loadCaseTable.CN;
loads.aero.CX = loadCaseTable.CX;

loads.total.Fz_lb = loadCaseTable.Fz_total_lb;
loads.total.Fx_lb = loadCaseTable.Fx_total_lb;
loads.total.My_lbft = loadCaseTable.My_total_lbft;
loads.total.My_lbin = 12.0 * loadCaseTable.My_total_lbft;

loads.semispan.Fz_lb = loadCaseTable.Fz_semispan_lb;
loads.semispan.Fx_lb = loadCaseTable.Fx_semispan_lb;
loads.semispan.My_lbft = loadCaseTable.My_semispan_lbft;
loads.semispan.My_lbin = 12.0 * loadCaseTable.My_semispan_lbft;

loads.root.Vz_lb = loadCaseTable.Vz_root_lb;
loads.root.Vx_lb = loadCaseTable.Vx_root_lb;
loads.root.Mx_lbft = loadCaseTable.Mx_root_lbft;
loads.root.Mz_lbft = loadCaseTable.Mz_root_lbft;
loads.root.My_lbft = loadCaseTable.My_root_lbft_or_semispan_value;
loads.root.Mx_lbin = 12.0 * loadCaseTable.Mx_root_lbft;
loads.root.Mz_lbin = 12.0 * loadCaseTable.Mz_root_lbft;
loads.root.My_lbin = 12.0 * loadCaseTable.My_root_lbft_or_semispan_value;

loads.units.force = "lb";
loads.units.length = "in";
loads.units.moment = "lb-in";
loads.units.sourceLength = "ft";
loads.notes = loadCaseTable.notes;

loads.hasSpanwise = false;
loads.spanwise = struct();
loads.geometry = struct();

if ~isempty(spanwiseLoads)
    requiredSpanwiseFields = ["y", "caseNames", "allCases", "distZ", "distX", ...
        "distMy", "idx_pos", "idx_neg", "GOV_POS", "GOV_NEG", "geometry"];
    missingSpanwiseFields = setdiff(requiredSpanwiseFields, string(fieldnames(spanwiseLoads)));
    if ~isempty(missingSpanwiseFields)
        error('criticalLoads:MissingSpanwiseFields', ...
            'Saved spanwise load struct is missing required fields: %s', ...
            strjoin(missingSpanwiseFields, ', '));
    end

    requiredGeometryFields = ["c", "cbar", "semiSpan", "b", "S", "AR"];
    missingGeometryFields = setdiff(requiredGeometryFields, string(fieldnames(spanwiseLoads.geometry)));
    if ~isempty(missingGeometryFields)
        error('criticalLoads:MissingGeometryFields', ...
            'Saved geometry struct is missing required fields: %s', ...
            strjoin(missingGeometryFields, ', '));
    end

    loads.hasSpanwise = true;
    loads.spanwise.y_ft = spanwiseLoads.y;
    loads.spanwise.y_in = 12.0 * spanwiseLoads.y;
    loads.spanwise.caseNames = spanwiseLoads.caseNames;
    loads.spanwise.allCases = spanwiseLoads.allCases;
    loads.spanwise.distZ = spanwiseLoads.distZ;
    loads.spanwise.distX = spanwiseLoads.distX;
    loads.spanwise.distMy = spanwiseLoads.distMy;

    loads.governing.idx_pos = spanwiseLoads.idx_pos;
    loads.governing.idx_neg = spanwiseLoads.idx_neg;
    loads.governing.positive = spanwiseLoads.GOV_POS;
    loads.governing.negative = spanwiseLoads.GOV_NEG;

    loads.geometry.ft = spanwiseLoads.geometry;
    loads.geometry.in.c = 12.0 * spanwiseLoads.geometry.c;
    loads.geometry.in.cbar = 12.0 * spanwiseLoads.geometry.cbar;
    loads.geometry.in.semiSpan = 12.0 * spanwiseLoads.geometry.semiSpan;
    loads.geometry.in.b = 12.0 * spanwiseLoads.geometry.b;
    loads.geometry.in.S = 144.0 * spanwiseLoads.geometry.S;
    loads.geometry.in.AR = spanwiseLoads.geometry.AR;
end

% Structural geometry, material properties, failure checks, and load-case
% stress scaling are handled by the final structural analysis modules.
end
