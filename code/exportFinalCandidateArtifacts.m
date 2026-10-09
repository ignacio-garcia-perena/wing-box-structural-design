function exportFinalCandidateArtifacts(design, spanwiseResult, analysisOptions, optimisationResults)
% exportFinalCandidateArtifacts Regenerate active final report artifacts.
%
% Inputs:
%   design              - Active final design struct.
%   spanwiseResult      - Detailed result from analyzeWingDesignSpanwise.
%   analysisOptions     - Final analysis options.
%   optimisationResults - Optional runFinalMonteCarlo results struct.

if nargin < 4
    optimisationResults = struct();
end
validateInputs(design, spanwiseResult, analysisOptions);

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
tablesFolder = fullfile(projectRoot, 'outputs', 'tables');
figuresFolder = fullfile(projectRoot, 'outputs', 'figures');
notesFolder = fullfile(projectRoot, 'docs', 'project_notes');
if ~isfolder(tablesFolder), mkdir(tablesFolder); end
if ~isfolder(figuresFolder), mkdir(figuresFolder); end
if ~isfolder(notesFolder), mkdir(notesFolder); end

writeCandidateSummary(design, spanwiseResult, analysisOptions, tablesFolder);
writeDesignGeometry(design, spanwiseResult, analysisOptions, tablesFolder);
writeMaterials(spanwiseResult, tablesFolder);
writeWeightSummary(spanwiseResult, tablesFolder);
writeFailureSummary(spanwiseResult, tablesFolder);
writePanelBucklingSummary(spanwiseResult, analysisOptions, tablesFolder);
writeSpanwiseSummary(spanwiseResult, tablesFolder);
writeSpanwisePanelBucklingByBay(spanwiseResult, tablesFolder);
writeLimitationsSummary(tablesFolder);
writeCandidateStatusNote(design, spanwiseResult, analysisOptions, optimisationResults, notesFolder);
writeCripplingStatusNote(spanwiseResult, notesFolder);

makeCrossSectionFigure(spanwiseResult, figuresFolder);
makeSpanwiseYieldFigure(spanwiseResult, figuresFolder);
makeSpanwisePanelCombinedFigure(spanwiseResult, figuresFolder);
makePanelBucklingRatiosFigure(spanwiseResult, figuresFolder);
makeWeightBreakdownFigure(spanwiseResult, figuresFolder);
end

function validateInputs(design, spanwiseResult, analysisOptions)
if nargin < 1 || isempty(design) || ~isstruct(design)
    error('exportFinalCandidateArtifacts:MissingDesign', 'design struct is required.');
end
if nargin < 2 || isempty(spanwiseResult) || ~isstruct(spanwiseResult) || ...
        ~isfield(spanwiseResult, 'summary')
    error('exportFinalCandidateArtifacts:MissingResult', ...
        'spanwiseResult from analyzeWingDesignSpanwise is required.');
end
if nargin < 3 || isempty(analysisOptions) || ~isstruct(analysisOptions)
    error('exportFinalCandidateArtifacts:MissingOptions', ...
        'analysisOptions struct is required.');
end
end

function writeCandidateSummary(design, result, analysisOptions, tablesFolder)
s = result.summary;
ribs = result.weight.components.ribs;
pb = s.panelBuckling;
combo = pb.maxCombinedActive;
row = table( ...
    "final_optimised_candidate", ...
    design.box.frontSpar_xc, design.box.rearSpar_xc, ...
    string(design.stringers.section_key), design.stringers.sec_no, ...
    design.sparCaps.sec_no, ...
    design.thickness.t3_upperMainSkin_in, design.thickness.t5_noseSkin_in, ...
    design.thickness.t1_frontSpar_in, design.thickness.t2_rearSpar_in, ...
    design.stringers.numTop, design.stringers.numBottom, ...
    design.stringers.numNoseTop, design.stringers.numNoseBottom, ...
    design.ribs.ribSpacing_in, design.columnBuckling.K_column, ...
    analysisOptions.loadFactor, logical(analysisOptions.usePanelIdealisation), ...
    logical(analysisOptions.enablePanelBuckling), ...
    logical(analysisOptions.panelBuckling.enableCombinedRatios), ...
    analysisOptions.panelBuckling.Kc, analysisOptions.panelBuckling.Kb, ...
    analysisOptions.panelBuckling.Ks, string(design.ribs.weightModel), ...
    design.ribs.weightPerRib_lb, string(ribs.material.name) + "-" + ...
    string(ribs.material.temper) + " " + string(ribs.material.product_form), ...
    ribs.grossArea_in2, ribs.solidAreaFactor, ribs.thickness_in, ...
    ribs.rib_weight_single_lb, ribs.count_fullspan, ...
    s.weight_excluding_ribs_misc_lb, ribs.fullspan_weight_lb, ...
    s.weight_including_ribs_excluding_misc_lb, s.max_yield_ratio, ...
    s.max_column_buckling_ratio, s.divergence_ratio, pb.max_Rc, pb.max_Rb, ...
    pb.max_Rs, pb.max_R_SC, pb.max_R_BS, pb.max_combined_active, ...
    getFieldOrNaN(combo, 'bay_index'), getFieldOrBlank(combo, 'case_name'), ...
    getFieldOrBlank(combo, 'panel_name'), getFieldOrBlank(combo, 'combined_mode'), ...
    logical(s.overall_pass), ...
    "Active final candidate selected as the lightest passing design from final Monte Carlo and regenerated with analyzeWingDesignSpanwise.", ...
    'VariableNames', { ...
    'design_label','frontSpar_xc','rearSpar_xc','stringer_section_key', ...
    'stringer_sec_no','spar_cap_sec_no','main_skin_t','noseSkin_t', ...
    't1_frontSpar_in','t2_rearSpar_in','numTop','numBottom', ...
    'numNoseTop','numNoseBottom','ribSpacing_in','K_column','loadFactor', ...
    'usePanelIdealisation','enablePanelBuckling','enableCombinedPanelRatios', ...
    'Kc','Kb','Ks','rib_weight_model','weightPerRib_lb','rib_material', ...
    'rib_grossArea_in2','rib_solidAreaFactor','rib_thickness_in', ...
    'rib_weight_single_lb','rib_count_fullspan','weight_excluding_ribs_misc_lb', ...
    'rib_weight_fullspan_lb','weight_including_ribs_excluding_misc_lb', ...
    'max_yield_ratio','max_column_buckling_ratio','divergence_ratio', ...
    'max_Rc','max_Rb','max_Rs','max_R_SC','max_R_BS', ...
    'max_active_combined_panel_ratio','controlling_combined_bay', ...
    'controlling_combined_case','controlling_combined_panel', ...
    'controlling_combined_mode','overall_pass','notes'});
writetable(row, fullfile(tablesFolder, 'final_candidate_summary.csv'));
end

function writeDesignGeometry(design, result, analysisOptions, tablesFolder)
wing = result.wing;
section = result.section;
ribs = result.weight.components.ribs;
rows = strings(0, 4);
rows = addRow(rows, "frontSpar_xc", wing.box.frontSpar_xc, "c", "");
rows = addRow(rows, "rearSpar_xc", wing.box.rearSpar_xc, "c", "");
rows = addRow(rows, "frontSpar_x_in", wing.box.frontSpar_x_in, "in", "");
rows = addRow(rows, "rearSpar_x_in", wing.box.rearSpar_x_in, "in", "");
rows = addRow(rows, "chord_in", wing.planform.chord_in, "in", "");
rows = addRow(rows, "semispan_in", wing.planform.semispan_in, "in", "");
rows = addRow(rows, "upperMainSkin_t_in", wing.thickness.t3_upperMainSkin_in, "in", "");
rows = addRow(rows, "lowerMainSkin_t_in", wing.thickness.t4_lowerMainSkin_in, "in", "");
rows = addRow(rows, "noseSkin_t_in", wing.thickness.t5_noseSkin_in, "in", "");
rows = addRow(rows, "upperNoseSkin_t_in", wing.thickness.upperNoseSkin_in, "in", "");
rows = addRow(rows, "lowerNoseSkin_t_in", wing.thickness.lowerNoseSkin_in, "in", "");
rows = addRow(rows, "frontSparWeb_t_in", wing.thickness.t1_frontSpar_in, "in", "");
rows = addRow(rows, "rearSparWeb_t_in", wing.thickness.t2_rearSpar_in, "in", "");
rows = addRow(rows, "stringer_section_key", design.stringers.section_key, "", "");
rows = addRow(rows, "stringer_sec_no", design.stringers.sec_no, "", "");
rows = addRow(rows, "stringer_area_in2", wing.stringers.areaTop_in2, "in^2", "");
rows = addRow(rows, "numTopStringers", wing.stringers.numTop, "count", "");
rows = addRow(rows, "numBottomStringers", wing.stringers.numBottom, "count", "");
rows = addRow(rows, "numNoseTopStringers", wing.stringers.numNoseTop, "count", "");
rows = addRow(rows, "numNoseBottomStringers", wing.stringers.numNoseBottom, "count", "");
rows = addRow(rows, "spar_cap_sec_no", wing.sparCaps.section.sec_no, "", "");
rows = addRow(rows, "spar_cap_area_in2", wing.sparCaps.section.area_in2, "in^2", "");
rows = addRow(rows, "ribSpacing_in", wing.ribs.spacing_in, "in", "");
rows = addRow(rows, "rib_count_semispan", ribs.count_semispan, "count", "");
rows = addRow(rows, "rib_count_fullspan", ribs.count_fullspan, "count", "");
rows = addRow(rows, "rib_weight_model", ribs.weightModel, "", "Geometric rib weight model uses NACA 2412 enclosed area and 2024-T3 sheet density.");
rows = addRow(rows, "rib_A_over_c2", wing.ribs.normalizedArea_A_over_c2, "", "");
rows = addRow(rows, "rib_chord_in", wing.ribs.chord_in, "in", "");
rows = addRow(rows, "rib_grossArea_in2", ribs.grossArea_in2, "in^2", "");
rows = addRow(rows, "rib_solidAreaFactor", ribs.solidAreaFactor, "", "");
rows = addRow(rows, "rib_thickness_in", ribs.thickness_in, "in", "");
rows = addRow(rows, "rib_weight_single_lb", ribs.rib_weight_single_lb, "lb", "");
rows = addRow(rows, "rib_weight_fullspan_lb", ribs.fullspan_weight_lb, "lb", "");
rows = addRow(rows, "panelIdealisationMethod", displayMethodName(analysisOptions.panelIdealisationMethod), "", "Modified Megson panel idealisation; equal split retained only as regression option.");
rows = addRow(rows, "panelIdealisationCentroidReference", section.idealization.centroidReference, "", "Boom-only centroid zc used for modified Megson z-ratio area assignment.");
rows = addRow(rows, "panelIdealisationFallbackPanels", section.idealization.nFallbackPanels, "count", "Fallback equal split used for centroid-straddling or virtual-endpoint panels.");
T = table(rows(:,1), rows(:,2), rows(:,3), rows(:,4), ...
    'VariableNames', {'item','value','units','notes'});
writetable(T, fullfile(tablesFolder, 'final_design_geometry.csv'));
end

function writeMaterials(result, tablesFolder)
skin = result.wing.material.skinWeb;
cap = result.wing.material.sparCapStringer;
T = table( ...
    ["skins_webs_ribs"; "stringers_spar_caps"], ...
    [skin.name; cap.name], [skin.temper; cap.temper], ...
    [skin.product_form; cap.product_form], ...
    [skin.E_psi; cap.E_psi], [skin.Ec_psi; cap.Ec_psi], ...
    [skin.G_psi; cap.G_psi], [skin.rho_lbin3; cap.rho_lbin3], ...
    [skin.Ftu_psi; cap.Ftu_psi], [skin.Fty_psi; cap.Fty_psi], ...
    [skin.Fcy_psi; cap.Fcy_psi], [skin.Fsu_psi; cap.Fsu_psi], ...
    [skin.source; cap.source], [skin.confidence; cap.confidence], ...
    'VariableNames', {'component','alloy','temper','product_form','E_psi', ...
    'Ec_psi','G_psi','density_lb_in3','Ftu_psi','Fty_psi','Fcy_psi', ...
    'Fsu_psi','source_reference','confidence'});
writetable(T, fullfile(tablesFolder, 'final_materials.csv'));
end

function writeWeightSummary(result, tablesFolder)
w = result.weight;
T = table( ...
    ["skins_webs"; "spar_caps"; "stringers"; "ribs"; ...
    "total_excluding_ribs_misc"; "total_including_ribs_excluding_misc"; ...
    "excluded_items"], ...
    [w.semispan.panel_lb; w.semispan.sparCap_lb; w.semispan.stringer_lb; ...
    w.semispan.rib_lb; w.semispan.total_excluding_ribs_and_misc_lb; ...
    w.semispan.total_including_ribs_excluding_misc_lb; NaN], ...
    [2*w.semispan.panel_lb; 2*w.semispan.sparCap_lb; 2*w.semispan.stringer_lb; ...
    w.fullspan.rib_lb; w.fullspan.total_excluding_ribs_and_misc_lb; ...
    w.fullspan.total_including_ribs_excluding_misc_lb; NaN], ...
    ["physical skin and web panels"; "A3.13 T-section spar caps"; ...
    "A3.10 Z-section stringers"; "geometric NACA 2412 rib weight model"; ...
    "skins/webs plus spar caps plus stringers"; ...
    "reported final analysis weight"; ...
    "fasteners joints doublers paint fuel access panels manufacturing allowances and miscellaneous items excluded"], ...
    'VariableNames', {'component','semispan_lb','fullspan_lb','notes'});
writetable(T, fullfile(tablesFolder, 'final_weight_summary.csv'));
end

function writeFailureSummary(result, tablesFolder)
s = result.summary;
pb = s.panelBuckling;
rows = strings(0, 13);
rows = addMetric(rows, "yield_ratio", s.max_yield_ratio, "<= 1", s.yield_pass, s.maxYield, "", "");
rows = addMetric(rows, "column_buckling_ratio", s.max_column_buckling_ratio, "<= 1", s.column_buckling_pass, s.maxColumnBuckling, "", "");
rows = addMetric(rows, "divergence_ratio", s.divergence_ratio, "> 1", s.divergence_pass, struct(), "", "");
rows = addMetric(rows, "panel_Rc", pb.max_Rc, "<= 1", pb.pass, pb.maxRc, "Rc", "");
rows = addMetric(rows, "panel_Rb", pb.max_Rb, "<= 1", pb.pass, pb.maxRb, "Rb", "");
rows = addMetric(rows, "panel_Rs", pb.max_Rs, "<= 1", pb.pass, pb.maxRs, "Rs", "");
rows = addMetric(rows, "panel_R_SC_global", pb.max_R_SC, "diagnostic only", NaN, pb.maxR_SC, "R_SC", "global diagnostic");
rows = addMetric(rows, "panel_R_BS_global", pb.max_R_BS, "diagnostic only", NaN, pb.maxR_BS, "R_BS", "global diagnostic");
rows = addMetric(rows, "panel_active_combined_ratio", pb.max_combined_active, "<= 1", pb.combined_pass, pb.maxCombinedActive, "active", getFieldOrBlank(pb.maxCombinedActive, 'combined_mode'));
rows = addMetric(rows, "overall_final_analysis_status", double(s.overall_pass), "true", s.overall_pass, struct(), "", "");
T = table(rows(:,1), rows(:,2), rows(:,3), rows(:,4), rows(:,5), ...
    rows(:,6), rows(:,7), rows(:,8), rows(:,9), rows(:,10), rows(:,11), ...
    rows(:,12), rows(:,13), ...
    'VariableNames', {'metric','value','passCriteria','pass','controlStation', ...
    'controlBay','controlBayStart_in','controlBayEnd_in','controlCase', ...
    'controlItem','controlMode','recoveryMethod','idealisationMethod'});
writetable(T, fullfile(tablesFolder, 'final_failure_summary.csv'));
end

function writePanelBucklingSummary(result, analysisOptions, tablesFolder)
pb = result.summary.panelBuckling;
cc = aggregateCurvedCorrection(result.bayResults);
rows = strings(0, 11);
method = "bay_averaged_panel_buckling_recovery_with_curvature_correction";
ideal = displayMethodName(analysisOptions.panelIdealisationMethod);
rows = addPB(rows, "Kc", analysisOptions.panelBuckling.Kc, "compression coefficient", struct(), "", method, ideal);
rows = addPB(rows, "Kb", analysisOptions.panelBuckling.Kb, "bending coefficient", struct(), "", method, ideal);
rows = addPB(rows, "Ks", analysisOptions.panelBuckling.Ks, "shear coefficient", struct(), "", method, ideal);
rows = addPB(rows, "curved_correction_enabled", string(analysisOptions.panelBuckling.enableCurvedPanelCorrection), "cylindrical curvature compression correction enabled", struct(), "", method, ideal);
rows = addPB(rows, "curved_correction_method", displayMethodName(analysisOptions.panelBuckling.curvedPanelCorrectionMethod), "sigma_cr_c_effective = sigma_cr_c_flat + sigma_cr_cylinder", struct(), "", method, ideal);
rows = addPB(rows, "curved_correction_apply_to", analysisOptions.panelBuckling.curvedPanelApplyTo, "curved correction applies to compression Rc only", struct(), "", method, ideal);
rows = addPB(rows, "curved_valid_range_applications", cc.validRange, "500 < R/t < 3000 and 0.1 < L/R < 2.5", struct(), "", method, ideal);
rows = addPB(rows, "curved_high_R_over_t_applications", cc.highR, "R/t >= 3000", struct(), "", method, ideal);
rows = addPB(rows, "curved_flat_only_outside_validity", cc.flatOnly, "flat compression critical stress retained outside curved-panel validity range", struct(), "", method, ideal);
rows = addPB(rows, "curved_not_curved_panel_count", cc.notCurved, "spar webs/non-airfoil panels", struct(), "", method, ideal);
rows = addPB(rows, "curved_invalid_geometry_fallback_count", cc.invalidGeometry, "curved panels with invalid radius/slope geometry", struct(), "", method, ideal);
rows = addPB(rows, "curved_panel_instances", cc.curvedInstances, "curved panel instances over all bays", struct(), "", method, ideal);
rows = addPB(rows, "curved_applied_instances", cc.applied, "curved correction applied instances over all bays", struct(), "", method, ideal);
rows = addPB(rows, "max_Rc", pb.max_Rc, "Rc = sigma_comp_avg / sigma_cr_c_effective", pb.maxRc, "", method, ideal);
rows = addPB(rows, "max_Rb", pb.max_Rb, "Rb = sigma_comp_avg / sigma_cr_b", pb.maxRb, "", method, ideal);
rows = addPB(rows, "max_Rs", pb.max_Rs, "Rs = abs(q_avg/t) / tau_cr_s", pb.maxRs, "", method, ideal);
rows = addPB(rows, "max_R_SC", pb.max_R_SC, "R_SC = Rs^1.5 + Rc", pb.maxR_SC, "global diagnostic", method, ideal);
rows = addPB(rows, "max_R_BS", pb.max_R_BS, "R_BS = Rb^2 + Rs^2", pb.maxR_BS, "global diagnostic", method, ideal);
rows = addPB(rows, "max_active_combined_ratio", pb.max_combined_active, "mapped active S+C or B+S combined panel ratio", pb.maxCombinedActive, getFieldOrBlank(pb.maxCombinedActive, 'combined_mode'), method, ideal);
T = table(rows(:,1), rows(:,2), rows(:,3), rows(:,4), rows(:,5), ...
    rows(:,6), rows(:,7), rows(:,8), rows(:,9), rows(:,10), rows(:,11), ...
    'VariableNames', {'metric','value','equation_or_definition', ...
    'controlling_bay','bay_start_in','bay_end_in','controlling_case', ...
    'controlling_panel','mode','recovery_method','panel_idealisation_method'});
writetable(T, fullfile(tablesFolder, 'final_panel_buckling_summary.csv'));
end

function writeSpanwiseSummary(result, tablesFolder)
n = numel(result.stationResults);
station_index = zeros(n,1); y_in = zeros(n,1); y_ft = zeros(n,1);
max_yield_ratio = zeros(n,1); yield_case = strings(n,1);
max_column_ratio = zeros(n,1); column_case = strings(n,1);
for i = 1:n
    st = result.stationResults(i);
    station_index(i) = st.station_index;
    y_in(i) = st.y_in;
    y_ft(i) = st.y_ft;
    max_yield_ratio(i) = st.yieldScreening.globalMaxRatio.yield_ratio;
    yield_case(i) = string(st.yieldScreening.globalMaxRatio.case_name);
    max_column_ratio(i) = st.columnBuckling.globalMaxRatio.buckling_ratio;
    column_case(i) = string(st.columnBuckling.globalMaxRatio.case_name);
end
T = table(station_index, y_in, y_ft, max_yield_ratio, yield_case, ...
    max_column_ratio, column_case);
writetable(T, fullfile(tablesFolder, 'final_spanwise_summary.csv'));
end

function writeSpanwisePanelBucklingByBay(result, tablesFolder)
n = numel(result.bayResults);
bay_index = zeros(n,1); y_start_in = zeros(n,1); y_end_in = zeros(n,1);
y_mid_in = zeros(n,1); max_active_combined = zeros(n,1);
control_case = strings(n,1); control_panel = strings(n,1); control_mode = strings(n,1);
max_Rc = zeros(n,1); max_Rb = zeros(n,1); max_Rs = zeros(n,1);
for i = 1:n
    bay = result.bayResults(i);
    pb = bay.panelBuckling;
    bay_index(i) = bay.bay_index;
    y_start_in(i) = bay.y_start_in;
    y_end_in(i) = bay.y_end_in;
    y_mid_in(i) = bay.y_mid_in;
    max_active_combined(i) = pb.maxCombinedActive.ratio;
    control_case(i) = string(pb.maxCombinedActive.case_name);
    control_panel(i) = string(pb.maxCombinedActive.panel_name);
    control_mode(i) = string(pb.maxCombinedActive.combined_mode);
    max_Rc(i) = pb.maxRc.ratio;
    max_Rb(i) = pb.maxRb.ratio;
    max_Rs(i) = pb.maxRs.ratio;
end
T = table(bay_index, y_start_in, y_end_in, y_mid_in, max_active_combined, ...
    control_case, control_panel, control_mode, max_Rc, max_Rb, max_Rs);
writetable(T, fullfile(tablesFolder, 'final_spanwise_panel_buckling_by_bay.csv'));
end

function writeLimitationsSummary(tablesFolder)
T = table( ...
    ["crippling"; "detailed_rib_stress_buckling"; "fasteners_joints_doublers_misc"; ...
    "panel_coefficients"; "postbuckling_tension_field"; ...
    "manufacturing_detail_sizing"; "certification_final_margins"], ...
    ["not implemented"; "not implemented"; "excluded from weight and sizing"; ...
    "estimated high-aspect-ratio constants"; "not implemented"; ...
    "not implemented"; "not implemented"], ...
    ["Crippling is outside the project scope and excluded from numerical pass/fail."; ...
    "Ribs are included in weight and bay spacing but detailed rib stress and buckling are not checked."; ...
    "Fasteners joints doublers paint fuel access panels manufacturing allowances and miscellaneous items are excluded."; ...
    "Kc Kb and Ks are constant project screening values rather than detailed boundary-condition coefficients."; ...
    "No postbuckling or tension-field behavior is credited."; ...
    "Manufacturing detail sizing is outside the final analysis model."; ...
    "Current status is final numerical-check screening, not certification substantiation."], ...
    'VariableNames', {'limitation','status','detail'});
writetable(T, fullfile(tablesFolder, 'final_limitations_summary.csv'));
end

function writeCandidateStatusNote(design, result, analysisOptions, optimisationResults, notesFolder)
s = result.summary;
pb = s.panelBuckling;
combo = pb.maxCombinedActive;
if isfield(optimisationResults, 'summary') && ~isempty(optimisationResults.summary)
    nDesigns = optimisationResults.summary.n_designs(1);
    nPassing = optimisationResults.summary.n_passing(1);
    passRate = optimisationResults.summary.pass_rate(1);
    rngSeed = optimisationResults.summary.rngSeed(1);
else
    nDesigns = NaN; nPassing = NaN; passRate = NaN; rngSeed = NaN;
end
text = [
    "# Final Candidate Status"
    ""
    "## Current Active Candidate"
    ""
    "The active final candidate was selected as the lightest passing design from the final Monte Carlo workflow and re-analysed with `analyzeWingDesignSpanwise`."
    ""
    "| Item | Value |"
    "|---|---:|"
    "| Candidate label | `final_optimised_candidate` |"
    "| Front spar x/c | `" + fmt(design.box.frontSpar_xc) + "` |"
    "| Rear spar x/c | `" + fmt(design.box.rearSpar_xc) + "` |"
    "| Stringer section | `" + string(design.stringers.section_key) + "` |"
    "| Spar cap section | `" + string(design.sparCaps.sec_no) + "` |"
    "| Main skin thickness | `" + fmt(design.thickness.t3_upperMainSkin_in) + " in` |"
    "| Nose skin thickness | `" + fmt(design.thickness.t5_noseSkin_in) + " in` |"
    "| Front spar web thickness | `" + fmt(design.thickness.t1_frontSpar_in) + " in` |"
    "| Rear spar web thickness | `" + fmt(design.thickness.t2_rearSpar_in) + " in` |"
    "| Top / bottom stringers | `" + string(design.stringers.numTop) + " / " + string(design.stringers.numBottom) + "` |"
    "| Nose top / bottom stringers | `" + string(design.stringers.numNoseTop) + " / " + string(design.stringers.numNoseBottom) + "` |"
    "| Rib spacing | `" + fmt(design.ribs.ribSpacing_in) + " in` |"
    "| Column `K_column` | `" + fmt(design.columnBuckling.K_column) + "` |"
    ""
    "## Monte Carlo Selection"
    ""
    "| Quantity | Value |"
    "|---|---:|"
    "| Designs evaluated | `" + string(nDesigns) + "` |"
    "| Passing designs | `" + string(nPassing) + "` |"
    "| Pass rate | `" + fmt(100*passRate) + "%` |"
    "| Stored RNG seed | `" + string(rngSeed) + "` |"
    ""
    "## Final Numerical Checks"
    ""
    "| Check or model feature | Status |"
    "|---|---|"
    "| Ultimate load factor | `" + fmt(analysisOptions.loadFactor) + "` applied to saved aerodynamic loads |"
    "| Spanwise/rib-station evaluation | Implemented |"
    "| Panel idealisation into boom areas | `modified_megson_panel_idealisation` |"
    "| von Mises/yield screening | Implemented |"
    "| Euler column buckling | Implemented with `K_column = 0.7` |"
    "| Divergence | Implemented |"
    "| Bay-averaged panel buckling recovery | Implemented and enabled |"
    "| Combined panel buckling ratios | Implemented and enabled |"
    "| Cylindrical curvature correction | Implemented and enabled |"
    "| Geometric rib weight | Implemented |"
    ""
    "Combined panel buckling screening uses `R_SC = Rs^1.5 + Rc < 1` and `R_BS = Rb^2 + Rs^2 < 1`."
    ""
    "The saved `My` input contains aerodynamic pitching-moment (`CM`) torsion only. The structural recovery includes vertical-load shear-center torsion through `qtot_Z`, where `torque_Z_per_lb = x_sc - 0.25c` and `computeShearFlows` uses `q = Vx*qtot_X + Vz*qtot_Z + My*qtot_Y`. Adding `Fz*(x_sc - x_ac)` directly to `My` would double-count this effect unless `qtot_Z` were redefined."
    ""
    "## Candidate Results"
    ""
    "| Result | Value |"
    "|---|---:|"
    "| Weight including ribs/excluding misc | `" + fmt(s.weight_including_ribs_excluding_misc_lb) + " lb` |"
    "| Yield ratio | `" + fmt(s.max_yield_ratio) + "` |"
    "| Column buckling ratio | `" + fmt(s.max_column_buckling_ratio) + "` |"
    "| Divergence ratio | `" + fmt(s.divergence_ratio) + "` |"
    "| Max `Rc` | `" + fmt(pb.max_Rc) + "` |"
    "| Max `Rb` | `" + fmt(pb.max_Rb) + "` |"
    "| Max `Rs` | `" + fmt(pb.max_Rs) + "` |"
    "| Max `R_SC` | `" + fmt(pb.max_R_SC) + "` |"
    "| Max `R_BS` | `" + fmt(pb.max_R_BS) + "` |"
    "| Active combined ratio | `" + fmt(pb.max_combined_active) + "` |"
    "| Controlling bay/case/panel/mode | `" + string(getFieldOrNaN(combo, 'bay_index')) + " / " + string(getFieldOrBlank(combo, 'case_name')) + " / " + string(getFieldOrBlank(combo, 'panel_name')) + " / " + string(getFieldOrBlank(combo, 'combined_mode')) + "` |"
    "| Overall final numerical-check pass | `" + string(s.overall_pass) + "` |"
    ""
    "## Checks Not Included"
    ""
    "- Crippling is outside numerical pass/fail by project scope."
    "- Detailed rib stress and rib buckling are not included."
    "- Fasteners, joints, doublers, paint, access panels, manufacturing allowances, and miscellaneous weight are excluded."
    "- Final certification margins are not implemented."
    ""
    "No formulas, material data, section data, saved aerodynamic loads, or validated analysis methods were changed when exporting this candidate."
    ];
writelines(text, fullfile(notesFolder, 'final_candidate_status.md'));
end

function writeCripplingStatusNote(result, notesFolder)
s = result.summary; pb = s.panelBuckling;
text = [
    "# Crippling Status"
    ""
    "Crippling is outside numerical pass/fail by project scope. The final optimised candidate is evaluated with yield/von Mises, Euler column buckling, divergence, spanwise/rib-station evaluation, modified Megson panel idealisation, bay-averaged panel buckling, combined panel buckling ratios, cylindrical curvature correction, and geometric rib weight."
    ""
    "| Metric | Value |"
    "|---|---:|"
    "| Weight including ribs/excluding misc | `" + fmt(s.weight_including_ribs_excluding_misc_lb) + " lb` |"
    "| Yield ratio | `" + fmt(s.max_yield_ratio) + "` |"
    "| Column buckling ratio | `" + fmt(s.max_column_buckling_ratio) + "` |"
    "| Divergence ratio | `" + fmt(s.divergence_ratio) + "` |"
    "| Max `Rc` | `" + fmt(pb.max_Rc) + "` |"
    "| Max `Rb` | `" + fmt(pb.max_Rb) + "` |"
    "| Max `Rs` | `" + fmt(pb.max_Rs) + "` |"
    "| Max active combined panel ratio | `" + fmt(pb.max_combined_active) + "` |"
    "| Overall final numerical-check pass | `" + string(s.overall_pass) + "` |"
    ""
    "Crippling applies to discrete compression members such as Z-stringers and T-spar caps, but no numerical crippling check is included because it is outside the project scope. Detailed rib stress and rib buckling are also not included."
    ];
writelines(text, fullfile(notesFolder, 'crippling_status.md'));
end

function makeCrossSectionFigure(result, figuresFolder)
fig = figure('Visible','off','Color','w');
hold on; axis equal; grid on;
panels = result.section.unitShear.stage3C.unitTotals.physicalPanels;
legendHandles = gobjects(0);
legendLabels = strings(0, 1);
skinHandle = gobjects(1);
frontSparHandle = gobjects(1);
rearSparHandle = gobjects(1);
for i = 1:numel(panels)
    hPanel = plot([panels(i).xStart_in, panels(i).xEnd_in], ...
        [panels(i).zStart_in, panels(i).zEnd_in], 'k-', 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    panelName = string(panels(i).name);
    panelGroup = string(panels(i).group);
    if isemptyHandle(frontSparHandle) && contains(panelName, "frontSpar")
        frontSparHandle = hPanel;
    elseif isemptyHandle(rearSparHandle) && contains(panelName, "rearSpar")
        rearSparHandle = hPanel;
    elseif isemptyHandle(skinHandle) && ...
            (contains(panelGroup, "skin") || contains(panelName, "top") || ...
            contains(panelName, "bottom") || contains(panelName, "Nose"))
        skinHandle = hPanel;
    end
end
if ~isemptyHandle(frontSparHandle)
    [legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, frontSparHandle, "Front spar");
end
if ~isemptyHandle(rearSparHandle)
    [legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, rearSparHandle, "Rear spar");
end
if ~isemptyHandle(skinHandle)
    [legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, skinHandle, "Skin/panels");
end
hSparCaps = plotBooms(result.section.sparCaps, 'rs');
hStringers = plotBooms(result.section.stringers.top, 'bo');
plotBooms(result.section.stringers.bottom, 'bo');
plotBooms(result.section.stringers.noseTop, 'bo');
plotBooms(result.section.stringers.noseBottom, 'bo');
if ~isemptyHandle(hStringers)
    [legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, hStringers, "Stringers");
end
if ~isemptyHandle(hSparCaps)
    [legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, hSparCaps, "Spar caps");
end
hCentroid = plot(result.section.centroid.x_in, result.section.centroid.z_in, ...
    'md', 'MarkerFaceColor','m', 'HandleVisibility','off');
hShearCenter = plot(result.section.shearCenter.x_sc_in, ...
    result.section.shearCenter.z_sc_in, 'kp', 'MarkerFaceColor','y', ...
    'HandleVisibility','off');
[legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, hCentroid, "Centroid");
[legendHandles, legendLabels] = addLegendItem(legendHandles, legendLabels, hShearCenter, "Shear center");
xlabel('x (in)'); ylabel('z (in)');
title('Final Wing-Box Cross Section');
legend(legendHandles, legendLabels, 'Location','bestoutside');
exportgraphics(fig, fullfile(figuresFolder, 'final_wingbox_cross_section.png'), 'Resolution', 200);
close(fig);
end

function h = plotBooms(booms, style)
h = gobjects(1);
if isempty(booms), return; end
x = [booms.x_in]; z = [booms.z_in];
h = plot(x, z, style, 'MarkerFaceColor', style(1), 'HandleVisibility','off');
end

function tf = isemptyHandle(h)
tf = isempty(h) || ~isgraphics(h);
end

function [handles, labels] = addLegendItem(handles, labels, h, label)
if isemptyHandle(h)
    return;
end
handles(end+1, 1) = h;
labels(end+1, 1) = string(label);
end

function makeSpanwiseYieldFigure(result, figuresFolder)
n = numel(result.stationResults);
y = zeros(n,1); r = zeros(n,1);
for i = 1:n
    y(i) = result.stationResults(i).y_in;
    r(i) = result.stationResults(i).yieldScreening.globalMaxRatio.yield_ratio;
end
fig = figure('Visible','off','Color','w');
plot(y, r, 'o-', 'LineWidth', 1.5); grid on; yline(1,'r--','Limit');
xlabel('Spanwise station y (in)'); ylabel('Max yield ratio');
title('Final Spanwise Yield Ratio');
exportgraphics(fig, fullfile(figuresFolder, 'final_spanwise_yield_ratio.png'), 'Resolution', 200);
close(fig);
end

function makeSpanwisePanelCombinedFigure(result, figuresFolder)
n = numel(result.bayResults);
y = zeros(n,1); r = zeros(n,1);
for i = 1:n
    y(i) = result.bayResults(i).y_mid_in;
    r(i) = result.bayResults(i).panelBuckling.maxCombinedActive.ratio;
end
fig = figure('Visible','off','Color','w');
plot(y, r, 'o-', 'LineWidth', 1.5); grid on; yline(1,'r--','Limit');
xlabel('Bay midpoint y (in)'); ylabel('Max active combined panel ratio');
title('Final Spanwise Combined Panel Ratio');
exportgraphics(fig, fullfile(figuresFolder, 'final_spanwise_panel_combined_ratio.png'), 'Resolution', 200);
close(fig);
end

function makePanelBucklingRatiosFigure(result, figuresFolder)
pb = result.summary.panelBuckling;
values = [pb.max_Rc, pb.max_Rb, pb.max_Rs, pb.max_combined_active];
fig = figure('Visible','off','Color','w');
bar(values); grid on; yline(1,'r--','Limit');
set(gca, 'XTickLabel', {'Rc','Rb','Rs','Active combined'});
ylabel('Ratio'); title('Final Panel Buckling Ratios');
exportgraphics(fig, fullfile(figuresFolder, 'final_panel_buckling_ratios.png'), 'Resolution', 200);
close(fig);
end

function makeWeightBreakdownFigure(result, figuresFolder)
w = result.weight;
values = [2*w.semispan.panel_lb, 2*w.semispan.sparCap_lb, ...
    2*w.semispan.stringer_lb, w.fullspan.rib_lb];
fig = figure('Visible','off','Color','w');
bar(values); grid on;
set(gca, 'XTickLabel', {'skins/webs','spar caps','stringers','ribs'});
ylabel('Weight (lb)'); title('Final Weight Breakdown');
exportgraphics(fig, fullfile(figuresFolder, 'final_weight_breakdown.png'), 'Resolution', 200);
close(fig);
end

function rows = addRow(rows, item, value, units, notes)
rows(end+1, :) = [string(item), valueToString(value), string(units), string(notes)];
end

function rows = addMetric(rows, metric, value, criteria, pass, entry, ratioName, mode)
rows(end+1, :) = [string(metric), valueToString(value), string(criteria), ...
    valueToString(pass), valueToString(getFieldOrNaN(entry, 'station_index')), ...
    valueToString(getFieldOrNaN(entry, 'bay_index')), ...
    valueToString(getFieldOrNaN(entry, 'bay_start_y_in')), ...
    valueToString(getFieldOrNaN(entry, 'bay_end_y_in')), ...
    string(getFieldOrBlank(entry, 'case_name')), ...
    string(firstNonBlank(getFieldOrBlank(entry, 'panel_name'), getFieldOrBlank(entry, 'name'))), ...
    string(firstNonBlank(mode, getFieldOrBlank(entry, 'combined_mode'))), ...
    "bay_averaged_panel_buckling_recovery_with_curvature_correction", ...
    "modified_megson_panel_idealisation"];
if strlength(string(ratioName)) > 0
    rows(end, 1) = string(metric);
end
end

function rows = addPB(rows, metric, value, definition, entry, mode, recovery, ideal)
rows(end+1, :) = [string(metric), valueToString(value), string(definition), ...
    valueToString(getFieldOrNaN(entry, 'bay_index')), ...
    valueToString(getFieldOrNaN(entry, 'bay_start_y_in')), ...
    valueToString(getFieldOrNaN(entry, 'bay_end_y_in')), ...
    string(getFieldOrBlank(entry, 'case_name')), ...
    string(getFieldOrBlank(entry, 'panel_name')), ...
    string(firstNonBlank(mode, getFieldOrBlank(entry, 'combined_mode'))), ...
    string(recovery), string(ideal)];
end

function label = displayMethodName(methodName)
methodName = string(methodName);
switch methodName
    case "modified_megson_z_ratio"
        label = "modified_megson_panel_idealisation";
    case "bay_averaged_four_stress_and_shear_flow"
        label = "bay_averaged_panel_buckling_recovery";
    case "bay_averaged_four_stress_and_shear_flow_with_curved_panel_compression_correction"
        label = "bay_averaged_panel_buckling_recovery_with_curvature_correction";
    case "cylindrical_compression_addition"
        label = "cylindrical_curvature_compression_correction";
    case "combined_loads_interaction"
        label = "combined_panel_buckling_ratios";
    otherwise
        label = methodName;
end
end

function cc = aggregateCurvedCorrection(bayResults)
cc.validRange = 0; cc.highR = 0; cc.flatOnly = 0; cc.notCurved = 0;
cc.invalidGeometry = 0; cc.curvedInstances = 0; cc.applied = 0;
for i = 1:numel(bayResults)
    if ~isfield(bayResults(i).panelBuckling, 'curvedCorrection')
        continue;
    end
    c = bayResults(i).panelBuckling.curvedCorrection;
    counts = getFieldOrStruct(c, 'counts');
    cc.validRange = cc.validRange + getCount(counts, 'curvature_valid_range');
    cc.highR = cc.highR + getCount(counts, 'curvature_high_R_over_t');
    cc.flatOnly = cc.flatOnly + getCount(counts, 'flat_only_outside_validity');
    cc.notCurved = cc.notCurved + getCount(counts, 'not_curved_panel');
    cc.invalidGeometry = cc.invalidGeometry + getCount(counts, 'fallback_invalid_geometry');
    if isfield(c, 'isCurvedPanel')
        cc.curvedInstances = cc.curvedInstances + nnz(c.isCurvedPanel);
    end
    if isfield(c, 'applied')
        cc.applied = cc.applied + nnz(c.applied);
    end
end
end

function value = getFieldOrStruct(s, fieldName)
if isstruct(s) && isfield(s, fieldName) && isstruct(s.(fieldName))
    value = s.(fieldName);
else
    value = struct();
end
end

function value = getCount(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = 0;
end
end

function value = getFieldOrNaN(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = NaN;
end
end

function value = getFieldOrBlank(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = "";
end
end

function out = firstNonBlank(a, b)
if strlength(string(a)) > 0
    out = a;
else
    out = b;
end
end

function s = valueToString(value)
if isstring(value) || ischar(value)
    s = string(value);
elseif islogical(value)
    s = string(value);
elseif isnumeric(value)
    if isempty(value) || any(isnan(value(:)))
        s = "NaN";
    else
        s = string(sprintf('%.12g', value(1)));
    end
else
    s = string(value);
end
end

function s = fmt(value)
if isnumeric(value)
    s = string(sprintf('%.6g', value));
else
    s = string(value);
end
end
