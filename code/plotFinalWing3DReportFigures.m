function visualisationSummary = plotFinalWing3DReportFigures()
% plotFinalWing3DReportFigures Generate final report 3D wing-box visuals.
%
% This is a report-visualisation helper only. It rebuilds the current saved
% final candidate geometry and saved spanwise load distributions, then
% exports MATLAB-rendered PNG figures and load-gradient summary tables.

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
figuresFolder = fullfile(projectRoot, 'outputs', 'figures');
tablesFolder = fullfile(projectRoot, 'outputs', 'tables');
if ~exist(figuresFolder, 'dir'), mkdir(figuresFolder); end
if ~exist(tablesFolder, 'dir'), mkdir(tablesFolder); end

candidate = readtable(fullfile(tablesFolder, ...
    'final_candidate_summary.csv'), 'TextType', 'string');
design = finalCandidateDesign(candidate);
analysisOptions = defaultFinalAnalysisOptions();
analysisOptions.suppressPracticalGeometryWarnings = true;

designOptions = loadDesignOptions();
wing = setGeometry(design, designOptions, analysisOptions);
section = computeSectionUnitLoads(wing, analysisOptions);

makeWingBox3DStructureFigure(wing, section, candidate, figuresFolder);
[gradientSummary, loadCaseLabel] = makeWing3DLoadGradientFigures( ...
    wing, analysisOptions, figuresFolder);
writetable(gradientSummary, fullfile(tablesFolder, ...
    'final_3d_visualisation_load_gradient_summary.csv'));
makeStructuralUtilisation3DFigure(wing, tablesFolder, figuresFolder);
makeLoadCaseDistributionFigures(analysisOptions, figuresFolder, tablesFolder);
makeFailureUtilisationBars(candidate, figuresFolder);
makeWeightBreakdownPie(tablesFolder, figuresFolder);
makeClearTorsionFigure(analysisOptions, figuresFolder);
makeAdaptiveMonteCarloReportFigures(tablesFolder, figuresFolder, candidate);

visualisationSummary = struct();
visualisationSummary.status = "complete";
visualisationSummary.structureFigure = fullfile(figuresFolder, ...
    'final_wingbox_3d_structure.png');
visualisationSummary.structureTopViewFigure = fullfile(figuresFolder, ...
    'final_wingbox_3d_structure_top_view.png');
visualisationSummary.loadGradientFigure = fullfile(figuresFolder, ...
    'final_wing_3d_load_gradient.png');
visualisationSummary.bendingGradientFigure = fullfile(figuresFolder, ...
    'final_wing_3d_bending_moment_gradient.png');
visualisationSummary.shearGradientFigure = fullfile(figuresFolder, ...
    'final_wing_3d_shear_force_gradient.png');
visualisationSummary.gradientComparisonFigure = fullfile(figuresFolder, ...
    'final_wing_3d_load_gradient_comparison.png');
visualisationSummary.structuralUtilisationFigure = fullfile(figuresFolder, ...
    'final_wing_3d_structural_utilisation_by_bay.png');
visualisationSummary.verticalLoadCasesFigure = fullfile(figuresFolder, ...
    'final_load_case_vertical_distributions.png');
visualisationSummary.dragTorsionLoadCasesFigure = fullfile(figuresFolder, ...
    'final_load_case_drag_torsion_distributions.png');
visualisationSummary.rootLoadSummaryFigure = fullfile(figuresFolder, ...
    'final_load_case_root_load_summary.png');
visualisationSummary.loadGradientTable = fullfile(tablesFolder, ...
    'final_3d_visualisation_load_gradient_summary.csv');
visualisationSummary.loadCaseLabel = loadCaseLabel;
visualisationSummary.notes = ...
    "Engineering visualisation only; not a CAD or manufacturing model.";
end

function design = finalCandidateDesign(candidate)
design = baselineDesign();
design.box.frontSpar_xc = candidate.frontSpar_xc(1);
design.box.rearSpar_xc = candidate.rearSpar_xc(1);
design.stringers.section_key = candidate.stringer_section_key(1);
design.stringers.sec_no = candidate.stringer_sec_no(1);
design.sparCaps.sec_no = candidate.spar_cap_sec_no(1);
design.thickness.t3_upperMainSkin_in = candidate.main_skin_t(1);
design.thickness.t4_lowerMainSkin_in = candidate.main_skin_t(1);
design.thickness.upperNoseSkin_in = candidate.noseSkin_t(1);
design.thickness.lowerNoseSkin_in = candidate.noseSkin_t(1);
design.thickness.t5_noseSkin_in = candidate.noseSkin_t(1);
design.thickness.t1_frontSpar_in = candidate.t1_frontSpar_in(1);
design.thickness.t2_rearSpar_in = candidate.t2_rearSpar_in(1);
design.stringers.numTop = candidate.numTop(1);
design.stringers.numBottom = candidate.numBottom(1);
design.stringers.numNoseTop = candidate.numNoseTop(1);
design.stringers.numNoseBottom = candidate.numNoseBottom(1);
design.ribs.ribSpacing_in = candidate.ribSpacing_in(1);
design.columnBuckling.K_column = candidate.K_column(1);
end

function makeWingBox3DStructureFigure(wing, section, candidate, figuresFolder)
chord = wing.planform.chord_in;
span = wing.planform.semispan_in;
frontX = wing.box.frontSpar_x_in;
rearX = wing.box.rearSpar_x_in;
ribs = wing.ribs.stations_semispan_in(:)';

xSkin = linspace(frontX, rearX, 60);
xNose = linspace(0, frontX, 45);
ySpan = linspace(0, span, 80);
[XTop, YTop] = meshgrid(xSkin, ySpan);
[XNose, YNose] = meshgrid(xNose, ySpan);
ZTop = arrayfun(@(x) getNaca2412Z(x / chord, true) * chord, XTop);
ZBot = arrayfun(@(x) getNaca2412Z(x / chord, false) * chord, XTop);
ZNoseTop = arrayfun(@(x) getNaca2412Z(x / chord, true) * chord, XNose);
ZNoseBot = arrayfun(@(x) getNaca2412Z(x / chord, false) * chord, XNose);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 850]);
hold on; grid on;

hUpper = surf(XTop, YTop, ZTop, 'FaceColor', [0.45 0.70 0.92], ...
    'FaceAlpha', 0.30, 'EdgeColor', 'none', ...
    'DisplayName', 'Upper skin');
hLower = surf(XTop, YTop, ZBot, 'FaceColor', [0.38 0.62 0.82], ...
    'FaceAlpha', 0.24, 'EdgeColor', 'none', ...
    'DisplayName', 'Lower skin');
hNoseUpper = surf(XNose, YNose, ZNoseTop, 'FaceColor', [0.58 0.78 0.96], ...
    'FaceAlpha', 0.22, 'EdgeColor', 'none', ...
    'DisplayName', 'Nose skin');
surf(XNose, YNose, ZNoseBot, 'FaceColor', [0.50 0.70 0.90], ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

hFront = drawSparSurface(frontX, span, chord, [0.20 0.26 0.34], ...
    'Front spar');
hRear = drawSparSurface(rearX, span, chord, [0.12 0.17 0.24], ...
    'Rear spar');

hRib = gobjects(1);
for i = 1:numel(ribs)
    h = drawRibPlane(ribs(i), frontX, rearX, chord);
    if i == 1
        hRib = h;
        set(hRib, 'DisplayName', 'Ribs');
    else
        set(h, 'HandleVisibility', 'off');
    end
end

[hStringers, hCaps] = drawSpanwiseMembers(section, span);
drawAirfoilOutlines(chord, span);

xlabel('Chordwise x (in)');
ylabel('Span y (in)');
zlabel('Airfoil z (in)');
title(sprintf(['Final 3D Wing-Box Structural Layout: front/rear = %.3f/%.3f, ', ...
    'stringer %s, spar cap %d'], candidate.frontSpar_xc(1), ...
    candidate.rearSpar_xc(1), candidate.stringer_section_key(1), ...
    candidate.spar_cap_sec_no(1)), 'Interpreter', 'none');
legend([hNoseUpper hUpper hLower hFront hRear hRib hStringers hCaps], ...
    {'Nose skin','Upper main skin','Lower main skin','Front spar','Rear spar','Ribs', ...
    'Stringers','Spar caps'}, 'Location', 'northeastoutside');
view([38 20]);
daspect([1 2.2 0.22]);
camlight('headlight');
lighting gouraud;
axis tight;
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_wingbox_3d_structure.png'), 'Resolution', 220);

view([0 90]);
title('Final 3D Wing-Box Structural Layout - Top View');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_wingbox_3d_structure_top_view.png'), 'Resolution', 220);
close(fig);
end

function makeFailureUtilisationBars(candidate, figuresFolder)
yieldUtil = candidate.max_yield_ratio(1);
columnUtil = candidate.max_column_buckling_ratio(1);
panelUtil = candidate.max_active_combined_panel_ratio(1);
divUtil = 1 / candidate.divergence_ratio(1);
values = [yieldUtil, columnUtil, panelUtil, divUtil];
labels = categorical({'Yield','Column buckling','Panel combined','Divergence'});
labels = reordercats(labels, cellstr(labels));

fig = figure('Visible','off','Color','w','Position',[100 100 950 560]);
bar(labels, values, 0.62, 'FaceColor', [0.22 0.47 0.78]);
hold on;
yline(1.0, 'r--', 'Pass/fail limit', 'LineWidth', 1.5);
grid on;
ylabel('Utilisation ratio');
ylim([0, max(1.12, 1.08 * max(values))]);
title('Final Failure Utilisation Summary');
text(1:numel(values), values + 0.025, compose('%.3f', values), ...
    'HorizontalAlignment','center', 'FontSize', 9);
annotation('textbox', [0.16 0.80 0.74 0.08], 'String', ...
    'Divergence plotted as utilisation = 1 / divergence ratio so all checks use limit = 1.', ...
    'EdgeColor','none', 'HorizontalAlignment','center', 'FontSize', 9);
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_failure_utilisation_bars.png'), 'Resolution', 220);
close(fig);
end

function makeWeightBreakdownPie(tablesFolder, figuresFolder)
weightFile = fullfile(tablesFolder, 'final_weight_summary.csv');
opts = detectImportOptions(weightFile, 'Delimiter', ',', ...
    'VariableNamingRule','preserve');
T = readtable(weightFile, opts);
names = string(T{:, 1});
weights = T{:, 3};
wanted = ["skins_webs", "spar_caps", "stringers", "ribs"];
plotNames = ["Skins/webs", "Spar caps", "Stringers", "Ribs"];
values = zeros(numel(wanted), 1);
for i = 1:numel(wanted)
    idx = find(names == wanted(i), 1);
    if ~isempty(idx)
        values(i) = weights(idx);
    end
end

fig = figure('Visible','off','Color','w','Position',[100 100 1100 780]);
totalWeight = sum(values, 'omitnan');
if totalWeight <= 0
    warning('Weight breakdown pie skipped because component weights are unavailable.');
    close(fig);
    return;
end
percentages = 100 * values / totalWeight;
labels = plotNames(:) + " " + compose('%.1f%%', percentages(:));
p = pie(values(:), cellstr(labels(:)));
title('Final Weight Breakdown by Component');
colormap([0.42 0.66 0.88; 0.76 0.28 0.22; 0.27 0.48 0.36; 0.86 0.68 0.28]);
ax = gca;
ax.Position = [0.18 0.25 0.64 0.66];
lgd = legend(plotNames + " (" + compose('%.1f lb', values') + ")", ...
    'Orientation','horizontal', 'Box','off', 'FontSize', 9);
lgd.Position = [0.20 0.055 0.60 0.06];
ax.Position = [0.18 0.25 0.64 0.66];
for i = 2:2:numel(p)
    p(i).FontSize = 10;
    p(i).FontWeight = 'bold';
end
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_weight_breakdown_pie.png'), 'Resolution', 220);
close(fig);
end

function makeLoadCaseDistributionFigures(analysisOptions, figuresFolder, tablesFolder)
spanData = load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'spanwise_loads.mat'));
spanwise = spanData.spanwiseLoads;
y_ft = spanwise.y(:);
caseNames = string(spanwise.caseNames(:));
preferredCases = ["PHAA", "PLAA", "NLAA_VC", "NG_VC"];
caseIdx = [];
for i = 1:numel(preferredCases)
    idx = find(caseNames == preferredCases(i), 1);
    if ~isempty(idx)
        caseIdx(end+1) = idx; %#ok<AGROW>
    end
end
if isempty(caseIdx)
    caseIdx = 1:min(4, numel(caseNames));
end
loadFactor = analysisOptions.loadFactor;

makeVerticalLoadCaseFigure(spanwise, y_ft, caseNames, caseIdx, ...
    loadFactor, figuresFolder);
makeDragTorsionLoadCaseFigure(spanwise, y_ft, caseNames, caseIdx, ...
    loadFactor, figuresFolder);
makeRootLoadSummaryFigure(spanwise, caseNames, caseIdx, loadFactor, ...
    figuresFolder, tablesFolder);
end

function makeVerticalLoadCaseFigure(spanwise, y_ft, caseNames, caseIdx, ...
    loadFactor, figuresFolder)
nRows = numel(caseIdx);
fig = figure('Visible','off','Color','w','Position',[80 80 1200 260*nRows]);
tiledlayout(nRows, 3, 'TileSpacing','compact', 'Padding','compact');
for row = 1:nRows
    iCase = caseIdx(row);
    label = caseNames(iCase);
    qz = loadFactor * spanwise.distZ{iCase}.q(:);
    Vz = loadFactor * spanwise.distZ{iCase}.V(:);
    Mx = loadFactor * spanwise.distZ{iCase}.M(:);

    nexttile;
    plot(y_ft, qz, 'LineWidth', 1.35);
    grid on;
    ylabel(label);
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('q_z(y) (lb/ft)');
    xlim([min(y_ft), max(y_ft)]);

    nexttile;
    plot(y_ft, Vz, 'LineWidth', 1.35);
    grid on;
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('V_z(y) (lb)');
    xlim([min(y_ft), max(y_ft)]);

    nexttile;
    plot(y_ft, Mx, 'LineWidth', 1.35);
    grid on;
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('M_x(y) (lb-ft)');
    xlim([min(y_ft), max(y_ft)]);
end
sgtitle(sprintf('Final Applied Vertical Load Distributions, Load Factor %.1f', loadFactor));
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_load_case_vertical_distributions.png'), 'Resolution', 220);
close(fig);
end

function makeDragTorsionLoadCaseFigure(spanwise, y_ft, caseNames, caseIdx, ...
    loadFactor, figuresFolder)
nRows = numel(caseIdx);
fig = figure('Visible','off','Color','w','Position',[80 80 1450 260*nRows]);
tiledlayout(nRows, 4, 'TileSpacing','compact', 'Padding','compact');
for row = 1:nRows
    iCase = caseIdx(row);
    label = caseNames(iCase);
    qx = loadFactor * spanwise.distX{iCase}.qx(:);
    Vx = loadFactor * spanwise.distX{iCase}.Vx(:);
    Mz = loadFactor * spanwise.distX{iCase}.Mz(:);
    My = loadFactor * spanwise.distMy{iCase}.My_dist(:);

    nexttile;
    plot(y_ft, qx, 'LineWidth', 1.35);
    grid on;
    ylabel(label);
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('q_x(y) (lb/ft)');
    xlim([min(y_ft), max(y_ft)]);

    nexttile;
    plot(y_ft, Vx, 'LineWidth', 1.35);
    grid on;
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('V_x(y) (lb)');
    xlim([min(y_ft), max(y_ft)]);

    nexttile;
    plot(y_ft, Mz, 'LineWidth', 1.35);
    grid on;
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('M_z(y) (lb-ft)');
    xlim([min(y_ft), max(y_ft)]);

    nexttile;
    plot(y_ft, My, 'LineWidth', 1.35);
    grid on;
    if row == nRows, xlabel('Span station, y (ft)'); end
    title('M_y(y) (lb-ft)');
    xlim([min(y_ft), max(y_ft)]);
end
sgtitle(sprintf('Final Applied Drag and Torsion Distributions, Load Factor %.1f', loadFactor));
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_load_case_drag_torsion_distributions.png'), 'Resolution', 220);
close(fig);
end

function makeRootLoadSummaryFigure(spanwise, caseNames, caseIdx, loadFactor, ...
    figuresFolder, tablesFolder)
nCases = numel(caseIdx);
labels = caseNames(caseIdx);
rootVz = zeros(nCases, 1);
rootMx = zeros(nCases, 1);
rootVx = zeros(nCases, 1);
rootMz = zeros(nCases, 1);
rootMy = zeros(nCases, 1);
for row = 1:nCases
    iCase = caseIdx(row);
    rootVz(row) = loadFactor * spanwise.distZ{iCase}.V_root;
    rootMx(row) = loadFactor * spanwise.distZ{iCase}.M_root;
    rootVx(row) = loadFactor * spanwise.distX{iCase}.Vx_root;
    rootMz(row) = loadFactor * spanwise.distX{iCase}.Mz_root;
    rootMy(row) = loadFactor * spanwise.distMy{iCase}.My_root;
end
summary = table(labels(:), rootVz, rootMx, rootVx, rootMz, rootMy, ...
    'VariableNames', {'case_name', 'root_Vz_lb', 'root_Mx_lbft', ...
    'root_Vx_lb', 'root_Mz_lbft', 'root_My_lbft'});
writetable(summary, fullfile(tablesFolder, 'final_load_case_root_load_summary.csv'));

fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
tiledlayout(2, 1, 'TileSpacing','compact', 'Padding','compact');
nexttile;
bar(categorical(labels), [rootVz, rootVx]);
grid on;
ylabel('Root shear (lb)');
legend({'V_z', 'V_x'}, 'Location','best');
title('Final Root Shear by Critical Load Case');

nexttile;
bar(categorical(labels), [rootMx, rootMz, rootMy]);
grid on;
ylabel('Root moment (lb-ft)');
legend({'M_x', 'M_z', 'M_y'}, 'Location','best');
title('Final Root Moments by Critical Load Case');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_load_case_root_load_summary.png'), 'Resolution', 220);
close(fig);
end

function makeClearTorsionFigure(analysisOptions, figuresFolder)
spanwiseFile = load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'spanwise_loads.mat'));
spanwise = spanwiseFile.spanwiseLoads;
y_ft = spanwise.y(:);
caseNames = string({spanwise.allCases.name});
nCases = numel(caseNames);
loadFactor = analysisOptions.loadFactor;
My = zeros(numel(y_ft), nCases);
for i = 1:nCases
    My(:, i) = loadFactor * spanwise.distMy{i}.My_dist(:);
end
[~, order] = sort(max(abs(My), [], 1), 'descend');
topN = min(5, nCases);
selected = order(1:topN);
envelope = max(abs(My), [], 2);

fig = figure('Visible','off','Color','w','Position',[100 100 1050 620]);
tiledlayout(2, 1, 'TileSpacing','compact');
nexttile;
styles = {'-o','-s','-^','-d','-v'};
hold on; grid on;
for k = 1:topN
    i = selected(k);
    plot(y_ft, My(:, i), styles{k}, 'LineWidth', 1.2, ...
        'MarkerIndices', 1:80:numel(y_ft), 'DisplayName', caseNames(i));
end
xlabel('Span y (ft)');
ylabel('Torsion M_y (lb-ft)');
title('Dominant Final Applied Torsion Cases');
legend('Location','eastoutside','Interpreter','none');

nexttile;
plot(y_ft, envelope, 'k-', 'LineWidth', 1.6);
grid on;
xlabel('Span y (ft)');
ylabel('Max |M_y| envelope (lb-ft)');
title('Absolute Torsion Envelope Over All Load Cases');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_spanwise_torsion_clear.png'), 'Resolution', 220);
close(fig);
end

function makeAdaptiveMonteCarloReportFigures(tablesFolder, figuresFolder, candidate)
results = readtable(fullfile(tablesFolder, ...
    'final_adaptive_monte_carlo_results.csv'), 'TextType','string');
stageSummary = readtable(fullfile(tablesFolder, ...
    'final_adaptive_monte_carlo_stage_summary.csv'));
stages = stageSummary.stage;
bestWeight = stageSummary.best_weight_including_ribs_excluding_misc_lb;
passRate = 100 * stageSummary.pass_rate;
passCount = stageSummary.n_passing;

fig = figure('Visible','off','Color','w','Position',[100 100 1150 720]);
tiledlayout(1, 3, 'TileSpacing','compact');
nexttile; plot(stages, bestWeight, '-o', 'LineWidth', 1.6); grid on;
xlabel('Stage'); ylabel('Best passing weight (lb)');
title('Best Weight');
nexttile; plot(stages, passRate, '-s', 'LineWidth', 1.6); grid on;
xlabel('Stage'); ylabel('Pass rate (%)');
title('Pass Rate');
nexttile; bar(stages, passCount, 'FaceColor', [0.24 0.50 0.75]); grid on;
xlabel('Stage'); ylabel('Passing designs');
title('Pass Count');
sgtitle('Adaptive Monte Carlo Stage Dashboard');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_adaptive_stage_dashboard.png'), 'Resolution', 220);
close(fig);

util = max([results.max_yield_ratio, ...
    results.max_column_buckling_ratio, ...
    results.max_active_combined_panel_ratio], [], 2);
passing = results.overall_pass == 1;
fig = figure('Visible','off','Color','w','Position',[100 100 950 620]);
hold on; grid on;
hFail = scatter(results.weight_including_ribs_excluding_misc_lb(~passing), ...
    util(~passing), 10, [0.72 0.72 0.72], 'filled', ...
    'DisplayName','Failing designs');
hPass = scatter(results.weight_including_ribs_excluding_misc_lb(passing), ...
    util(passing), 16, [0.20 0.52 0.32], 'filled', ...
    'DisplayName','Passing designs');
hSelected = scatter(candidate.weight_including_ribs_excluding_misc_lb(1), ...
    max([candidate.max_yield_ratio(1), ...
    candidate.max_column_buckling_ratio(1), ...
    candidate.max_active_combined_panel_ratio(1)]), ...
    95, [0.85 0.10 0.10], 'p', 'filled', ...
    'DisplayName','Selected final design');
xLimits = xlim;
hLimit = plot(xLimits, [1.0 1.0], 'r--', 'LineWidth', 1.5, ...
    'DisplayName','Pass/fail limit');
xlim(xLimits);
xText = xLimits(1) + 0.02 * diff(xLimits);
text(xText, 1.02, 'Pass/fail limit', 'Color', [0.70 0 0], ...
    'FontSize', 9, 'FontWeight','bold', 'VerticalAlignment','bottom');
xlabel('Weight including ribs, excluding misc (lb)');
ylabel('Maximum utilisation among yield, column, and active panel checks');
title('Adaptive Monte Carlo Weight vs. Governing Utilisation');
legend([hPass hFail hSelected hLimit], ...
    {'Passing designs','Failing designs','Selected final design','Pass/fail limit'}, ...
    'Location','best');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_adaptive_weight_vs_utilisation.png'), 'Resolution', 220);
close(fig);

fig = figure('Visible','off','Color','w','Position',[100 100 950 620]);
hold on; grid on;
scatter(results.frontSpar_xc(~passing), ...
    results.weight_including_ribs_excluding_misc_lb(~passing), ...
    10, [0.70 0.70 0.70], 'filled', 'DisplayName','Failing designs');
scatter(results.frontSpar_xc(passing), ...
    results.weight_including_ribs_excluding_misc_lb(passing), ...
    16, [0.20 0.47 0.78], 'filled', 'DisplayName','Passing designs');
scatter(candidate.frontSpar_xc(1), ...
    candidate.weight_including_ribs_excluding_misc_lb(1), ...
    95, [0.85 0.10 0.10], 'p', 'filled', ...
    'DisplayName','Selected final design');
xlabel('Front spar x/c');
ylabel('Weight including ribs, excluding misc (lb)');
title('Adaptive Monte Carlo Design-Space Projection');
legend('Location','best');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_adaptive_design_space_projection.png'), 'Resolution', 220);
close(fig);
end

function h = drawSparSurface(xIn, span, chord, faceColor, displayName)
y = [0 span; 0 span];
zTop = getNaca2412Z(xIn / chord, true) * chord;
zBot = getNaca2412Z(xIn / chord, false) * chord;
x = xIn * ones(2, 2);
z = [zTop zTop; zBot zBot];
h = surf(x, y, z, 'FaceColor', faceColor, 'FaceAlpha', 0.62, ...
    'EdgeColor', faceColor * 0.55, 'LineWidth', 0.7, ...
    'DisplayName', displayName);
end

function h = drawRibPlane(yIn, frontX, rearX, chord)
x = [frontX rearX rearX frontX];
z = [getNaca2412Z(frontX / chord, true), ...
    getNaca2412Z(rearX / chord, true), ...
    getNaca2412Z(rearX / chord, false), ...
    getNaca2412Z(frontX / chord, false)] * chord;
y = yIn * ones(size(x));
h = patch(x, y, z, [0.78 0.78 0.78], 'FaceAlpha', 0.20, ...
    'EdgeColor', [0.55 0.55 0.55], 'LineWidth', 0.7);
end

function [hStringers, hCaps] = drawSpanwiseMembers(section, span)
hStringers = gobjects(1);
hCaps = gobjects(1);

stringerGroups = ["top", "bottom", "noseTop", "noseBottom"];
firstStringer = true;
for groupName = stringerGroups
    if isfield(section.stringers, groupName)
        booms = section.stringers.(groupName);
        for i = 1:numel(booms)
            h = plot3([booms(i).x_in booms(i).x_in], [0 span], ...
                [booms(i).z_in booms(i).z_in], '-', ...
                'Color', [0.02 0.08 0.16], 'LineWidth', 1.3);
            if firstStringer
                hStringers = h;
                set(hStringers, 'DisplayName', 'Stringers');
                firstStringer = false;
            else
                set(h, 'HandleVisibility', 'off');
            end
        end
    end
end

firstCap = true;
for i = 1:numel(section.sparCaps)
    h = plot3([section.sparCaps(i).x_in section.sparCaps(i).x_in], ...
        [0 span], [section.sparCaps(i).z_in section.sparCaps(i).z_in], ...
        '-', 'Color', [0.75 0.10 0.10], 'LineWidth', 2.2);
    if firstCap
        hCaps = h;
        set(hCaps, 'DisplayName', 'Spar caps');
        firstCap = false;
    else
        set(h, 'HandleVisibility', 'off');
    end
end
end

function drawAirfoilOutlines(chord, span)
xc = linspace(0, 1, 150);
x = xc * chord;
zu = arrayfun(@(v) getNaca2412Z(v, true), xc) * chord;
zl = arrayfun(@(v) getNaca2412Z(v, false), xc) * chord;
for y = [0 span]
    plot3(x, y * ones(size(x)), zu, 'k-', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    plot3(x, y * ones(size(x)), zl, 'k-', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end
end

function [gradientSummary, loadCaseLabel] = makeWing3DLoadGradientFigures( ...
    wing, analysisOptions, figuresFolder)
spanwiseFile = load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'spanwise_loads.mat'));
spanwise = spanwiseFile.spanwiseLoads;
loadFactor = analysisOptions.loadFactor;

y_ft = spanwise.y(:);
y_in = y_ft * 12;
caseNames = string({spanwise.allCases.name})';
nCases = numel(caseNames);
Mx = zeros(numel(y_ft), nCases);
Mz = zeros(numel(y_ft), nCases);
Vx = zeros(numel(y_ft), nCases);
Vz = zeros(numel(y_ft), nCases);
My = zeros(numel(y_ft), nCases);
for i = 1:nCases
    Mx(:, i) = loadFactor * spanwise.distZ{i}.M(:);
    Mz(:, i) = loadFactor * spanwise.distX{i}.Mz(:);
    Vx(:, i) = loadFactor * spanwise.distX{i}.Vx(:);
    Vz(:, i) = loadFactor * spanwise.distZ{i}.V(:);
    My(:, i) = loadFactor * spanwise.distMy{i}.My_dist(:);
end

bendingEnvelope = max(sqrt(Mx.^2 + Mz.^2), [], 2);
shearEnvelope = max(sqrt(Vx.^2 + Vz.^2), [], 2);
[~, bendCaseIdx] = max(sqrt(Mx.^2 + Mz.^2), [], 2);
[~, shearCaseIdx] = max(sqrt(Vx.^2 + Vz.^2), [], 2);
torsionEnvelope = max(abs(My), [], 2);

gradientSummary = table(y_in, y_ft, bendingEnvelope, ...
    caseNames(bendCaseIdx), shearEnvelope, caseNames(shearCaseIdx), ...
    torsionEnvelope, repmat(loadFactor, numel(y_in), 1), ...
    repmat("envelope_over_all_load_cases", numel(y_in), 1), ...
    'VariableNames', {'span_station_in', 'span_station_ft', ...
    'bending_moment_envelope_lbft', 'bending_controlling_case', ...
    'shear_force_envelope_lb', 'shear_controlling_case', ...
    'torsion_envelope_lbft', 'loadFactor', 'value_type'});
loadCaseLabel = "maximum envelope over saved load cases";

makeLoadSurfaceFigure(wing, y_in, bendingEnvelope, ...
    'Final Wing 3D Load Gradient: Bending Moment Envelope', ...
    'Bending moment envelope (lb-ft)', ...
    fullfile(figuresFolder, 'final_wing_3d_load_gradient.png'));
makeLoadSurfaceFigure(wing, y_in, bendingEnvelope, ...
    'Final Wing 3D Bending Moment Gradient', ...
    'Bending moment envelope (lb-ft)', ...
    fullfile(figuresFolder, 'final_wing_3d_bending_moment_gradient.png'));
makeLoadSurfaceFigure(wing, y_in, shearEnvelope, ...
    'Final Wing 3D Shear Force Gradient', ...
    'Shear force envelope (lb)', ...
    fullfile(figuresFolder, 'final_wing_3d_shear_force_gradient.png'));
makeLoadGradientComparisonFigure(wing, y_in, bendingEnvelope, ...
    shearEnvelope, figuresFolder);
end

function makeLoadSurfaceFigure(wing, y_in, scalar, titleText, cbarLabel, fileName)
chord = wing.planform.chord_in;
frontX = wing.box.frontSpar_x_in;
rearX = wing.box.rearSpar_x_in;
x = linspace(frontX, rearX, 50);
[X, Y] = meshgrid(x, y_in);
ZTop = arrayfun(@(xVal) getNaca2412Z(xVal / chord, true) * chord, X);
C = repmat(scalar(:), 1, numel(x));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 760]);
surf(X, Y, ZTop, C, 'EdgeColor', 'none', 'FaceAlpha', 0.92);
hold on;
plot3([frontX frontX], [0 max(y_in)], ...
    [getNaca2412Z(frontX/chord, true) getNaca2412Z(frontX/chord, true)]*chord, ...
    'k-', 'LineWidth', 2);
plot3([rearX rearX], [0 max(y_in)], ...
    [getNaca2412Z(rearX/chord, true) getNaca2412Z(rearX/chord, true)]*chord, ...
    'k-', 'LineWidth', 2);
grid on;
colormap(turbo);
cb = colorbar;
cb.Label.String = cbarLabel;
title(titleText);
xlabel('Chordwise x (in)');
ylabel('Span y (in)');
zlabel('Upper surface z (in)');
view([38 24]);
daspect([1 2.0 0.20]);
camlight('headlight');
lighting gouraud;
axis tight;
exportgraphics(fig, fileName, 'Resolution', 220);
close(fig);
end

function makeLoadGradientComparisonFigure(wing, y_in, bendingEnvelope, ...
    shearEnvelope, figuresFolder)
chord = wing.planform.chord_in;
frontX = wing.box.frontSpar_x_in;
rearX = wing.box.rearSpar_x_in;
x = linspace(frontX, rearX, 50);
[X, Y] = meshgrid(x, y_in);
ZTop = arrayfun(@(xVal) getNaca2412Z(xVal / chord, true) * chord, X);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
surf(X, Y, ZTop, repmat(bendingEnvelope(:), 1, numel(x)), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.94);
grid on; colormap(gca, turbo);
cb = colorbar; cb.Label.String = 'Bending moment envelope (lb-ft)';
title('Bending Moment Envelope');
xlabel('Chordwise x (in)'); ylabel('Span y (in)'); zlabel('z (in)');
view([38 24]); daspect([1 2.0 0.20]); camlight('headlight'); lighting gouraud; axis tight;

nexttile;
surf(X, Y, ZTop, repmat(shearEnvelope(:), 1, numel(x)), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.94);
grid on; colormap(gca, parula);
cb = colorbar; cb.Label.String = 'Shear force envelope (lb)';
title('Shear Force Envelope');
xlabel('Chordwise x (in)'); ylabel('Span y (in)'); zlabel('z (in)');
view([38 24]); daspect([1 2.0 0.20]); camlight('headlight'); lighting gouraud; axis tight;
sgtitle('Final Applied Load Gradient Comparison, Load Factor = 1.5');
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_wing_3d_load_gradient_comparison.png'), 'Resolution', 220);
close(fig);
end

function makeStructuralUtilisation3DFigure(wing, tablesFolder, figuresFolder)
bayTable = readtable(fullfile(tablesFolder, ...
    'final_spanwise_panel_buckling_by_bay.csv'));
if ~all(ismember({'y_start_in','y_end_in','max_active_combined'}, ...
        bayTable.Properties.VariableNames))
    return;
end

chord = wing.planform.chord_in;
frontX = wing.box.frontSpar_x_in;
rearX = wing.box.rearSpar_x_in;
x = linspace(frontX, rearX, 44);
y = unique([bayTable.y_start_in; bayTable.y_end_in])';
if numel(y) < 2
    return;
end
[X, Y] = meshgrid(x, y);
ZTop = arrayfun(@(xVal) getNaca2412Z(xVal / chord, true) * chord, X);
C = zeros(size(Y));
for i = 1:height(bayTable)
    idx = Y >= bayTable.y_start_in(i) & Y <= bayTable.y_end_in(i);
    C(idx) = bayTable.max_active_combined(i);
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1180 720]);
surf(X, Y, ZTop, C, 'EdgeColor', [0.55 0.55 0.55], 'FaceAlpha', 0.95);
hold on;
for i = 1:height(bayTable)
    y0 = bayTable.y_start_in(i);
    plot3([frontX rearX], [y0 y0], ...
        [getNaca2412Z(frontX/chord,true)*chord, getNaca2412Z(rearX/chord,true)*chord], ...
        'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
grid on; colormap(turbo);
cb = colorbar; cb.Label.String = 'Active combined panel utilisation';
caxis([0, max(1, max(bayTable.max_active_combined))]);
title('Final Structural Panel Utilisation by Bay');
xlabel('Chordwise x (in)'); ylabel('Span y (in)'); zlabel('Upper surface z (in)');
view([38 24]); daspect([1 2.0 0.20]); camlight('headlight'); lighting gouraud; axis tight;
exportgraphics(fig, fullfile(figuresFolder, ...
    'final_wing_3d_structural_utilisation_by_bay.png'), 'Resolution', 220);
close(fig);
end

function zc = getNaca2412Z(xc, isUpper)
xc = min(max(xc, 0), 1);
m = 0.02;
p = 0.40;
t = 0.12;
yt = 5 * t * (0.2969 * sqrt(xc) - 0.1260 * xc - ...
    0.3516 * xc.^2 + 0.2843 * xc.^3 - 0.1015 * xc.^4);
if xc < p
    yc = m / p^2 * (2 * p * xc - xc.^2);
    dyc = 2 * m / p^2 * (p - xc);
else
    yc = m / (1 - p)^2 * ((1 - 2*p) + 2*p*xc - xc.^2);
    dyc = 2 * m / (1 - p)^2 * (p - xc);
end
theta = atan(dyc);
if isUpper
    zc = yc + yt * cos(theta);
else
    zc = yc - yt * cos(theta);
end
end
