function exportAdaptiveMonteCarloResults(adaptiveResults, outputRoot)
% exportAdaptiveMonteCarloResults Write staged adaptive Monte Carlo outputs.

if nargin < 1 || isempty(adaptiveResults) || ~isstruct(adaptiveResults)
    error('exportAdaptiveMonteCarloResults:MissingResults', ...
        'adaptiveResults struct is required.');
end

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
if nargin < 2 || isempty(outputRoot)
    tablesFolder = fullfile(projectRoot, 'outputs', 'tables');
    figuresFolder = fullfile(projectRoot, 'outputs', 'figures');
else
    tablesFolder = char(outputRoot);
    figuresFolder = char(outputRoot);
end
if ~isfolder(tablesFolder), mkdir(tablesFolder); end
if ~isfolder(figuresFolder), mkdir(figuresFolder); end

resultsFile = fullfile(tablesFolder, 'final_adaptive_monte_carlo_results.csv');
stageSummaryFile = fullfile(tablesFolder, 'final_adaptive_monte_carlo_stage_summary.csv');
bestByStageFile = fullfile(tablesFolder, 'final_adaptive_monte_carlo_best_by_stage.csv');
rangeHistoryFile = fullfile(tablesFolder, 'final_adaptive_range_history.csv');

writetable(adaptiveResults.allDesigns, resultsFile);
writetable(adaptiveResults.stageSummary, stageSummaryFile);
writetable(adaptiveResults.bestByStage, bestByStageFile);
writetable(adaptiveResults.rangeHistory, rangeHistoryFile);

makeConvergenceFigure(adaptiveResults.stageSummary, figuresFolder);
makePassRateFigure(adaptiveResults.stageSummary, figuresFolder);
makePassCountFigure(adaptiveResults.stageSummary, figuresFolder);

fprintf('Adaptive Monte Carlo outputs exported:\n');
fprintf('- %s\n', resultsFile);
fprintf('- %s\n', stageSummaryFile);
fprintf('- %s\n', bestByStageFile);
fprintf('- %s\n', rangeHistoryFile);
fprintf('- %s\n', fullfile(figuresFolder, 'final_adaptive_convergence.png'));
fprintf('- %s\n', fullfile(figuresFolder, 'final_adaptive_pass_rate.png'));
fprintf('- %s\n', fullfile(figuresFolder, 'final_adaptive_pass_count_by_stage.png'));
end

function makeConvergenceFigure(stageSummary, figuresFolder)
fig = figure('Visible','off','Color','w');
plot(stageSummary.stage, stageSummary.best_weight_including_ribs_excluding_misc_lb, ...
    'o-', 'LineWidth', 1.5);
grid on;
xlabel('Adaptive Monte Carlo stage');
ylabel('Best passing weight including ribs/excluding misc (lb)');
title('Adaptive Monte Carlo Convergence');
exportgraphics(fig, fullfile(figuresFolder, 'final_adaptive_convergence.png'), ...
    'Resolution', 200);
close(fig);
end

function makePassRateFigure(stageSummary, figuresFolder)
fig = figure('Visible','off','Color','w');
plot(stageSummary.stage, 100 * stageSummary.pass_rate, 's-', 'LineWidth', 1.5);
grid on;
xlabel('Adaptive Monte Carlo stage');
ylabel('Pass rate (%)');
title('Adaptive Monte Carlo Pass Rate');
exportgraphics(fig, fullfile(figuresFolder, 'final_adaptive_pass_rate.png'), ...
    'Resolution', 200);
close(fig);
end

function makePassCountFigure(stageSummary, figuresFolder)
fig = figure('Visible','off','Color','w');
bar(stageSummary.stage, stageSummary.n_passing);
grid on;
xlabel('Adaptive Monte Carlo stage');
ylabel('Passing wings');
title('Adaptive Monte Carlo Passing Designs by Stage');
exportgraphics(fig, fullfile(figuresFolder, 'final_adaptive_pass_count_by_stage.png'), ...
    'Resolution', 200);
close(fig);
end
