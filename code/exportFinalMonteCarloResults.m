function exportFinalMonteCarloResults(results)
% exportFinalMonteCarloResults Write final Monte Carlo CSV artifacts.

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
tablesFolder = fullfile(projectRoot, 'outputs', 'tables');
if ~isfolder(tablesFolder)
    mkdir(tablesFolder);
end

resultsFile = fullfile(tablesFolder, 'final_monte_carlo_results.csv');
summaryFile = fullfile(tablesFolder, 'final_monte_carlo_summary.csv');
bestFile = fullfile(tablesFolder, 'final_monte_carlo_best_design.csv');

writetable(results.allDesigns, resultsFile);
writetable(results.summary, summaryFile);
writetable(results.bestDesign, bestFile);

fprintf('Final Monte Carlo outputs exported:\n');
fprintf('- %s\n', resultsFile);
fprintf('- %s\n', summaryFile);
fprintf('- %s\n', bestFile);
end
