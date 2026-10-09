function updateInfo = updateFinalCandidateArtifactsIfImproved(adaptiveResults, analysisOptions, updateOnlyIfImproved)
% updateFinalCandidateArtifactsIfImproved Protect the stored final candidate.
%
% This function only controls report-artifact overwrite behavior after an
% adaptive Monte Carlo run. It does not alter optimisation, analysis, or
% failure-screening logic.

if nargin < 3
    updateOnlyIfImproved = false;
end

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptFolder, '..');
tablesFolder = fullfile(projectRoot, 'outputs', 'tables');
archiveFolder = fullfile(projectRoot, 'outputs', 'archive', 'final_candidates');
if ~isfolder(tablesFolder), mkdir(tablesFolder); end
if ~isfolder(archiveFolder), mkdir(archiveFolder); end

runId = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
candidateFile = fullfile(tablesFolder, 'final_candidate_summary.csv');
historyFile = fullfile(tablesFolder, 'final_adaptive_run_history.csv');
bestOverallFile = fullfile(tablesFolder, 'final_best_overall_candidate.csv');

[previousCandidate, previousStatus] = readCandidateSummary(candidateFile);
previousWeight = tableWeight(previousCandidate);
hasPrevious = isfinite(previousWeight);

hasNewCandidate = isfield(adaptiveResults, 'bestDesign') && ...
    istable(adaptiveResults.bestDesign) && height(adaptiveResults.bestDesign) > 0 && ...
    isfield(adaptiveResults, 'bestDesignStruct') && ~isempty(adaptiveResults.bestDesignStruct) && ...
    isfield(adaptiveResults, 'bestSpanwiseResult') && ~isempty(adaptiveResults.bestSpanwiseResult);

newWeight = NaN;
newPass = false;
if hasNewCandidate
    newWeight = adaptiveResults.bestDesign.weight_including_ribs_excluding_misc_lb(1);
    if ismember("overall_pass", string(adaptiveResults.bestDesign.Properties.VariableNames))
        newPass = logical(adaptiveResults.bestDesign.overall_pass(1));
    elseif isfield(adaptiveResults.bestSpanwiseResult, 'summary') && ...
            isfield(adaptiveResults.bestSpanwiseResult.summary, 'overall_pass')
        newPass = logical(adaptiveResults.bestSpanwiseResult.summary.overall_pass);
    end
end

if isfield(adaptiveResults, 'stage4BestDesign') && istable(adaptiveResults.stage4BestDesign) && ...
        height(adaptiveResults.stage4BestDesign) > 0
    writetable(adaptiveResults.stage4BestDesign, ...
        fullfile(archiveFolder, runId + "_stage4_best_candidate.csv"));
end
if hasNewCandidate
    writetable(adaptiveResults.bestDesign, ...
        fullfile(archiveFolder, runId + "_global_best_candidate_used_for_update.csv"));
end
if isfield(adaptiveResults, 'stageSummary') && istable(adaptiveResults.stageSummary)
    writetable(adaptiveResults.stageSummary, ...
        fullfile(archiveFolder, runId + "_adaptive_stage_summary.csv"));
end

shouldUpdate = false;
reason = "";
if ~hasNewCandidate || ~newPass
    reason = "No passing candidate from current run.";
elseif ~updateOnlyIfImproved
    shouldUpdate = true;
    reason = "UPDATE_FINAL_ONLY_IF_IMPROVED is false; active final candidate overwritten by global best across all stages.";
elseif ~hasPrevious
    shouldUpdate = true;
    reason = "No readable previous active final candidate; global best across all stages becomes active.";
elseif newWeight < previousWeight
    shouldUpdate = true;
    reason = "Global best across all stages is lighter than previous active final candidate.";
else
    reason = "New run global best is heavier than the stored final candidate; active final candidate not overwritten.";
end

if shouldUpdate
    fprintf('\nRegenerating active final candidate artifacts from global best across all stages...\n');
    exportFinalCandidateArtifacts(adaptiveResults.bestDesignStruct, ...
        adaptiveResults.bestSpanwiseResult, analysisOptions, adaptiveResults);
    exportAdaptiveMonteCarloResults(adaptiveResults);
    adaptiveOutputsFolder = tablesFolder;
    [activeCandidateAfter, afterStatus] = readCandidateSummary(candidateFile);
else
    fprintf('\n%s\n', reason);
    diagnosticsFolder = fullfile(projectRoot, 'outputs', 'diagnostics', ...
        'archive', runId);
    exportAdaptiveMonteCarloResults(adaptiveResults, diagnosticsFolder);
    adaptiveOutputsFolder = diagnosticsFolder;
    fprintf('Current run adaptive outputs saved for diagnostics only: %s\n', ...
        diagnosticsFolder);
    activeCandidateAfter = previousCandidate;
    afterStatus = previousStatus;
end

afterWeight = tableWeight(activeCandidateAfter);
writeBestOverallCandidate(bestOverallFile, activeCandidateAfter, activeCandidateAfter);
appendRunHistory(historyFile, adaptiveResults, runId, updateOnlyIfImproved, ...
    previousWeight, afterWeight, shouldUpdate, reason, previousStatus, afterStatus, ...
    newWeight);

fprintf('Final candidate update policy: update only if improved = %s\n', ...
    string(logical(updateOnlyIfImproved)));
fprintf('Previous active final candidate weight: %s lb\n', formatWeight(previousWeight));
fprintf('Current run candidate used for update: %s\n', candidateUsedForUpdate(adaptiveResults));
fprintf('Current run candidate-used weight: %s lb\n', formatWeight(newWeight));
fprintf('Active final candidate weight after run: %s lb\n', formatWeight(afterWeight));
fprintf('Active final candidate updated: %s\n', string(shouldUpdate));
fprintf('Run history updated: %s\n', historyFile);
    fprintf('Best overall candidate table updated: %s\n', bestOverallFile);

updateInfo = struct();
updateInfo.run_id = runId;
updateInfo.updateOnlyIfImproved = logical(updateOnlyIfImproved);
updateInfo.previousWeight = previousWeight;
updateInfo.newWeight = newWeight;
updateInfo.afterWeight = afterWeight;
updateInfo.updated = shouldUpdate;
updateInfo.reason = reason;
updateInfo.previousStatus = previousStatus;
updateInfo.afterStatus = afterStatus;
updateInfo.adaptiveOutputsFolder = string(adaptiveOutputsFolder);
end

function [candidate, status] = readCandidateSummary(candidateFile)
candidate = table();
status = "missing";
if ~isfile(candidateFile)
    return;
end
try
    candidate = readtable(candidateFile, 'TextType', 'string');
    if height(candidate) > 0
        candidate = candidate(1, :);
        status = "readable";
    else
        status = "empty";
    end
catch ME
    candidate = table();
    status = "unreadable: " + string(ME.identifier);
end
end

function weight = tableWeight(T)
weight = NaN;
if istable(T) && height(T) > 0 && ...
        ismember("weight_including_ribs_excluding_misc_lb", string(T.Properties.VariableNames))
    value = T.weight_including_ribs_excluding_misc_lb(1);
    if isnumeric(value)
        weight = double(value);
    else
        weight = str2double(string(value));
    end
end
end

function writeBestOverallCandidate(bestOverallFile, activeCandidateAfter, fallbackCandidate)
if istable(activeCandidateAfter) && height(activeCandidateAfter) > 0
    writetable(activeCandidateAfter, bestOverallFile);
elseif istable(fallbackCandidate) && height(fallbackCandidate) > 0
    writetable(fallbackCandidate, bestOverallFile);
else
    T = table("no_active_candidate", "No readable active final candidate is available.", ...
        'VariableNames', {'status', 'notes'});
    writetable(T, bestOverallFile);
end
end

function appendRunHistory(historyFile, adaptiveResults, runId, updateOnlyIfImproved, ...
    previousWeight, afterWeight, wasUpdated, reason, previousStatus, afterStatus, ...
    stage4BestWeight)

nPerStageText = "";
seedMode = "";
if isfield(adaptiveResults, 'optimisationOptions')
    opt = adaptiveResults.optimisationOptions;
    if isfield(opt, 'adaptive') && isfield(opt.adaptive, 'nPerStage')
        nPerStageText = string(mat2str(opt.adaptive.nPerStage));
    end
    if isfield(opt, 'seedMode')
        seedMode = string(opt.seedMode);
    end
end

storedSeed = NaN;
globalBestWeight = NaN;
if isfield(adaptiveResults, 'stageSummary') && istable(adaptiveResults.stageSummary) && ...
        height(adaptiveResults.stageSummary) > 0
    storedSeed = adaptiveResults.stageSummary.rngSeed(end);
end
if isfield(adaptiveResults, 'bestOverallDesign') && istable(adaptiveResults.bestOverallDesign) && ...
        height(adaptiveResults.bestOverallDesign) > 0
    globalBestWeight = adaptiveResults.bestOverallDesign.weight_including_ribs_excluding_misc_lb(1);
end
if isfield(adaptiveResults, 'stage4BestDesign') && istable(adaptiveResults.stage4BestDesign) && ...
        height(adaptiveResults.stage4BestDesign) > 0
    stage4BestWeight = adaptiveResults.stage4BestDesign.weight_including_ribs_excluding_misc_lb(1);
end
candidateUsed = candidateUsedForUpdate(adaptiveResults);

newRow = table( ...
    string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), runId, ...
    nPerStageText, seedMode, storedSeed, stage4BestWeight, globalBestWeight, ...
    candidateUsed, previousWeight, afterWeight, logical(wasUpdated), logical(updateOnlyIfImproved), ...
    string(previousStatus), string(afterStatus), string(reason), ...
    'VariableNames', {'timestamp','run_id','N_PER_STAGE','seed_mode', ...
    'stored_rng_seed','stage4_best_weight_including_ribs_excluding_misc_lb', ...
    'global_best_across_stages_weight_including_ribs_excluding_misc_lb', ...
    'candidate_used_for_update', ...
    'active_final_candidate_weight_before_run', ...
    'active_final_candidate_weight_after_run', 'active_candidate_updated', ...
    'update_final_only_if_improved', 'previous_candidate_status', ...
    'after_candidate_status', 'notes'});

if isfile(historyFile)
    try
        existing = readtable(historyFile, 'TextType', 'string');
        if ~ismember("candidate_used_for_update", string(existing.Properties.VariableNames))
            existing.candidate_used_for_update = repmat("stage4_best_reference", height(existing), 1);
        end
        if all(ismember(string(newRow.Properties.VariableNames), string(existing.Properties.VariableNames)))
            existing = existing(:, newRow.Properties.VariableNames);
            T = [existing; newRow];
        else
            T = newRow;
        end
    catch
        T = newRow;
    end
else
    T = newRow;
end
writetable(T, historyFile);
end

function text = formatWeight(weight)
if isfinite(weight)
    text = string(sprintf('%.4f', weight));
else
    text = "NaN";
end
end

function value = candidateUsedForUpdate(adaptiveResults)
if isfield(adaptiveResults, 'candidateUsedForUpdate')
    value = string(adaptiveResults.candidateUsedForUpdate);
else
    value = "global_best_across_all_stages";
end
end
