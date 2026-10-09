function [nextRanges, historyRows] = refineFinalMonteCarloRanges(currentRanges, eliteDesigns, initialRanges, refineOptions, stageNumber)
% refineFinalMonteCarloRanges Narrow Monte Carlo ranges from elite designs.
%
% The refinement operates on the current candidate lists by index and never
% creates values outside the initial Stage 1 candidate lists.

if nargin < 5
    stageNumber = NaN;
end
validateInputs(currentRanges, eliteDesigns, initialRanges, refineOptions);

nextRanges = currentRanges;
historyRows = repmat(makeHistoryRow(), 0, 1);

numericVars = ["frontSpar_xc", "rearSpar_xc", "main_skin_t", ...
    "noseSkin_t", "t1_frontSpar_in", "t2_rearSpar_in", "numTop", ...
    "numBottom", "numNoseTop", "numNoseBottom", "ribSpacing_in"];

for iVar = 1:numel(numericVars)
    varName = numericVars(iVar);
    [newValues, row] = refineIndexedRange( ...
        currentRanges.(varName), eliteDesigns.(varName), ...
        initialRanges.(varName), refineOptions.rangePaddingSteps, ...
        stageNumber, varName, 2);
    nextRanges.(varName) = newValues;
    historyRows(end+1, 1) = row; %#ok<AGROW>
end

[newValues, row] = refineIndexedRange( ...
    currentRanges.stringer_section_key, eliteDesigns.stringer_section_key, ...
    initialRanges.stringer_section_key, refineOptions.rangePaddingSteps, ...
    stageNumber, "stringer_section_key", refineOptions.minSectionCandidates);
nextRanges.stringer_section_key = string(newValues(:));
historyRows(end+1, 1) = row;

[newValues, row] = refineIndexedRange( ...
    currentRanges.spar_cap_sec_no, eliteDesigns.spar_cap_sec_no, ...
    initialRanges.spar_cap_sec_no, refineOptions.rangePaddingSteps, ...
    stageNumber, "spar_cap_sec_no", refineOptions.minSectionCandidates);
nextRanges.spar_cap_sec_no = newValues;
historyRows(end+1, 1) = row;

nextRanges.minSparSeparation_xc = currentRanges.minSparSeparation_xc;
nextRanges.K_column = currentRanges.K_column;
end

function validateInputs(currentRanges, eliteDesigns, initialRanges, refineOptions)
if ~isstruct(currentRanges) || ~isstruct(initialRanges)
    error('refineFinalMonteCarloRanges:InvalidRanges', ...
        'currentRanges and initialRanges must be structs.');
end
if ~istable(eliteDesigns) || height(eliteDesigns) < 1
    error('refineFinalMonteCarloRanges:NoEliteDesigns', ...
        'eliteDesigns must be a nonempty table.');
end
required = ["rangePaddingSteps", "minSectionCandidates"];
missing = setdiff(required, string(fieldnames(refineOptions)));
if ~isempty(missing)
    error('refineFinalMonteCarloRanges:MissingOptions', ...
        'refineOptions is missing fields: %s', strjoin(missing, ', '));
end
end

function [newValues, row] = refineIndexedRange(parentValues, eliteValues, initialValues, paddingSteps, stageNumber, varName, minCount)
parentValues = normalizeVector(parentValues);
initialValues = normalizeVector(initialValues);
eliteValues = normalizeVector(eliteValues);

if isempty(parentValues)
    error('refineFinalMonteCarloRanges:EmptyParentRange', ...
        'Parent range for %s is empty.', varName);
end

eliteMask = false(size(parentValues));
for iElite = 1:numel(eliteValues)
    matches = valuesEqual(parentValues, eliteValues(iElite));
    eliteMask = eliteMask | matches;
end
eliteIdx = find(eliteMask);
if isempty(eliteIdx)
    error('refineFinalMonteCarloRanges:EliteOutsideParentRange', ...
        'No elite values for %s were found in the current parent range.', varName);
end

iMin = max(1, min(eliteIdx) - paddingSteps);
iMax = min(numel(parentValues), max(eliteIdx) + paddingSteps);
[iMin, iMax] = enforceMinimumCount(iMin, iMax, numel(parentValues), minCount);
candidateValues = parentValues(iMin:iMax);

initialMask = false(size(candidateValues));
for iValue = 1:numel(candidateValues)
    initialMask(iValue) = any(valuesEqual(initialValues, candidateValues(iValue)));
end
newValues = candidateValues(initialMask);
if isempty(newValues)
    error('refineFinalMonteCarloRanges:NoValuesAfterInitialClip', ...
        'Refined range for %s has no values after clipping to initial range.', varName);
end

row = makeHistoryRow();
row.stage = stageNumber;
row.variable = string(varName);
row.range_before = rangeToString(parentValues);
row.elite_min = valueToString(parentValues(min(eliteIdx)));
row.elite_max = valueToString(parentValues(max(eliteIdx)));
row.range_after = rangeToString(newValues);
row.n_parent_values = numel(parentValues);
row.n_after_values = numel(newValues);
end

function [iMin, iMax] = enforceMinimumCount(iMin, iMax, nParent, minCount)
if minCount <= 1 || nParent <= 1
    return;
end
targetCount = min(minCount, nParent);
while (iMax - iMin + 1) < targetCount
    canGrowLeft = iMin > 1;
    canGrowRight = iMax < nParent;
    if canGrowLeft
        iMin = iMin - 1;
    end
    if (iMax - iMin + 1) >= targetCount
        break;
    end
    if canGrowRight
        iMax = iMax + 1;
    end
    if ~canGrowLeft && ~canGrowRight
        break;
    end
end
end

function row = makeHistoryRow()
row = struct();
row.stage = NaN;
row.variable = "";
row.range_before = "";
row.elite_min = "";
row.elite_max = "";
row.range_after = "";
row.n_parent_values = NaN;
row.n_after_values = NaN;
end

function values = normalizeVector(values)
if isstring(values) || ischar(values) || iscellstr(values)
    values = string(values);
    values = values(:);
else
    values = values(:)';
end
end

function tf = valuesEqual(values, value)
if isstring(values) || ischar(values) || isstring(value) || ischar(value)
    tf = string(values) == string(value);
else
    tf = values == value;
end
end

function text = rangeToString(values)
values = normalizeVector(values);
if isstring(values)
    text = strjoin(values(:)', " ");
else
    text = string(mat2str(values));
end
end

function text = valueToString(value)
if isstring(value) || ischar(value)
    text = string(value);
else
    text = string(num2str(value, '%.15g'));
end
end
