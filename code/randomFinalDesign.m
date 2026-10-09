function design = randomFinalDesign(designOptions, optimisationOptions)
% randomFinalDesign Generate one final-design random candidate design.

if nargin < 1 || isempty(designOptions) || ~isstruct(designOptions)
    error('randomFinalDesign:MissingDesignOptions', ...
        'designOptions struct from loadDesignOptions is required.');
end
if nargin < 2 || isempty(optimisationOptions)
    optimisationOptions = defaultFinalOptimisationOptions(designOptions);
end
validateOptions(optimisationOptions);

ranges = optimisationOptions.ranges;
design = baselineDesign();

[frontSpar_xc, rearSpar_xc] = randomValidSparPair(ranges);
design.box.frontSpar_xc = frontSpar_xc;
design.box.rearSpar_xc = rearSpar_xc;

design.stringers.section_key = randomChoice(ranges.stringer_section_key);
design.stringers.sec_no = stringerSecNoFromKey(designOptions, design.stringers.section_key);
design.sparCaps.sec_no = randomChoice(ranges.spar_cap_sec_no);

mainSkin_t = randomChoice(ranges.main_skin_t);
noseSkin_t = randomChoice(ranges.noseSkin_t);
design.thickness.t3_upperMainSkin_in = mainSkin_t;
design.thickness.t4_lowerMainSkin_in = mainSkin_t;
design.thickness.upperNoseSkin_in = noseSkin_t;
design.thickness.lowerNoseSkin_in = noseSkin_t;
design.thickness.t5_noseSkin_in = noseSkin_t;
design.thickness.t1_frontSpar_in = randomChoice(ranges.t1_frontSpar_in);
design.thickness.t2_rearSpar_in = randomChoice(ranges.t2_rearSpar_in);

design.stringers.numTop = randomChoice(ranges.numTop);
design.stringers.numBottom = randomChoice(ranges.numBottom);
design.stringers.numNoseTop = randomChoice(ranges.numNoseTop);
design.stringers.numNoseBottom = randomChoice(ranges.numNoseBottom);

design.ribs.ribSpacing_in = randomChoice(ranges.ribSpacing_in);
design.ribs.weightModel = "geometric_airfoil_area";
design.columnBuckling.K_column = ranges.K_column;

validateGeneratedDesign(design, designOptions, ranges);

design.generation.method = "randomFinalDesign";
design.generation.objective = optimisationOptions.objective;
design.generation.ranges = ranges;
design.generation.notes = [
    "Generated for final-design Monte Carlo using full spanwise final analysis."
    "Rib weight uses geometric_airfoil_area from baselineDesign."
    "K_column fixed at 0.7."
    ];
end

function validateOptions(optimisationOptions)
requiredFields = ["frontSpar_xc", "rearSpar_xc", "minSparSeparation_xc", ...
    "stringer_section_key", "spar_cap_sec_no", "main_skin_t", "noseSkin_t", ...
    "t1_frontSpar_in", "t2_rearSpar_in", "numTop", "numBottom", ...
    "numNoseTop", "numNoseBottom", "ribSpacing_in", "K_column"];
missing = setdiff(requiredFields, string(fieldnames(optimisationOptions.ranges)));
if ~isempty(missing)
    error('randomFinalDesign:MissingRangeFields', ...
        'optimisationOptions.ranges is missing fields: %s', strjoin(missing, ', '));
end
if isempty(optimisationOptions.ranges.stringer_section_key) || ...
        isempty(optimisationOptions.ranges.spar_cap_sec_no)
    error('randomFinalDesign:EmptySectionRanges', ...
        'Stringer and spar-cap candidate lists must not be empty.');
end
end

function [frontSpar_xc, rearSpar_xc] = randomValidSparPair(ranges)
maxAttempts = 1000;
for iAttempt = 1:maxAttempts
    frontSpar_xc = randomChoice(ranges.frontSpar_xc);
    rearSpar_xc = randomChoice(ranges.rearSpar_xc);
    if rearSpar_xc - frontSpar_xc >= ranges.minSparSeparation_xc
        return;
    end
end
error('randomFinalDesign:SparLocationGenerationFailed', ...
    'Failed to generate valid spar locations after %d attempts.', maxAttempts);
end

function value = randomChoice(values)
idx = randi(numel(values));
value = values(idx);
end

function secNo = stringerSecNoFromKey(designOptions, sectionKey)
matches = string(designOptions.stringers.section_key) == string(sectionKey);
if nnz(matches) ~= 1
    error('randomFinalDesign:StringerKeyMissing', ...
        'Expected exactly one stringer option with section_key = %s; found %d.', ...
        string(sectionKey), nnz(matches));
end
secNo = designOptions.stringers.sec_no(matches);
end

function validateGeneratedDesign(design, designOptions, ranges)
if design.box.rearSpar_xc - design.box.frontSpar_xc < ranges.minSparSeparation_xc
    error('randomFinalDesign:InvalidSparSeparation', ...
        'Generated spar separation is less than %.3f c.', ...
        ranges.minSparSeparation_xc);
end
if ~any(string(designOptions.stringers.section_key) == string(design.stringers.section_key))
    error('randomFinalDesign:StringerNotInOptions', ...
        'Generated stringer section_key is not present in designOptions.');
end
if ~any(designOptions.sparCaps.sec_no == design.sparCaps.sec_no)
    error('randomFinalDesign:SparCapNotInOptions', ...
        'Generated spar cap sec_no is not present in designOptions.');
end
counts = [design.stringers.numTop, design.stringers.numBottom, ...
    design.stringers.numNoseTop, design.stringers.numNoseBottom];
if any(~isfinite(counts)) || any(counts < 1) || any(mod(counts, 1) ~= 0)
    error('randomFinalDesign:InvalidStringerCounts', ...
        'All stringer counts must be positive finite integers.');
end
thicknesses = [design.thickness.t1_frontSpar_in, ...
    design.thickness.t2_rearSpar_in, design.thickness.t3_upperMainSkin_in, ...
    design.thickness.t4_lowerMainSkin_in, design.thickness.upperNoseSkin_in, ...
    design.thickness.lowerNoseSkin_in, design.thickness.t5_noseSkin_in];
if any(~isfinite(thicknesses)) || any(thicknesses <= 0)
    error('randomFinalDesign:InvalidThicknesses', ...
        'All generated thicknesses must be positive finite values.');
end
end
