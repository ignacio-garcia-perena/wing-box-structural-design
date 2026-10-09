function wing = setGeometry(design, designOptions, analysisOptions)
% setGeometry Define the baseline wing-box geometry.
%
% Inputs:
%   design        - Optional design struct. If omitted, baselineDesign
%                   supplies the current deterministic Bruhn baseline.
%   designOptions - Optional design option tables from loadDesignOptions.
%                   If omitted, they are loaded from data/*.csv.
%   analysisOptions - Optional analysis settings. If
%                     suppressPracticalGeometryWarnings = true, broad
%                     practical spar-location warnings are muted while hard
%                     geometry validation remains active.
%
% Outputs:
%   wing - Struct containing geometry, material groups, and layout settings
%          for a two-cell NACA 2412 wing box.

if nargin < 1 || isempty(design)
    design = baselineDesign();
end
if nargin < 2 || isempty(designOptions)
    designOptions = loadDesignOptions();
end
if nargin < 3 || isempty(analysisOptions)
    analysisOptions = struct();
end

selected = selectDesignRows(design, designOptions, analysisOptions);

wing = struct();

wing.airfoil.name = "NACA 2412";

wing.units.length = "in";
wing.units.area = "in^2";
wing.units.stress = "psi";
wing.units.density = "lb/in^3";

wing.planform.chord_in = 5.0 * 12.0;       % Project geometry.
wing.planform.semispan_in = 17.0 * 12.0;   % Project geometry.

% Stage: spar locations are design-driven. baselineDesign()
% preserves the original deterministic baseline values.
wing.box.frontSpar_xc = design.box.frontSpar_xc;
wing.box.rearSpar_xc = design.box.rearSpar_xc;
wing.box.frontSpar_x_in = wing.box.frontSpar_xc * wing.planform.chord_in;
wing.box.rearSpar_x_in = wing.box.rearSpar_xc * wing.planform.chord_in;

% Note: NACA 2412 surface coordinates are generated downstream by the
% section-analysis helpers used by computeSectionUnitLoads.
wing.box.upperSkin = [];
wing.box.lowerSkin = [];
wing.box.spars = [];
wing.box.cells = [];

% Stage: spar-cap section comes from designOptions, but the same section
% is still applied to all four spar caps.
sparCapArea_in2 = selected.sparCap.area_in2;
wing.sparCaps.A1_in2 = sparCapArea_in2;
wing.sparCaps.A2_in2 = sparCapArea_in2;
wing.sparCaps.A3_in2 = sparCapArea_in2;
wing.sparCaps.A4_in2 = sparCapArea_in2;
wing.sparCaps.areas_in2 = [wing.sparCaps.A1_in2, wing.sparCaps.A2_in2, ...
    wing.sparCaps.A3_in2, wing.sparCaps.A4_in2];
wing.sparCaps.section = selected.sparCap;
wing.sparCaps.section.warning = ...
    "T-section axis labels should be confirmed before final buckling sizing.";

% Stage: thickness values are design-driven. baselineDesign()
% preserves the original deterministic baseline values.
wing.thickness.t1_frontSpar_in = design.thickness.t1_frontSpar_in;
wing.thickness.t2_rearSpar_in = design.thickness.t2_rearSpar_in;
wing.thickness.t3_upperMainSkin_in = design.thickness.t3_upperMainSkin_in;
wing.thickness.t4_lowerMainSkin_in = design.thickness.t4_lowerMainSkin_in;
wing.thickness.t5_noseSkin_in = design.thickness.t5_noseSkin_in;
wing.thickness.upperNoseSkin_in = design.thickness.upperNoseSkin_in;
wing.thickness.lowerNoseSkin_in = design.thickness.lowerNoseSkin_in;

% Stage: stringer counts are design-driven. baselineDesign()
% preserves the original deterministic baseline counts.
wing.stringers.numTop = design.stringers.numTop;
wing.stringers.numBottom = design.stringers.numBottom;
wing.stringers.numNoseTop = design.stringers.numNoseTop;
wing.stringers.numNoseBottom = design.stringers.numNoseBottom;
% Stage: stringer section comes from designOptions, but the same section
% is still applied to all stringer groups.
stringerArea_in2 = selected.stringer.area_in2;
wing.stringers.areaTop_in2 = stringerArea_in2;
wing.stringers.areaBottom_in2 = stringerArea_in2;
wing.stringers.areaNoseTop_in2 = stringerArea_in2;
wing.stringers.areaNoseBottom_in2 = stringerArea_in2;
wing.stringers.placementRule = ...
    "Even x/c spacing is used to initialise stringer placement before geometry refinement.";
wing.stringers.type = selected.stringer.section_type;
wing.stringers.geometry = selected.stringer;

wing.ribs.spacing_in = design.ribs.ribSpacing_in;
wing.ribs.weightModel = string(design.ribs.weightModel);
wing.ribs.weightPerRib_lb = design.ribs.weightPerRib_lb;
wing.ribs.countingConvention = string(design.ribs.countingConvention);
wing.ribs.notes = string(design.ribs.notes);
wing.ribs.material = selected.ribMaterial;
wing.ribs.rib_density_lb_in3 = selected.ribMaterial.rho_lbin3;
if isfield(design.ribs, 'airfoil')
    wing.ribs.airfoil = string(design.ribs.airfoil);
else
    wing.ribs.airfoil = wing.airfoil.name;
end
wing.ribs.normalizedArea_A_over_c2 = getOptionalRibScalar( ...
    design.ribs, "normalizedArea_A_over_c2", NaN);
wing.ribs.chord_in = getOptionalRibScalar(design.ribs, "chord_in", ...
    wing.planform.chord_in);
wing.ribs.grossArea_in2 = getOptionalRibScalar(design.ribs, ...
    "grossArea_in2", NaN);
wing.ribs.solidAreaFactor = getOptionalRibScalar(design.ribs, ...
    "solidAreaFactor", NaN);
wing.ribs.thickness_in = getOptionalRibScalar(design.ribs, ...
    "thickness_in", NaN);
wing.ribs.stations_semispan_in = makeRibStations( ...
    wing.planform.semispan_in, wing.ribs.spacing_in);
wing.ribs.count_semispan = numel(wing.ribs.stations_semispan_in);
wing.ribs.count_fullspan = 2 * wing.ribs.count_semispan - 1;
wing.columnBuckling.K_column = design.columnBuckling.K_column;

% Material properties are read from the verified project material table.
wing.material.skinWeb = selected.skinWebMaterial;
wing.material.sparCapStringer = selected.sparCapStringerMaterial;

% The final model uses separate material groups for skins/webs/ribs and
% spar caps/stringers. The grouped material structs above are used by the
% analysis modules.
wing.material.name = "Mixed baseline: 2024-T3 skin/web, 7075-T6 spar cap/stringer";
wing.material.temper = "";
wing.material.source = "Verified project material table; see data/material_properties.csv.";
wing.material.allowable_basis = "Component-specific allowables are stored in wing.material.skinWeb and wing.material.sparCapStringer.";
wing.material.E_psi = NaN;
wing.material.G_psi = NaN;
wing.material.nu = NaN;
wing.material.rho_lbin3 = NaN;
wing.material.sigma_yield_psi = NaN;
wing.material.sigma_ult_psi = NaN;
wing.material.isPlaceholder = false;

wing.analysis.status = "design_driven_geometry_model";
wing.analysis.design_label = string(design.label);
wing.analysis.source = "Bruhn section/material baseline with design-driven section, material, rib spacing, and K selections.";
wing.analysis.warning = ...
    "Final optimisation varies spar locations, thicknesses, stringer counts, rib spacing, and section selections within the configured design ranges.";
wing.analysis.designWarnings = [
    "Material values with confidence other than HIGH should be source-checked."
    "T-section axis labels should be confirmed before final buckling sizing."
    "Rib spacing and K_column are configured design inputs."
    "Stage O3 applies thickness variables from design.thickness."
    "Stage O4 applies stringer count variables from design.stringers."
    "Stage O5 applies spar location variables from design.box."
    ];

validateGeometry(wing);
end

function validateGeometry(wing)
% validateGeometry Check baseline geometry fields without requiring final
% rib spacing or material values.

if ~(isfinite(wing.planform.chord_in) && wing.planform.chord_in > 0)
    error('setGeometry:InvalidChord', 'wing.planform.chord_in must be positive and finite.');
end

if ~(isfinite(wing.box.frontSpar_xc) && isfinite(wing.box.rearSpar_xc) && ...
        wing.box.frontSpar_xc > 0 && wing.box.frontSpar_xc < wing.box.rearSpar_xc && ...
        wing.box.rearSpar_xc < 1)
    error('setGeometry:InvalidSparLocations', ...
        'Expected 0 < frontSpar_xc < rearSpar_xc < 1.');
end

if any(~isfinite(wing.sparCaps.areas_in2)) || any(wing.sparCaps.areas_in2 <= 0)
    error('setGeometry:InvalidSparCapAreas', 'All spar cap areas must be positive and finite.');
end

thicknessValues = [wing.thickness.t1_frontSpar_in, wing.thickness.t2_rearSpar_in, ...
    wing.thickness.t3_upperMainSkin_in, wing.thickness.t4_lowerMainSkin_in, ...
    wing.thickness.t5_noseSkin_in, wing.thickness.upperNoseSkin_in, ...
    wing.thickness.lowerNoseSkin_in];
if any(~isfinite(thicknessValues)) || any(thicknessValues <= 0)
    error('setGeometry:InvalidThicknesses', 'All thickness values must be positive and finite.');
end

stringerCounts = [wing.stringers.numTop, wing.stringers.numBottom, ...
    wing.stringers.numNoseTop, wing.stringers.numNoseBottom];
if any(~isfinite(stringerCounts)) || any(stringerCounts < 0) || any(mod(stringerCounts, 1) ~= 0)
    error('setGeometry:InvalidStringerCounts', ...
        'Stringer counts must be nonnegative finite integers.');
end

stringerAreas = [wing.stringers.areaTop_in2, wing.stringers.areaBottom_in2, ...
    wing.stringers.areaNoseTop_in2, wing.stringers.areaNoseBottom_in2];
if any(~isfinite(stringerAreas)) || any(stringerAreas <= 0)
    error('setGeometry:InvalidStringerAreas', 'Stringer areas must be positive and finite.');
end

if ~isfield(wing, 'ribs') || ~isfield(wing.ribs, 'spacing_in') || ...
        ~(isfinite(wing.ribs.spacing_in) && wing.ribs.spacing_in > 0)
    error('setGeometry:InvalidRibSpacing', ...
        'wing.ribs.spacing_in must be positive and finite when design input is used.');
end

if ~isfield(wing, 'columnBuckling') || ...
        ~isfield(wing.columnBuckling, 'K_column') || ...
        ~(isfinite(wing.columnBuckling.K_column) && wing.columnBuckling.K_column > 0)
    error('setGeometry:InvalidKColumn', ...
        'wing.columnBuckling.K_column must be positive and finite.');
end
end

function selected = selectDesignRows(design, designOptions, analysisOptions)
if nargin < 3
    analysisOptions = struct();
end
validateDesignInputs(design, designOptions);
validateDesignBox(design.box, analysisOptions);
validateDesignRibs(design.ribs);
validateDesignThicknesses(design.thickness);
validateDesignStringerCounts(design.stringers);

stringerRow = lookupStringerSection(designOptions.stringers, design.stringers);
sparCapRow = lookupBySecNo(designOptions.sparCaps, ...
    design.sparCaps.sec_no, 'spar cap');
skinWebRow = lookupMaterial(designOptions.materials, ...
    design.material.skinWeb, 'skin/web');
sparCapStringerRow = lookupMaterial(designOptions.materials, ...
    design.material.sparCapStringer, 'spar cap/stringer');
if isfield(design.ribs, 'material')
    ribMaterialRow = lookupMaterial(designOptions.materials, ...
        design.ribs.material, 'rib');
else
    ribMaterialRow = skinWebRow;
end

selected = struct();
selected.stringer = makeStringerSection(stringerRow);
selected.sparCap = makeSparCapSection(sparCapRow);
selected.skinWebMaterial = makeMaterialGroup(skinWebRow);
selected.sparCapStringerMaterial = makeMaterialGroup(sparCapStringerRow);
selected.ribMaterial = makeMaterialGroup(ribMaterialRow);

validateSelectedRows(selected);
warnForMaterialConfidence(selected.skinWebMaterial, 'skin/web');
warnForMaterialConfidence(selected.sparCapStringerMaterial, 'spar cap/stringer');
warnForMaterialConfidence(selected.ribMaterial, 'rib');
end

function validateDesignBox(box, analysisOptions)
if nargin < 2
    analysisOptions = struct();
end
requiredFields = ["frontSpar_xc", "rearSpar_xc"];
for iField = 1:numel(requiredFields)
    fieldName = requiredFields(iField);
    if ~isfield(box, fieldName)
        error('setGeometry:MissingBoxField', ...
            'design.box.%s is required.', fieldName);
    end

    value = box.(fieldName);
    if ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('setGeometry:InvalidBoxValue', ...
            'design.box.%s must be a finite numeric scalar.', fieldName);
    end
end

frontSpar_xc = box.frontSpar_xc;
rearSpar_xc = box.rearSpar_xc;
if ~(frontSpar_xc > 0 && frontSpar_xc < rearSpar_xc && rearSpar_xc < 1)
    error('setGeometry:InvalidDesignSparLocations', ...
        'Expected 0 < design.box.frontSpar_xc < design.box.rearSpar_xc < 1.');
end
if rearSpar_xc - frontSpar_xc < 0.20
    error('setGeometry:InsufficientSparSpacing', ...
        'Expected design.box.rearSpar_xc - design.box.frontSpar_xc >= 0.20.');
end

if shouldShowPracticalGeometryWarnings(analysisOptions) && ...
        (frontSpar_xc < 0.10 || frontSpar_xc > 0.30)
    warning('setGeometry:FrontSparOutsidePracticalRange', ...
        'design.box.frontSpar_xc = %.4f is outside the broad practical range [0.10, 0.30].', ...
        frontSpar_xc);
end
if shouldShowPracticalGeometryWarnings(analysisOptions) && ...
        (rearSpar_xc < 0.60 || rearSpar_xc > 0.85)
    warning('setGeometry:RearSparOutsidePracticalRange', ...
        'design.box.rearSpar_xc = %.4f is outside the broad practical range [0.60, 0.85].', ...
        rearSpar_xc);
end
end

function tf = shouldShowPracticalGeometryWarnings(analysisOptions)
tf = true;
if isstruct(analysisOptions) && ...
        isfield(analysisOptions, 'suppressPracticalGeometryWarnings') && ...
        logical(analysisOptions.suppressPracticalGeometryWarnings)
    tf = false;
end
end

function validateDesignStringerCounts(stringers)
requiredFields = [
    "numTop"
    "numBottom"
    "numNoseTop"
    "numNoseBottom"
    ];

for iField = 1:numel(requiredFields)
    fieldName = requiredFields(iField);
    if ~isfield(stringers, fieldName)
        error('setGeometry:MissingStringerCountField', ...
            'design.stringers.%s is required.', fieldName);
    end

    value = stringers.(fieldName);
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
            value >= 0 && mod(value, 1) == 0)
        error('setGeometry:InvalidStringerCount', ...
            'design.stringers.%s must be a nonnegative finite integer scalar.', ...
            fieldName);
    end
end

if stringers.numTop < 1
    error('setGeometry:InvalidTopStringerCount', ...
        'design.stringers.numTop must be at least 1.');
end
if stringers.numBottom < 1
    error('setGeometry:InvalidBottomStringerCount', ...
        'design.stringers.numBottom must be at least 1.');
end
end

function validateDesignThicknesses(thickness)
requiredFields = [
    "t1_frontSpar_in"
    "t2_rearSpar_in"
    "t3_upperMainSkin_in"
    "t4_lowerMainSkin_in"
    "t5_noseSkin_in"
    "upperNoseSkin_in"
    "lowerNoseSkin_in"
    ];

for iField = 1:numel(requiredFields)
    fieldName = requiredFields(iField);
    if ~isfield(thickness, fieldName)
        error('setGeometry:MissingThicknessField', ...
            'design.thickness.%s is required.', fieldName);
    end

    value = thickness.(fieldName);
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
        error('setGeometry:InvalidThicknessValue', ...
            'design.thickness.%s must be a positive finite numeric scalar.', ...
            fieldName);
    end

    if value < 0.005 || value > 0.250
        warning('setGeometry:ThicknessOutsidePracticalRange', ...
            'design.thickness.%s = %.4f in is outside the broad practical range [0.005, 0.250] in.', ...
            fieldName, value);
    end
end
end

function validateDesignInputs(design, designOptions)
if ~isstruct(design)
    error('setGeometry:InvalidDesign', 'design must be a struct.');
end
if ~isstruct(designOptions)
    error('setGeometry:InvalidDesignOptions', 'designOptions must be a struct.');
end
if ~isfield(designOptions, 'stringers') || ~isfield(designOptions, 'sparCaps') || ...
        ~isfield(designOptions, 'materials')
    error('setGeometry:MissingDesignOptionTables', ...
        'designOptions must include stringers, sparCaps, and materials tables.');
end

requiredDesignFields = [
    isfield(design, 'stringers') && isfield(design.stringers, 'sec_no')
    isfield(design, 'sparCaps') && isfield(design.sparCaps, 'sec_no')
    isfield(design, 'material') && isfield(design.material, 'skinWeb')
    isfield(design, 'material') && isfield(design.material, 'sparCapStringer')
    isfield(design, 'ribs') && isfield(design.ribs, 'ribSpacing_in')
    isfield(design, 'columnBuckling') && isfield(design.columnBuckling, 'K_column')
    ];
if any(~requiredDesignFields)
    error('setGeometry:MissingDesignFields', ...
        'design is missing required Stage O2 fields.');
end

if ~(isnumeric(design.columnBuckling.K_column) && isscalar(design.columnBuckling.K_column) && ...
        isfinite(design.columnBuckling.K_column) && design.columnBuckling.K_column > 0)
    error('setGeometry:InvalidDesignKColumn', ...
        'design.columnBuckling.K_column must be positive and finite.');
end
end

function validateDesignRibs(ribs)
requiredFields = ["ribSpacing_in", "weightModel", "weightPerRib_lb", ...
    "countingConvention", "notes"];
for iField = 1:numel(requiredFields)
    fieldName = requiredFields(iField);
    if ~isfield(ribs, fieldName)
        error('setGeometry:MissingRibField', ...
            'design.ribs.%s is required for final-design rib weight accounting.', ...
            fieldName);
    end
end

if ~(isnumeric(ribs.ribSpacing_in) && isscalar(ribs.ribSpacing_in) && ...
        isfinite(ribs.ribSpacing_in) && ribs.ribSpacing_in > 0)
    error('setGeometry:InvalidDesignRibSpacing', ...
        'design.ribs.ribSpacing_in must be positive and finite.');
end
weightModel = string(ribs.weightModel);
if weightModel ~= "constant_per_rib" && weightModel ~= "geometric_airfoil_area"
    error('setGeometry:UnsupportedRibWeightModel', ...
        'Supported rib weight models are "constant_per_rib" and "geometric_airfoil_area".');
end
if ~(isnumeric(ribs.weightPerRib_lb) && isscalar(ribs.weightPerRib_lb) && ...
        isfinite(ribs.weightPerRib_lb) && ribs.weightPerRib_lb > 0)
    error('setGeometry:InvalidRibWeightPerRib', ...
        'design.ribs.weightPerRib_lb must be a positive finite scalar.');
end

if weightModel == "geometric_airfoil_area"
    requiredGeometricFields = ["material", "airfoil", "normalizedArea_A_over_c2", ...
        "chord_in", "grossArea_in2", "solidAreaFactor", "thickness_in"];
    for iField = 1:numel(requiredGeometricFields)
        fieldName = requiredGeometricFields(iField);
        if ~isfield(ribs, fieldName)
            error('setGeometry:MissingGeometricRibField', ...
                'design.ribs.%s is required for geometric_airfoil_area rib weight.', ...
                fieldName);
        end
    end

    scalarFields = ["normalizedArea_A_over_c2", "chord_in", ...
        "grossArea_in2", "solidAreaFactor", "thickness_in"];
    for iField = 1:numel(scalarFields)
        fieldName = scalarFields(iField);
        value = ribs.(fieldName);
        if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
            error('setGeometry:InvalidGeometricRibValue', ...
                'design.ribs.%s must be a positive finite scalar.', fieldName);
        end
    end
    if ribs.solidAreaFactor > 1
        error('setGeometry:InvalidRibSolidAreaFactor', ...
            'design.ribs.solidAreaFactor must be no greater than 1.');
    end
end
end

function value = getOptionalRibScalar(ribs, fieldName, defaultValue)
fieldName = char(fieldName);
if isfield(ribs, fieldName)
    value = ribs.(fieldName);
else
    value = defaultValue;
end
end

function ribStations_semispan_in = makeRibStations(semispan_in, ribSpacing_in)
ribStations_semispan_in = 0:ribSpacing_in:semispan_in;
if isempty(ribStations_semispan_in) || ribStations_semispan_in(1) ~= 0
    error('setGeometry:InvalidRibStations', ...
        'Rib station generation must include the root station y = 0.');
end
if ribStations_semispan_in(end) < semispan_in
    ribStations_semispan_in(end + 1) = semispan_in;
end
if any(~isfinite(ribStations_semispan_in)) || any(diff(ribStations_semispan_in) <= 0)
    error('setGeometry:InvalidRibStations', ...
        'Generated rib stations must be finite and strictly increasing.');
end
end

function row = lookupBySecNo(optionTable, secNo, label)
if ~(isnumeric(secNo) && isscalar(secNo) && isfinite(secNo))
    error('setGeometry:InvalidSectionNumber', ...
        'Selected %s sec_no must be a finite numeric scalar.', label);
end

matches = optionTable.sec_no == secNo;
if nnz(matches) ~= 1
    error('setGeometry:SectionLookupFailed', ...
        'Expected exactly one %s option with sec_no = %g; found %d.', ...
        label, secNo, nnz(matches));
end
row = optionTable(matches, :);
end

function row = lookupStringerSection(stringerTable, stringerDesign)
if isfield(stringerDesign, 'section_key') && strlength(string(stringerDesign.section_key)) > 0
    sectionKey = string(stringerDesign.section_key);
    matches = string(stringerTable.section_key) == sectionKey;
    if nnz(matches) ~= 1
        error('setGeometry:StringerKeyLookupFailed', ...
            'Expected exactly one stringer option with section_key = %s; found %d.', ...
            sectionKey, nnz(matches));
    end
    row = stringerTable(matches, :);

    if isfield(stringerDesign, 'sec_no') && isnumeric(stringerDesign.sec_no) && ...
            isscalar(stringerDesign.sec_no) && row.sec_no ~= stringerDesign.sec_no
        error('setGeometry:StringerKeySecNoMismatch', ...
            'design.stringers.section_key = %s maps to sec_no = %g, but design.stringers.sec_no = %g.', ...
            sectionKey, row.sec_no, stringerDesign.sec_no);
    end
    return;
end

if ~isfield(stringerDesign, 'sec_no')
    error('setGeometry:MissingStringerSelection', ...
        'design.stringers.section_key or design.stringers.sec_no is required.');
end

if ~(isnumeric(stringerDesign.sec_no) && isscalar(stringerDesign.sec_no) && ...
        isfinite(stringerDesign.sec_no))
    error('setGeometry:InvalidSectionNumber', ...
        'Selected stringer sec_no must be a finite numeric scalar.');
end

matches = stringerTable.sec_no == stringerDesign.sec_no;
if nnz(matches) > 1
    matchingKeys = strjoin(string(stringerTable.section_key(matches)), ', ');
    error('setGeometry:AmbiguousStringerSectionNumber', ...
        ['Stringer sec_no = %g matches multiple source-table rows: %s. ', ...
        'Set design.stringers.section_key to select the intended source explicitly.'], ...
        stringerDesign.sec_no, matchingKeys);
elseif nnz(matches) ~= 1
    error('setGeometry:SectionLookupFailed', ...
        'Expected exactly one stringer option with sec_no = %g; found %d.', ...
        stringerDesign.sec_no, nnz(matches));
end
row = stringerTable(matches, :);
end

function row = lookupMaterial(materials, materialDesign, label)
requiredFields = ["alloy", "temper", "product_form"];
for iField = 1:numel(requiredFields)
    fieldName = requiredFields(iField);
    if ~isfield(materialDesign, fieldName)
        error('setGeometry:MissingMaterialDesignField', ...
            'design.material.%s is missing field %s.', label, fieldName);
    end
end

matches = string(materials.alloy) == string(materialDesign.alloy) & ...
    string(materials.temper) == string(materialDesign.temper) & ...
    string(materials.product_form) == string(materialDesign.product_form);
if nnz(matches) ~= 1
    error('setGeometry:MaterialLookupFailed', ...
        'Expected exactly one %s material row for %s-%s %s; found %d.', ...
        label, string(materialDesign.alloy), string(materialDesign.temper), ...
        string(materialDesign.product_form), nnz(matches));
end
row = materials(matches, :);
end

function section = makeStringerSection(row)
section = struct();
section.section_key = string(row.section_key);
section.sec_no = row.sec_no;
section.section_type = string(row.section_type);
section.A_in = row.A_in;
section.B_in = row.B_in;
section.T_in = row.T_in;
section.R_in = row.R_in;
section.area_in2 = row.area_in2;
section.Ixx_in4 = row.Ixx_in4;
section.Iyy_in4 = row.Iyy_in4;
section.rho_xx_in = row.rho_xx_in;
section.rho_yy_in = row.rho_yy_in;
section.source = string(row.source_table) + " " + string(row.section_key);
section.confidence = string(row.confidence);
section.notes = string(row.notes);
end

function section = makeSparCapSection(row)
section = struct();
section.sec_no = row.sec_no;
section.section_type = string(row.section_type);
section.A_in = row.A_in;
section.B_in = row.B_in;
section.T_in = row.T_in;
section.R_in = row.R_in;
section.area_in2 = row.area_in2;
section.zbar_in = row.zbar_in;
section.Iyy_in4 = row.Iyy_in4;
section.Izz_in4 = row.Izz_in4;
section.rho_xx_in = row.rho_xx_in;
section.rho_yy_in = row.rho_yy_in;
section.source = string(row.source_table) + " sec_no " + string(row.sec_no);
section.confidence = string(row.confidence);
section.notes = string(row.notes);
end

function material = makeMaterialGroup(row)
material = struct();
material.name = string(row.alloy);
material.temper = string(row.temper);
material.product_form = string(row.product_form);
material.E_psi = row.E_psi;
if ismember("Ec_psi", string(row.Properties.VariableNames))
    material.Ec_psi = row.Ec_psi;
else
    material.Ec_psi = row.E_psi;
end
material.G_psi = row.G_psi;
material.nu = row.nu;
material.rho_lbin3 = row.density_lb_in3;
material.Ftu_psi = row.Ftu_psi;
material.Fty_psi = row.Fty_psi;
material.Fcy_psi = row.Fcy_psi;
material.Fsu_psi = row.Fsu_psi;
material.source = string(row.source_file) + ":" + string(row.page_or_image) + ...
    ":" + string(row.table_title);
material.confidence = string(row.confidence);
material.notes = string(row.notes);
end

function validateSelectedRows(selected)
stringerValues = [selected.stringer.A_in, selected.stringer.B_in, ...
    selected.stringer.T_in, selected.stringer.R_in, ...
    selected.stringer.area_in2, selected.stringer.Ixx_in4, ...
    selected.stringer.Iyy_in4, selected.stringer.rho_xx_in, ...
    selected.stringer.rho_yy_in];
if any(~isfinite(stringerValues)) || any(stringerValues <= 0)
    error('setGeometry:InvalidSelectedStringer', ...
        'Selected stringer section values must be positive and finite.');
end

sparCapValues = [selected.sparCap.A_in, selected.sparCap.B_in, ...
    selected.sparCap.T_in, selected.sparCap.R_in, ...
    selected.sparCap.area_in2, selected.sparCap.zbar_in, ...
    selected.sparCap.Iyy_in4, selected.sparCap.Izz_in4, ...
    selected.sparCap.rho_xx_in, selected.sparCap.rho_yy_in];
if any(~isfinite(sparCapValues)) || any(sparCapValues <= 0)
    error('setGeometry:InvalidSelectedSparCap', ...
        'Selected spar-cap section values must be positive and finite.');
end

validateMaterial(selected.skinWebMaterial, 'skin/web');
validateMaterial(selected.sparCapStringerMaterial, 'spar cap/stringer');
end

function validateMaterial(material, label)
values = [material.E_psi, material.Ec_psi, material.G_psi, material.nu, ...
    material.rho_lbin3, material.Ftu_psi, material.Fty_psi, ...
    material.Fcy_psi, material.Fsu_psi];
if any(~isfinite(values)) || any(values <= 0)
    error('setGeometry:InvalidSelectedMaterial', ...
        'Selected %s material values must be positive and finite.', label);
end
end

function warnForMaterialConfidence(material, label)
if upper(string(material.confidence)) ~= "HIGH"
    warning('setGeometry:MaterialConfidence', ...
        'Selected %s material confidence is %s; final material allowables to be confirmed before final report.', ...
        label, string(material.confidence));
end
end


