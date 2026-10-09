function weight = computeWeight(wing, section)
% computeWeight Estimate wing-box structural weight.
%
% Inputs:
%   wing    - Geometry/material struct from setGeometry.
%   section - Struct from computeSectionUnitLoads.
%
% Outputs:
%   weight - Struct containing component and total weights for the current
%            wing-box design.
%
% Scope:
%   Structural weight for panels/webs, spar caps, stringers, and ribs.
%
% Excluded:
%   detailed rib stress sizing, fasteners, joints, doublers, paint, fuel,
%   access panels, manufacturing allowances, buckling, divergence,
%   optimisation, and final failure margins.

validateInputs(wing, section);

semispan_in = wing.planform.semispan_in;
skinWebDensity = wing.material.skinWeb.rho_lbin3;
sparCapStringerDensity = wing.material.sparCapStringer.rho_lbin3;

physicalPanels = section.unitShear.stage3C.unitTotals.physicalPanels;
panelComponents = repmat(makePanelWeightEntry(), numel(physicalPanels), 1);
for iPanel = 1:numel(physicalPanels)
    panel = physicalPanels(iPanel);
    area_in2 = panel.ds_in * panel.thickness_in;
    volume_in3 = area_in2 * semispan_in;
    panelComponents(iPanel) = makePanelWeightEntry(panel, area_in2, ...
        volume_in3, skinWebDensity, volume_in3 * skinWebDensity);
end

sparCaps = section.sparCaps(:);
sparCapComponents = repmat(makeBoomWeightEntry(), numel(sparCaps), 1);
for iCap = 1:numel(sparCaps)
    cap = sparCaps(iCap);
    volume_in3 = cap.area_in2 * semispan_in;
    sparCapComponents(iCap) = makeBoomWeightEntry(cap, iCap, ...
        volume_in3, sparCapStringerDensity, volume_in3 * sparCapStringerDensity, ...
        "sparCap");
end

stringers = collectStringers(section);
stringerComponents = repmat(makeBoomWeightEntry(), numel(stringers), 1);
for iStringer = 1:numel(stringers)
    stringer = stringers(iStringer);
    volume_in3 = stringer.area_in2 * semispan_in;
    stringerComponents(iStringer) = makeBoomWeightEntry(stringer, iStringer, ...
        volume_in3, sparCapStringerDensity, volume_in3 * sparCapStringerDensity, ...
        "stringer");
end

panelWeight_lb = sum([panelComponents.weight_lb]);
sparCapWeight_lb = sum([sparCapComponents.weight_lb]);
stringerWeight_lb = sum([stringerComponents.weight_lb]);
totalExcludingRibsAndMisc_lb = panelWeight_lb + sparCapWeight_lb + stringerWeight_lb;
ribComponents = makeRibWeightEntry(wing.ribs);
totalIncludingRibsExcludingMisc_fullspan_lb = ...
    2.0 * totalExcludingRibsAndMisc_lb + ribComponents.fullspan_weight_lb;
totalIncludingRibsExcludingMisc_semispan_lb = ...
    totalExcludingRibsAndMisc_lb + ribComponents.semispan_weight_lb;

weight = struct();
weight.components.panels = panelComponents;
weight.components.sparCaps = sparCapComponents;
weight.components.stringers = stringerComponents;
weight.components.ribs = ribComponents;
weight.semispan.panel_lb = panelWeight_lb;
weight.semispan.sparCap_lb = sparCapWeight_lb;
weight.semispan.stringer_lb = stringerWeight_lb;
weight.semispan.rib_lb = ribComponents.semispan_weight_lb;
weight.semispan.total_excluding_ribs_and_misc_lb = totalExcludingRibsAndMisc_lb;
weight.semispan.total_including_ribs_excluding_misc_lb = ...
    totalIncludingRibsExcludingMisc_semispan_lb;
weight.fullspan.total_excluding_ribs_and_misc_lb = 2.0 * totalExcludingRibsAndMisc_lb;
weight.fullspan.rib_lb = ribComponents.fullspan_weight_lb;
weight.fullspan.total_including_ribs_excluding_misc_lb = ...
    totalIncludingRibsExcludingMisc_fullspan_lb;
weight.materials.skinWebDensity_lb_in3 = skinWebDensity;
weight.materials.sparCapStringerDensity_lb_in3 = sparCapStringerDensity;
weight.materials.skinWeb = wing.material.skinWeb;
weight.materials.sparCapStringer = wing.material.sparCapStringer;
weight.panelCount = numel(panelComponents);
weight.sparCapCount = numel(sparCapComponents);
weight.stringerCount = numel(stringerComponents);
weight.units.weight = "lb";
weight.units.length = "in";
weight.units.volume = "in^3";
weight.units.density = "lb/in^3";
weight.method = [
    "panel_volume_in3 = ds_in * thickness_in * semispan_in"
    "boom_volume_in3 = area_in2 * semispan_in"
    "fullspan_total_excluding_ribs_and_misc_lb = 2 * semispan_total_excluding_ribs_and_misc_lb"
    "constant_per_rib fullspan rib weight = count_fullspan * weightPerRib_lb"
    "geometric_airfoil_area rib weight = count_fullspan * grossArea_in2 * solidAreaFactor * thickness_in * density_lb_in3"
    ];
weight.assumptions = [
    "Structural weight for the current wing-box model."
    "Constant cross-section over the rectangular semispan."
    "Physical panel list is used to avoid double-counting the shared front spar."
    "Spar caps and stringers are treated as constant-area booms running over the semispan."
    "Full wing weight assumes symmetry."
    "Rib station count includes root and tip on each semispan; full wing counts the root rib once."
    "Rib weight uses the selected rib weight model from wing.ribs."
    ];
weight.warnings = [
    "Rib weight is included using the selected geometric or constant-per-rib model."
    "Detailed rib stress and buckling checks are not included."
    "Fasteners, joints, doublers, paint, fuel, access panels, and manufacturing allowances are excluded."
    "Material densities are read from the verified project material table."
    "Miscellaneous non-structural items are excluded from this weight."
    ];
weight.status = "weight_including_ribs_excluding_misc";
weight.exclusions = [
    "geometry-based rib sizing"
    "fasteners"
    "joints"
    "doublers"
    "paint"
    "fuel"
    "access panels"
    "manufacturing allowances"
    "buckling"
    "divergence"
    "optimisation"
    "final failure margins"
    ];

validateWeight(weight);
end

function entry = makeRibWeightEntry(ribs)
entry = struct();
entry.weightModel = string(ribs.weightModel);
entry.weightPerRib_lb = ribs.weightPerRib_lb;
entry.count_semispan = ribs.count_semispan;
entry.count_fullspan = ribs.count_fullspan;
entry.stations_semispan_in = ribs.stations_semispan_in;
entry.countingConvention = string(ribs.countingConvention);
entry.notes = string(ribs.notes);
entry.rib_weight_single_lb = NaN;
entry.grossArea_in2 = NaN;
entry.solidAreaFactor = NaN;
entry.thickness_in = NaN;
entry.density_lb_in3 = NaN;
entry.material = struct();

switch entry.weightModel
    case "constant_per_rib"
        entry.rib_weight_single_lb = entry.weightPerRib_lb;
        entry.fullspan_weight_lb = entry.count_fullspan * entry.weightPerRib_lb;
        entry.status = "constant_per_rib_weight_included";
    case "geometric_airfoil_area"
        requiredFields = ["grossArea_in2", "solidAreaFactor", ...
            "thickness_in", "rib_density_lb_in3", "material"];
        ribFields = string(fieldnames(ribs));
        missing = setdiff(requiredFields, ribFields);
        if ~isempty(missing)
            error('computeWeight:MissingGeometricRibFields', ...
                'wing.ribs is missing geometric rib fields: %s', strjoin(missing, ', '));
        end
        if ~(isfinite(ribs.grossArea_in2) && ribs.grossArea_in2 > 0 && ...
                isfinite(ribs.solidAreaFactor) && ribs.solidAreaFactor > 0 && ...
                ribs.solidAreaFactor <= 1 && isfinite(ribs.thickness_in) && ...
                ribs.thickness_in > 0 && isfinite(ribs.rib_density_lb_in3) && ...
                ribs.rib_density_lb_in3 > 0)
            error('computeWeight:InvalidGeometricRibFields', ...
                'Geometric rib fields must be positive finite values, with solidAreaFactor <= 1.');
        end
        entry.grossArea_in2 = ribs.grossArea_in2;
        entry.solidAreaFactor = ribs.solidAreaFactor;
        entry.thickness_in = ribs.thickness_in;
        entry.density_lb_in3 = ribs.rib_density_lb_in3;
        entry.material = ribs.material;
        entry.rib_weight_single_lb = entry.grossArea_in2 * ...
            entry.solidAreaFactor * entry.thickness_in * entry.density_lb_in3;
        entry.fullspan_weight_lb = entry.count_fullspan * entry.rib_weight_single_lb;
        entry.status = "geometric_airfoil_area_rib_weight_included";
    otherwise
        error('computeWeight:UnsupportedRibWeightModel', ...
            'Supported rib weight models are constant_per_rib and geometric_airfoil_area.');
end
entry.semispan_weight_lb = entry.fullspan_weight_lb / 2.0;
end

function validateInputs(wing, section)
if nargin < 1 || isempty(wing) || ~isstruct(wing)
    error('computeWeight:MissingWing', 'A wing geometry/material struct is required.');
end
if nargin < 2 || isempty(section) || ~isstruct(section)
    error('computeWeight:MissingSection', 'A section struct from computeSectionUnitLoads is required.');
end

if ~isfield(wing, 'planform') || ~isfield(wing.planform, 'semispan_in') || ...
        ~(isfinite(wing.planform.semispan_in) && wing.planform.semispan_in > 0)
    error('computeWeight:InvalidSemispan', 'wing.planform.semispan_in must be positive and finite.');
end
if ~isfield(wing, 'ribs') || ~isfield(wing.ribs, 'weightModel') || ...
        ~isfield(wing.ribs, 'weightPerRib_lb') || ...
        ~isfield(wing.ribs, 'stations_semispan_in') || ...
        ~isfield(wing.ribs, 'count_semispan') || ...
        ~isfield(wing.ribs, 'count_fullspan')
    error('computeWeight:MissingRibFields', ...
        'wing.ribs must include Stage F2 constant-per-rib weight fields.');
end
if string(wing.ribs.weightModel) ~= "constant_per_rib" && ...
        string(wing.ribs.weightModel) ~= "geometric_airfoil_area"
    error('computeWeight:UnsupportedRibWeightModel', ...
        'Supported wing.ribs.weightModel values are "constant_per_rib" and "geometric_airfoil_area".');
end
if ~(isfinite(wing.ribs.weightPerRib_lb) && wing.ribs.weightPerRib_lb > 0)
    error('computeWeight:InvalidRibWeightPerRib', ...
        'wing.ribs.weightPerRib_lb must be positive and finite.');
end
if any(~isfinite(wing.ribs.stations_semispan_in)) || ...
        any(diff(wing.ribs.stations_semispan_in) <= 0)
    error('computeWeight:InvalidRibStations', ...
        'wing.ribs.stations_semispan_in must be finite and strictly increasing.');
end
if ~(isfinite(wing.ribs.count_semispan) && wing.ribs.count_semispan > 0 && ...
        mod(wing.ribs.count_semispan, 1) == 0 && ...
        isfinite(wing.ribs.count_fullspan) && wing.ribs.count_fullspan > 0 && ...
        mod(wing.ribs.count_fullspan, 1) == 0)
    error('computeWeight:InvalidRibCounts', ...
        'wing.ribs count_semispan and count_fullspan must be positive finite integers.');
end
if wing.ribs.count_fullspan ~= 2 * wing.ribs.count_semispan - 1
    error('computeWeight:InvalidFullspanRibCount', ...
        'Expected wing.ribs.count_fullspan = 2 * count_semispan - 1.');
end
if ~isfield(wing, 'material') || ~isfield(wing.material, 'skinWeb') || ...
        ~isfield(wing.material, 'sparCapStringer')
    error('computeWeight:MissingMaterialGroups', ...
        'wing.material.skinWeb and wing.material.sparCapStringer are required.');
end
if ~isfield(wing.material.skinWeb, 'rho_lbin3') || ...
        ~(isfinite(wing.material.skinWeb.rho_lbin3) && wing.material.skinWeb.rho_lbin3 > 0)
    error('computeWeight:InvalidSkinWebDensity', ...
        'wing.material.skinWeb.rho_lbin3 must be positive and finite.');
end
if ~isfield(wing.material.sparCapStringer, 'rho_lbin3') || ...
        ~(isfinite(wing.material.sparCapStringer.rho_lbin3) && wing.material.sparCapStringer.rho_lbin3 > 0)
    error('computeWeight:InvalidSparCapStringerDensity', ...
        'wing.material.sparCapStringer.rho_lbin3 must be positive and finite.');
end

if ~isfield(section, 'unitShear') || ~isfield(section.unitShear, 'stage3C') || ...
        ~isfield(section.unitShear.stage3C, 'unitTotals') || ...
        ~isfield(section.unitShear.stage3C.unitTotals, 'physicalPanels')
    error('computeWeight:MissingPhysicalPanels', ...
        'section.unitShear.stage3C.unitTotals.physicalPanels is required.');
end
physicalPanels = section.unitShear.stage3C.unitTotals.physicalPanels;
if isempty(physicalPanels)
    error('computeWeight:EmptyPhysicalPanels', 'Physical panel list must not be empty.');
end
requiredPanelFields = ["name", "group", "ds_in", "thickness_in"];
panelFields = string(fieldnames(physicalPanels));
missingPanelFields = setdiff(requiredPanelFields, panelFields);
if ~isempty(missingPanelFields)
    error('computeWeight:MissingPanelFields', ...
        'Physical panel list is missing required fields: %s', ...
        strjoin(missingPanelFields, ', '));
end
if any(~isfinite([physicalPanels.ds_in])) || any([physicalPanels.ds_in] <= 0) || ...
        any(~isfinite([physicalPanels.thickness_in])) || any([physicalPanels.thickness_in] <= 0)
    error('computeWeight:InvalidPanelGeometry', ...
        'Physical panel ds_in and thickness_in values must be positive and finite.');
end
if nnz([physicalPanels.group] == "frontSparSharedPhysical") ~= 1
    error('computeWeight:SharedFrontSparCount', ...
        'Physical panel list must contain exactly one shared front spar panel.');
end
if any([physicalPanels.group] == "frontSpar") || any([physicalPanels.group] == "frontSparCell2")
    error('computeWeight:SharedFrontSparDoubleCountRisk', ...
        'Physical panel list must not contain duplicate front spar loop panels.');
end

if ~isfield(section, 'sparCaps') || isempty(section.sparCaps)
    error('computeWeight:MissingSparCaps', 'section.sparCaps is required.');
end
if any(~isfinite([section.sparCaps.area_in2])) || any([section.sparCaps.area_in2] <= 0)
    error('computeWeight:InvalidSparCapAreas', 'Spar cap areas must be positive and finite.');
end

stringers = collectStringers(section);
if isempty(stringers)
    error('computeWeight:MissingStringers', 'At least one stringer is required for this baseline.');
end
if any(~isfinite([stringers.area_in2])) || any([stringers.area_in2] <= 0)
    error('computeWeight:InvalidStringerAreas', 'Stringer areas must be positive and finite.');
end
end

function stringers = collectStringers(section)
if ~isfield(section, 'stringers') || ...
        ~isfield(section.stringers, 'top') || ...
        ~isfield(section.stringers, 'bottom') || ...
        ~isfield(section.stringers, 'noseTop') || ...
        ~isfield(section.stringers, 'noseBottom')
    error('computeWeight:MissingStringerGroups', ...
        'section.stringers must include top, bottom, noseTop, and noseBottom.');
end

stringers = [
    section.stringers.top(:)
    section.stringers.bottom(:)
    section.stringers.noseTop(:)
    section.stringers.noseBottom(:)
    ];
end

function entry = makePanelWeightEntry(panel, area_in2, volume_in3, density_lb_in3, weight_lb)
if nargin == 0
    panel = makePanelPlaceholder();
    area_in2 = NaN;
    volume_in3 = NaN;
    density_lb_in3 = NaN;
    weight_lb = NaN;
end

entry = struct();
entry.name = panel.name;
entry.group = panel.group;
entry.ds_in = panel.ds_in;
entry.thickness_in = panel.thickness_in;
entry.area_in2 = area_in2;
entry.semispan_volume_in3 = volume_in3;
entry.density_lb_in3 = density_lb_in3;
entry.weight_lb = weight_lb;
end

function entry = makeBoomWeightEntry(boom, index, volume_in3, density_lb_in3, weight_lb, componentType)
if nargin == 0
    boom = makeBoomPlaceholder();
    index = NaN;
    volume_in3 = NaN;
    density_lb_in3 = NaN;
    weight_lb = NaN;
    componentType = "";
end

entry = struct();
entry.name = boom.name;
entry.component_type = string(componentType);
entry.index = index;
entry.x_in = boom.x_in;
entry.z_in = boom.z_in;
entry.area_in2 = boom.area_in2;
entry.semispan_volume_in3 = volume_in3;
entry.density_lb_in3 = density_lb_in3;
entry.weight_lb = weight_lb;
end

function panel = makePanelPlaceholder()
panel = struct();
panel.name = "";
panel.group = "";
panel.ds_in = NaN;
panel.thickness_in = NaN;
end

function boom = makeBoomPlaceholder()
boom = struct();
boom.name = "";
boom.x_in = NaN;
boom.z_in = NaN;
boom.area_in2 = NaN;
end

function validateWeight(weight)
values = [
    weight.semispan.panel_lb
    weight.semispan.sparCap_lb
    weight.semispan.stringer_lb
    weight.semispan.rib_lb
    weight.semispan.total_excluding_ribs_and_misc_lb
    weight.semispan.total_including_ribs_excluding_misc_lb
    weight.fullspan.rib_lb
    weight.fullspan.total_excluding_ribs_and_misc_lb
    weight.fullspan.total_including_ribs_excluding_misc_lb
    ];
if any(~isfinite(values)) || any(values <= 0)
    error('computeWeight:InvalidComputedWeight', ...
        'Computed weights must be positive and finite.');
end
expectedFullspan = weight.fullspan.total_excluding_ribs_and_misc_lb + ...
    weight.fullspan.rib_lb;
if abs(weight.fullspan.total_including_ribs_excluding_misc_lb - expectedFullspan) > 1e-9
    error('computeWeight:InvalidIncludingRibsTotal', ...
        'Fullspan total including ribs must equal excluding-ribs total plus fullspan rib weight.');
end
end

