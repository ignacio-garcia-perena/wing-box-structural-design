function designOptions = loadDesignOptions(projectRoot)
% loadDesignOptions Load Bruhn-derived design option lookup tables.
%
% Inputs:
%   projectRoot - Optional project root folder. If omitted, it is inferred
%                 from this file location.
%
% Outputs:
%   designOptions - Struct containing stringer, spar-cap, and material
%                   option tables plus status and warnings.
%
% Scope:
%   Read, validate, and package design option data only.
%
% Excluded:
%   geometry selection, yield checks, buckling, divergence, weight, and
%   optimisation.

if nargin < 1 || isempty(projectRoot)
    scriptFolder = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(scriptFolder, '..');
end

dataFolder = fullfile(projectRoot, 'data');
stringerFile = fullfile(dataFolder, 'stringer_options.csv');
sparCapFile = fullfile(dataFolder, 'spar_cap_options.csv');
materialFile = fullfile(dataFolder, 'material_properties.csv');

if ~isfile(stringerFile)
    error('loadDesignOptions:MissingStringerOptions', ...
        'Missing stringer options table: %s', stringerFile);
end
if ~isfile(sparCapFile)
    error('loadDesignOptions:MissingSparCapOptions', ...
        'Missing spar-cap options table: %s', sparCapFile);
end
if ~isfile(materialFile)
    error('loadDesignOptions:MissingMaterialOptions', ...
        'Missing material options table: %s', materialFile);
end

stringers = readtable(stringerFile, 'TextType', 'string');
sparCaps = readtable(sparCapFile, 'TextType', 'string');
materials = readtable(materialFile, 'TextType', 'string');
if ~ismember("Ec_psi", string(materials.Properties.VariableNames))
    materials.Ec_psi = materials.E_psi;
end

requiredStringerColumns = [
    "section_key"
    "sec_no"
    "section_type"
    "A_in"
    "B_in"
    "T_in"
    "R_in"
    "area_in2"
    "Ixx_in4"
    "Iyy_in4"
    "rho_xx_in"
    "rho_yy_in"
    "source_table"
    "confidence"
    "notes"
    ];

requiredSparCapColumns = [
    "sec_no"
    "section_type"
    "A_in"
    "B_in"
    "T_in"
    "R_in"
    "area_in2"
    "zbar_in"
    "Iyy_in4"
    "Izz_in4"
    "rho_xx_in"
    "rho_yy_in"
    "source_table"
    "confidence"
    "notes"
    ];

requiredMaterialColumns = [
    "alloy"
    "temper"
    "product_form"
    "E_psi"
    "Ec_psi"
    "G_psi"
    "nu"
    "density_lb_in3"
    "Ftu_psi"
    "Fty_psi"
    "Fcy_psi"
    "Fsu_psi"
    "source_file"
    "page_or_image"
    "table_title"
    "confidence"
    "notes"
    ];

validateColumns(stringers, requiredStringerColumns, 'stringer_options.csv');
validateColumns(sparCaps, requiredSparCapColumns, 'spar_cap_options.csv');
validateColumns(materials, requiredMaterialColumns, 'material_properties.csv');

warnings = strings(0, 1);
if height(materials) == 0
    warnings(end + 1, 1) = ...
        "material_properties.csv has zero rows; material properties must be supplied before analysis.";
else
    materialConfidence = upper(strtrim(string(materials.confidence)));
    activeOrUnclassifiedRows = materialConfidence ~= "HIGH" & ...
        materialConfidence ~= "UNUSED_REFERENCE";
    if any(activeOrUnclassifiedRows)
    warnings(end + 1, 1) = ...
        "One or more active or unclassified material rows are not marked HIGH confidence; review material allowables before using those rows.";
    end
end

designOptions = struct();
designOptions.stringers = stringers;
designOptions.sparCaps = sparCaps;
designOptions.materials = materials;
designOptions.files.stringers = stringerFile;
designOptions.files.sparCaps = sparCapFile;
designOptions.files.materials = materialFile;
designOptions.counts.stringers = height(stringers);
designOptions.counts.sparCaps = height(sparCaps);
designOptions.counts.materials = height(materials);
designOptions.status = "loaded_design_options_only";
designOptions.warnings = warnings;
designOptions.exclusions = [
    "geometry selection"
    "yield checks"
    "buckling"
    "divergence"
    "weight"
    "optimisation"
    ];

fprintf('Loaded design option tables.\n');
fprintf('Stringer options: %d\n', designOptions.counts.stringers);
fprintf('Spar cap options: %d\n', designOptions.counts.sparCaps);
fprintf('Material options: %d\n', designOptions.counts.materials);
for iWarning = 1:numel(warnings)
    warning('loadDesignOptions:DesignOptionWarning', '%s', warnings(iWarning));
end
end

function validateColumns(tableData, requiredColumns, sourceName)
columnNames = string(tableData.Properties.VariableNames);
missingColumns = setdiff(requiredColumns, columnNames);
if ~isempty(missingColumns)
    error('loadDesignOptions:MissingColumns', ...
        '%s is missing required columns: %s', ...
        sourceName, strjoin(missingColumns, ', '));
end
end

