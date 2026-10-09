function section = idealizePanelsToBooms(section, analysisOptions)
% idealizePanelsToBooms Optionally add skin/web panel area to boom areas.
%
% Panel area is assigned to adjacent boom endpoints before centroid,
% inertia, bending, and unit shear-flow coefficients are computed. Physical
% weight is not changed here.

method = string(analysisOptions.panelIdealisationMethod);
enabled = logical(analysisOptions.usePanelIdealisation);

section.idealization = makeEmptyIdealization(enabled, method);
initializeBoomAreaFields();

if ~enabled
    section.idealization.warnings = ...
        "Panel idealisation disabled; section uses original boom-only areas.";
    return;
end

supportedMethods = [
    "split_panel_area_equally_to_end_booms"
    "modified_megson_z_ratio"
    ];
if ~ismember(method, supportedMethods)
    error('idealizePanelsToBooms:UnsupportedMethod', ...
        'Unsupported panel idealisation method: %s', method);
end

boomOnlyCentroid = computeBoomOnlyCentroid();
section.idealization.centroidReference = "boom_only_zc";
section.idealization.boomOnlyCentroid_x_in = boomOnlyCentroid.x_in;
section.idealization.boomOnlyCentroid_z_in = boomOnlyCentroid.z_in;

physicalPanels = collectPhysicalPanelsForIdealization(section);
assignments = repmat(makeAssignment(), 0, 1);
warnings = strings(0, 1);
totalPanelArea = 0.0;
totalAssignedArea = 0.0;

for iPanel = 1:numel(physicalPanels)
    panel = physicalPanels(iPanel);
    panelArea = panel.length_in * panel.thickness_in;
    if ~(isfinite(panelArea) && panelArea > 0)
        error('idealizePanelsToBooms:InvalidPanelArea', ...
            'Panel %s has invalid idealised panel area.', panel.name);
    end

    totalPanelArea = totalPanelArea + panelArea;

    startName = string(panel.startBoomName);
    endName = string(panel.endBoomName);
    startExists = boomExists(startName);
    endExists = boomExists(endName);

    if startExists && endExists
        startBoom = getBoomByName(startName);
        endBoom = getBoomByName(endName);
        assignmentRule = "half_panel_area";
        areaStart = 0.5 * panelArea;
        areaEnd = 0.5 * panelArea;
        fallbackFlag = false;
        warningText = "";
        sigma2_over_sigma1 = NaN;
        sigma1_over_sigma2 = NaN;

        if method == "modified_megson_z_ratio"
            [areaStart, areaEnd, assignmentRule, fallbackFlag, warningText, ...
                sigma2_over_sigma1, sigma1_over_sigma2] = ...
                megsonZRatioAreas(panelArea, startBoom.z_in, endBoom.z_in, ...
                boomOnlyCentroid.z_in, panel.name);
            if strlength(warningText) > 0
                warnings(end + 1, 1) = warningText; %#ok<AGROW>
            end
        end

        [assignments, assigned] = assignToBoom(assignments, panel, startName, ...
            "start", areaStart, panelArea, assignmentRule, startBoom.z_in, ...
            endBoom.z_in, boomOnlyCentroid.z_in, areaStart, areaEnd, ...
            sigma2_over_sigma1, sigma1_over_sigma2, fallbackFlag, warningText);
        totalAssignedArea = totalAssignedArea + assigned;
        [assignments, assigned] = assignToBoom(assignments, panel, endName, ...
            "end", areaEnd, panelArea, assignmentRule, startBoom.z_in, ...
            endBoom.z_in, boomOnlyCentroid.z_in, areaStart, areaEnd, ...
            sigma2_over_sigma1, sigma1_over_sigma2, fallbackFlag, warningText);
        totalAssignedArea = totalAssignedArea + assigned;
    elseif startExists && ~endExists
        startBoom = getBoomByName(startName);
        [assignments, assigned] = assignToBoom(assignments, panel, startName, ...
            "start", panelArea, panelArea, "full_panel_area_reassigned_from_virtual_endpoint", ...
            startBoom.z_in, NaN, boomOnlyCentroid.z_in, panelArea, ...
            0.0, NaN, NaN, true, "virtual end endpoint");
        totalAssignedArea = totalAssignedArea + assigned;
        warnings(end + 1, 1) = sprintf( ...
            "Panel %s has virtual end endpoint %s; full panel area assigned to real start boom %s.", ...
            panel.name, endName, startName);
    elseif ~startExists && endExists
        endBoom = getBoomByName(endName);
        [assignments, assigned] = assignToBoom(assignments, panel, endName, ...
            "end", panelArea, panelArea, "full_panel_area_reassigned_from_virtual_endpoint", ...
            NaN, endBoom.z_in, boomOnlyCentroid.z_in, 0.0, ...
            panelArea, NaN, NaN, true, "virtual start endpoint");
        totalAssignedArea = totalAssignedArea + assigned;
        warnings(end + 1, 1) = sprintf( ...
            "Panel %s has virtual start endpoint %s; full panel area assigned to real end boom %s.", ...
            panel.name, startName, endName);
    else
        error('idealizePanelsToBooms:PanelHasNoRealBoomEndpoints', ...
            'Panel %s has no real boom endpoint for panel idealisation.', panel.name);
    end
end

section.idealization.totalPhysicalPanelArea_in2 = totalPanelArea;
section.idealization.totalPhysicalPanelAreaAssigned_in2 = totalPanelArea;
section.idealization.totalPanelAreaAssigned_in2 = totalAssignedArea;
section.idealization.totalBoomAreaIncrease_in2 = computeTotalBoomIncrease();
section.idealization.assignmentTable = struct2table(assignments);
section.idealization.rows = section.idealization.assignmentTable;
section.idealization.nPanelsIdealized = numel(physicalPanels);
assignmentTableForCounts = section.idealization.assignmentTable;
if isempty(assignmentTableForCounts)
    section.idealization.nZRatioPanels = 0;
    section.idealization.nFallbackPanels = 0;
else
    zRatioRows = assignmentTableForCounts.assignment_rule == "modified_megson_z_ratio" & ...
        ~assignmentTableForCounts.fallback_flag;
    fallbackRows = assignmentTableForCounts.fallback_flag;
    section.idealization.nZRatioPanels = numel(unique( ...
        assignmentTableForCounts.panel_name(zRatioRows)));
    section.idealization.nFallbackPanels = numel(unique( ...
        assignmentTableForCounts.panel_name(fallbackRows)));
end
section.idealization.exclusions = [
    "frontSpar_noseCell loop duplicate is omitted so the physical front spar web is counted once"
    "no physical weight is changed by panel idealisation"
    ];
section.idealization.warnings = warnings;
section.panelIdealisation = section.idealization;

if abs(section.idealization.totalBoomAreaIncrease_in2 - totalAssignedArea) > 1e-8
    error('idealizePanelsToBooms:AreaAssignmentMismatch', ...
        'Total boom area increase does not match assigned panel area.');
end

    function initializeBoomAreaFields()
        section.sparCaps = initializeBoomArray(section.sparCaps);
        section.stringers.top = initializeBoomArray(section.stringers.top);
        section.stringers.bottom = initializeBoomArray(section.stringers.bottom);
        section.stringers.noseTop = initializeBoomArray(section.stringers.noseTop);
        section.stringers.noseBottom = initializeBoomArray(section.stringers.noseBottom);
    end

    function tf = boomExists(boomName)
        [tf, ~, ~] = findBoomLocation(boomName);
    end

    function [assignmentsOut, assignedArea] = assignToBoom(assignmentsIn, panel, boomName, endpoint, areaAdd, panelArea, rule, ...
            z1_in, z2_in, zc_in, B1_in2, B2_in2, sigma2_over_sigma1, ...
            sigma1_over_sigma2, fallbackFlag, warningText)
        if ~(isfinite(areaAdd) && areaAdd > 0)
            error('idealizePanelsToBooms:InvalidAssignedArea', ...
                'Invalid assigned area for panel %s endpoint %s.', panel.name, endpoint);
        end
        [found, groupName, idx] = findBoomLocation(boomName);
        if ~found
            error('idealizePanelsToBooms:BoomNotFound', ...
                'Boom %s was not found for panel idealisation.', boomName);
        end
        boom = getBoom(groupName, idx);
        boom.idealized_area_contribution_in2 = boom.idealized_area_contribution_in2 + areaAdd;
        boom.effective_area_in2 = boom.original_area_in2 + boom.idealized_area_contribution_in2;
        boom.area_in2 = boom.effective_area_in2;
        boom.sourcePanels(end + 1, 1) = panel.name;
        setBoom(groupName, idx, boom);

        assignment = makeAssignment();
        assignment.panel_name = panel.name;
        assignment.panel_cell = panel.cell;
        assignment.startBoomName = string(panel.startBoomName);
        assignment.endBoomName = string(panel.endBoomName);
        assignment.assigned_boom_name = boomName;
        assignment.assigned_boom_group = string(groupName);
        assignment.endpoint = string(endpoint);
        assignment.panel_length_in = panel.length_in;
        assignment.panel_thickness_in = panel.thickness_in;
        assignment.panel_area_in2 = panelArea;
        assignment.b_in = panel.length_in;
        assignment.t_in = panel.thickness_in;
        assignment.z1_in = z1_in;
        assignment.z2_in = z2_in;
        assignment.zc_in = zc_in;
        assignment.B1_in2 = B1_in2;
        assignment.B2_in2 = B2_in2;
        assignment.sigma2_over_sigma1 = sigma2_over_sigma1;
        assignment.sigma1_over_sigma2 = sigma1_over_sigma2;
        assignment.assigned_area_in2 = areaAdd;
        assignment.assignment_rule = string(rule);
        assignment.method_used = string(rule);
        assignment.fallback_flag = logical(fallbackFlag);
        assignment.fallback_warning = string(warningText);
        assignmentsOut = [assignmentsIn; assignment]; %#ok<AGROW>
        assignedArea = areaAdd;
    end

    function [found, groupName, idx] = findBoomLocation(boomName)
        groups = ["sparCaps", "top", "bottom", "noseTop", "noseBottom"];
        found = false;
        groupName = "";
        idx = NaN;
        for iGroup = 1:numel(groups)
            booms = getBoomArray(groups(iGroup));
            if isempty(booms)
                continue;
            end
            matches = [booms.name] == string(boomName);
            if nnz(matches) == 1
                found = true;
                groupName = groups(iGroup);
                idx = find(matches, 1);
                return;
            elseif nnz(matches) > 1
                error('idealizePanelsToBooms:DuplicateBoomName', ...
                    'Boom name %s appears more than once.', boomName);
            end
        end
    end

    function booms = getBoomArray(groupName)
        switch string(groupName)
            case "sparCaps"
                booms = section.sparCaps;
            case "top"
                booms = section.stringers.top;
            case "bottom"
                booms = section.stringers.bottom;
            case "noseTop"
                booms = section.stringers.noseTop;
            case "noseBottom"
                booms = section.stringers.noseBottom;
            otherwise
                error('idealizePanelsToBooms:UnknownBoomGroup', ...
                    'Unknown boom group %s.', groupName);
        end
    end

    function boom = getBoom(groupName, idx)
        booms = getBoomArray(groupName);
        boom = booms(idx);
    end

    function boom = getBoomByName(boomName)
        [found, groupName, idx] = findBoomLocation(boomName);
        if ~found
            error('idealizePanelsToBooms:BoomNotFound', ...
                'Boom %s was not found for panel idealisation.', boomName);
        end
        boom = getBoom(groupName, idx);
    end

    function setBoom(groupName, idx, boom)
        switch string(groupName)
            case "sparCaps"
                section.sparCaps(idx) = boom;
            case "top"
                section.stringers.top(idx) = boom;
            case "bottom"
                section.stringers.bottom(idx) = boom;
            case "noseTop"
                section.stringers.noseTop(idx) = boom;
            case "noseBottom"
                section.stringers.noseBottom(idx) = boom;
            otherwise
                error('idealizePanelsToBooms:UnknownBoomGroup', ...
                    'Unknown boom group %s.', groupName);
        end
    end

    function totalIncrease = computeTotalBoomIncrease()
        allBooms = [
            section.sparCaps(:)
            section.stringers.top(:)
            section.stringers.bottom(:)
            section.stringers.noseTop(:)
            section.stringers.noseBottom(:)
            ];
        totalIncrease = sum([allBooms.idealized_area_contribution_in2]);
    end

    function centroid = computeBoomOnlyCentroid()
        allBooms = [
            section.sparCaps(:)
            section.stringers.top(:)
            section.stringers.bottom(:)
            section.stringers.noseTop(:)
            section.stringers.noseBottom(:)
            ];
        totalOriginalArea = sum([allBooms.original_area_in2]);
        centroid.x_in = sum([allBooms.x_in] .* [allBooms.original_area_in2]) / totalOriginalArea;
        centroid.z_in = sum([allBooms.z_in] .* [allBooms.original_area_in2]) / totalOriginalArea;
    end
end

function [B1, B2, rule, fallbackFlag, warningText, sigma2_over_sigma1, sigma1_over_sigma2] = ...
    megsonZRatioAreas(panelArea, z1, z2, zc, panelName)
tol = 1.0e-8;
rule = "modified_megson_z_ratio";
fallbackFlag = false;
warningText = "";
sigma2_over_sigma1 = NaN;
sigma1_over_sigma2 = NaN;

dz1 = z1 - zc;
dz2 = z2 - zc;
if abs(dz1) < tol || abs(dz2) < tol
    [B1, B2, rule, fallbackFlag, warningText] = equalSplitFallback( ...
        panelArea, panelName, "endpoint near boom-only centroid");
    return;
end

sigma2_over_sigma1 = dz2 / dz1;
sigma1_over_sigma2 = dz1 / dz2;
B1 = panelArea / 6.0 * (2.0 + sigma2_over_sigma1);
B2 = panelArea / 6.0 * (2.0 + sigma1_over_sigma2);

if dz1 * dz2 < 0
    [B1, B2, rule, fallbackFlag, warningText] = equalSplitFallback( ...
        panelArea, panelName, "endpoints straddle boom-only centroid");
    return;
end
if ~isfinite(B1) || ~isfinite(B2) || B1 <= 0 || B2 <= 0
    [B1, B2, rule, fallbackFlag, warningText] = equalSplitFallback( ...
        panelArea, panelName, "non-finite or non-positive z-ratio boom area");
end
end

function [B1, B2, rule, fallbackFlag, warningText] = equalSplitFallback(panelArea, panelName, reason)
B1 = 0.5 * panelArea;
B2 = 0.5 * panelArea;
rule = "fallback_equal_split";
fallbackFlag = true;
warningText = "Panel " + string(panelName) + ...
    " used equal-split fallback for z-ratio panel idealisation: " + string(reason) + ".";
end

function idealization = makeEmptyIdealization(enabled, method)
idealization = struct();
idealization.enabled = enabled;
idealization.method = method;
idealization.centroidReference = "";
idealization.boomOnlyCentroid_x_in = NaN;
idealization.boomOnlyCentroid_z_in = NaN;
idealization.totalPhysicalPanelArea_in2 = 0.0;
idealization.totalPhysicalPanelAreaAssigned_in2 = 0.0;
idealization.totalPanelAreaAssigned_in2 = 0.0;
idealization.totalBoomAreaIncrease_in2 = 0.0;
idealization.assignmentTable = table();
idealization.rows = table();
idealization.nPanelsIdealized = 0;
idealization.nZRatioPanels = 0;
idealization.nFallbackPanels = 0;
idealization.exclusions = strings(0, 1);
idealization.warnings = strings(0, 1);
idealization.units.area = "in^2";
end

function booms = initializeBoomArray(booms)
for iBoom = 1:numel(booms)
    if ~isfield(booms(iBoom), 'original_area_in2') || ...
            ~isfinite(booms(iBoom).original_area_in2)
        booms(iBoom).original_area_in2 = booms(iBoom).area_in2;
    end
    booms(iBoom).idealized_area_contribution_in2 = 0.0;
    booms(iBoom).effective_area_in2 = booms(iBoom).original_area_in2;
    booms(iBoom).area_in2 = booms(iBoom).original_area_in2;
    booms(iBoom).sourcePanels = strings(0, 1);
end
end

function panels = collectPhysicalPanelsForIdealization(section)
panels = [
    section.webs.top(:)
    section.webs.bottom(:)
    section.webs.rearSpar(:)
    section.webs.frontSpar(:)
    section.webs.lowerNose(:)
    section.webs.upperNose(:)
    ];
end

function assignment = makeAssignment()
assignment = struct();
assignment.panel_name = "";
assignment.panel_cell = NaN;
assignment.startBoomName = "";
assignment.endBoomName = "";
assignment.assigned_boom_name = "";
assignment.assigned_boom_group = "";
assignment.endpoint = "";
assignment.panel_length_in = NaN;
assignment.panel_thickness_in = NaN;
assignment.panel_area_in2 = NaN;
assignment.b_in = NaN;
assignment.t_in = NaN;
assignment.z1_in = NaN;
assignment.z2_in = NaN;
assignment.zc_in = NaN;
assignment.B1_in2 = NaN;
assignment.B2_in2 = NaN;
assignment.sigma2_over_sigma1 = NaN;
assignment.sigma1_over_sigma2 = NaN;
assignment.assigned_area_in2 = NaN;
assignment.assignment_rule = "";
assignment.method_used = "";
assignment.fallback_flag = false;
assignment.fallback_warning = "";
end
