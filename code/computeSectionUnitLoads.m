function section = computeSectionUnitLoads(wing, analysisOptions)
% computeSectionUnitLoads Thin-walled section model of the two-cell wing box.
%
%   section = computeSectionUnitLoads(wing, analysisOptions) builds the
%   idealised cross-section (booms + walls on the NACA 2412 contour) and
%   returns every unit-load coefficient used downstream:
%     - centroid and centroidal second moments of the boom areas;
%     - direct-stress coefficient of every boom for unit Mx and Mz;
%     - shear-flow coefficient of every wall for unit Vx, unit Vz (acting
%       through the quarter-chord load reference) and unit My;
%     - shear centre location.
%
% Theory (standard thin-walled beam theory, e.g. T.H.G. Megson, "Aircraft
% Structures for Engineering Students", and E.F. Bruhn, "Analysis and
% Design of Flight Vehicle Structures"):
%   1. Idealised section: skins and spar webs carry shear only; direct
%      stress is carried by booms (spar caps, stringers and, optionally,
%      equivalent panel areas assigned by idealizePanelsToBooms).
%   2. Unsymmetric bending: sigma = f(x, z, Mx, Mz, Ixx, Izz, Ixz).
%   3. Shear of a closed multicell section: each cell is "cut" at its
%      first wall, the open-section flow is accumulated from boom to boom
%      (q_k = q_{k-1} - dq_boom), and one constant flow per cell restores
%      zero twist (sum over the cell of q ds/t = 0). Walls shared by two
%      cells carry the difference of the two cell constants.
%   4. Shear centre: moment equilibrium of the unit shear flows about a
%      reference point (front-spar station on the chord line), using
%      2 x swept area for each wall.
%   5. Pure torsion: Bredt-Batho multicell torsion; all cells share the
%      same rate of twist and sum(2 A_i q_i) = T.
%
% The output field names are the interface consumed by the rest of the
% analysis (stresses, shear flows, buckling, divergence, weight, export).

if nargin < 2 || isempty(analysisOptions)
    analysisOptions = defaultAnalysisOptions();
end
checkInputs(wing, analysisOptions);

c = wing.planform.chord_in;
xcF = wing.box.frontSpar_xc;
xcR = wing.box.rearSpar_xc;
xf = xcF * c;             % front-spar station (in)
xr = xcR * c;             % rear-spar station (in)
xRef = wing.box.frontSpar_x_in;   % moment reference: front spar on the chord line

%% 1. Geometry: booms (nodes) and walls (edges) on the contour
section = struct();
section.units = struct('length', "in", 'area', "in^2", 'inertia', "in^4");
section.analysis.stage = "stage_3C2_geometry_centroid_inertia_unit_bending_unit_shear_coefficients";
section.analysis.source = "Two-cell thin-walled idealised section model (booms and panels on the NACA 2412 contour)";
section.analysis.exclusions = ["load-case scaling"; "failure checks"; "buckling"; ...
    "von Mises"; "divergence"; "weight"; "optimisation"];
section.metadata.panelIdealization = ...
    "Optional panel idealisation into boom areas is implemented when analysisOptions.usePanelIdealisation = true; final-design method is modified Megson z-ratio panel idealisation. When disabled, the section uses boom-only areas. Detailed plate stress recovery remains outside computeSectionUnitLoads.";

capArea = [wing.sparCaps.A1_in2, wing.sparCaps.A2_in2, wing.sparCaps.A3_in2, wing.sparCaps.A4_in2];
section.sparCaps = [
    contourBoom("A1_front_upper", xf, xcF, c, true,  capArea(1))
    contourBoom("A2_front_lower", xf, xcF, c, false, capArea(2))
    contourBoom("A3_rear_upper",  xr, xcR, c, true,  capArea(3))
    contourBoom("A4_rear_lower",  xr, xcR, c, false, capArea(4))
    ];
section.stringers.top        = stringerRow("top",        wing.stringers.numTop,        wing.stringers.areaTop_in2,        xf,  xr,  c, true);
section.stringers.bottom     = stringerRow("bottom",     wing.stringers.numBottom,     wing.stringers.areaBottom_in2,     xr,  xf,  c, false);
section.stringers.noseTop    = stringerRow("noseTop",    wing.stringers.numNoseTop,    wing.stringers.areaNoseTop_in2,    0.0, xf,  c, true);
section.stringers.noseBottom = stringerRow("noseBottom", wing.stringers.numNoseBottom, wing.stringers.areaNoseBottom_in2, xf,  0.0, c, false);

leUpper = contourBoom("leadingEdgeUpper", 0.0, 0.0, c, true,  0.0);
leLower = contourBoom("leadingEdgeLower", 0.0, 0.0, c, false, 0.0);
cap = section.sparCaps;
t = wing.thickness;
section.webs.top            = wallChain("top",       1, [cap(1); section.stringers.top(:);        cap(3)], t.t3_upperMainSkin_in);
section.webs.bottom         = wallChain("bottom",    1, [cap(4); section.stringers.bottom(:);     cap(2)], t.t4_lowerMainSkin_in);
section.webs.rearSpar       = wall("rearSpar",           1, cap(3), cap(4), t.t2_rearSpar_in);
section.webs.frontSpar      = wall("frontSpar_mainCell", 1, cap(2), cap(1), t.t1_frontSpar_in);
section.webs.lowerNose      = wallChain("lowerNose", 2, [cap(2); section.stringers.noseBottom(:); leLower], t.lowerNoseSkin_in);
section.webs.upperNose      = wallChain("upperNose", 2, [leUpper; section.stringers.noseTop(:);   cap(1)], t.upperNoseSkin_in);
section.webs.frontSparCell2 = wall("frontSpar_noseCell", 2, cap(1), cap(2), t.t1_frontSpar_in);

% Optional equivalent-boom areas for the skin and web panels.
section = idealizePanelsToBooms(section, analysisOptions);

%% 2. Section properties of the boom set
B = boomTable(section);
A = B.area;
section.centroid.x_in = sum(B.x .* A) / sum(A);
section.centroid.z_in = sum(B.z .* A) / sum(A);
xb = B.x - section.centroid.x_in;
zb = B.z - section.centroid.z_in;
section.inertia.Ix_in4 = sum(A .* zb.^2);
section.inertia.Iz_in4 = sum(A .* xb.^2);
section.inertia.Ixz_in4 = sum(A .* xb .* zb);
if section.idealization.enabled
    section.centroid.method = "Centroid uses effective boom areas with optional panel idealisation.";
    section.inertia.method = "Inertias use effective boom areas with optional panel idealisation.";
else
    section.centroid.method = "Boom-only centroid (area-weighted boom positions); panel idealisation disabled.";
    section.inertia.method = "Boom-only centroidal inertias (parallel-axis sums over booms); panel idealisation disabled.";
end
I = [section.inertia.Ix_in4, section.inertia.Iz_in4, section.inertia.Ixz_in4];
detI = I(1) * I(2) - I(3)^2;
if ~(isfinite(detI) && detI > 0)
    error('computeSectionUnitLoads:InvalidBendingDenominator', ...
        'Ix*Iz - Ixz^2 must be positive and finite before computing bending coefficients.');
end

%% 3. Direct-stress coefficients for unit bending moments
section.unitBending = bendingCoefficients(B, xb, zb, I, detI);

%% 4. Shear-flow solution of the two-cell section
% Each cell is described as a closed loop of walls traversed in one
% direction. The front spar belongs to both loops (with opposite
% traversal), which is how the shared wall enters the compatibility
% equations.
loops = {
    [loopWalls(section.webs.top, "top", true, c, xRef, B); ...
     loopWalls(section.webs.rearSpar, "rearSpar", NaN, c, xRef, B); ...
     loopWalls(section.webs.bottom, "bottom", false, c, xRef, B); ...
     loopWalls(section.webs.frontSpar, "frontSpar", NaN, c, xRef, B)]
    [loopWalls(section.webs.lowerNose, "lowerNose", false, c, xRef, B); ...
     loopWalls(section.webs.upperNose, "upperNose", true, c, xRef, B); ...
     loopWalls(section.webs.frontSparCell2, "frontSparCell2", NaN, c, xRef, B)]
    };
sharedGroup = ["frontSpar"; "frontSparCell2"];   % shared wall as seen from each loop
cellNames = ["mainBox"; "noseCell"];

for k = 1:numel(loops)
    loops{k} = openSectionFlows(loops{k}, section.centroid, I);
end
section.unitShear.stage3A = openFlowSummary(loops, section);

[qs, compat] = cellConstants(loops, sharedGroup, cellNames);
section.unitShear.stage3B.compatibility = compat;

[preTorsion, sharedWall, physical, cellArea] = closedSectionFlows(loops, qs, sharedGroup, section.unitShear.stage3A.panelCount);
[section.shearCenter, preTorsion.momentSums] = shearCentre(loops, qs, cellArea, sharedGroup, xRef);
section.unitShear.stage3B.preTorsion = preTorsion;
section.unitShear.stage3B.sharedWallPhysical = sharedWall;

section.unitShear.stage3C.torsion = unitTorsion(physical, compat, cellArea, preTorsion);
section.unitShear.stage3C.unitTotals = unitTotals(section, c, preTorsion);

checkSection(section, wing);
end

%% ---------------------------------------------------------------- geometry
function boom = contourBoom(name, x_in, xc, c, isUpper, area_in2)
% Boom on the airfoil contour at chordwise station x_in (= xc * chord).
boom = struct();
boom.name = string(name);
boom.x_in = x_in;
boom.z_in = nacaSurfaceOrdinate(xc, double(isUpper)) * c;
boom.area_in2 = area_in2;
boom.original_area_in2 = area_in2;
boom.idealized_area_contribution_in2 = 0.0;
boom.effective_area_in2 = area_in2;
boom.sourcePanels = strings(0, 1);
end

function row = stringerRow(prefix, n, area_in2, xFrom, xTo, c, isUpper)
% n equally spaced stringers strictly between xFrom and xTo.
row = repmat(contourBoom("", 0.0, 0.0, c, isUpper, area_in2), n, 1);
step = (xTo - xFrom) / (n + 1);
for k = 1:n
    x = xFrom + step * k;
    row(k) = contourBoom(sprintf('%s_%02d', prefix, k), x, x / c, c, isUpper, area_in2);
end
end

function walls = wallChain(prefix, cellId, nodes, thickness_in)
nWalls = numel(nodes) - 1;
walls = repmat(wall("", cellId, nodes(1), nodes(end), thickness_in), nWalls, 1);
for k = 1:nWalls
    walls(k) = wall(sprintf('%s_%02d', prefix, k), cellId, nodes(k), nodes(k + 1), thickness_in);
end
end

function w = wall(name, cellId, fromNode, toNode, thickness_in)
w = struct();
w.name = string(name);
w.cell = cellId;
w.startBoomName = string(fromNode.name);
w.endBoomName = string(toNode.name);
w.xStart_in = fromNode.x_in;
w.xEnd_in = toNode.x_in;
w.zStart_in = fromNode.z_in;
w.zEnd_in = toNode.z_in;
w.thickness_in = thickness_in;
w.dx_in = w.xEnd_in - w.xStart_in;
w.dz_in = w.zEnd_in - w.zStart_in;
w.length_in = hypot(w.dx_in, w.dz_in);
w.areaSource = "geometry_only_stage_1";
end

function B = boomTable(section)
% Flat list of the area-carrying booms in a fixed order.
all = [section.sparCaps(:); section.stringers.top(:); section.stringers.bottom(:); ...
    section.stringers.noseTop(:); section.stringers.noseBottom(:)];
B.name = [all.name]';
B.x = [all.x_in]';
B.z = [all.z_in]';
B.area = [all.area_in2]';
B.group = [repmat("sparCaps", numel(section.sparCaps), 1)
    repmat("topStringers", numel(section.stringers.top), 1)
    repmat("bottomStringers", numel(section.stringers.bottom), 1)
    repmat("noseTopStringers", numel(section.stringers.noseTop), 1)
    repmat("noseBottomStringers", numel(section.stringers.noseBottom), 1)];
end

%% ----------------------------------------------------------------- bending
function ub = bendingCoefficients(B, xb, zb, I, detI)
n = numel(B.x);
sMx = zeros(n, 1);
sMz = zeros(n, 1);
for k = 1:n
    sMx(k) = unsymmetricBendingStress(xb(k), zb(k), 1.0, 0.0, I(1), I(2), I(3));
    sMz(k) = unsymmetricBendingStress(xb(k), zb(k), 0.0, 1.0, I(1), I(2), I(3));
end
entries = struct('name', cellstr(B.name), 'group', cellstr(B.group), ...
    'x_in', num2cell(B.x), 'z_in', num2cell(B.z), 'dx_in', num2cell(xb), ...
    'dz_in', num2cell(zb), 'sigY_Mx_per_lbin', num2cell(sMx), ...
    'sigY_Mz_per_lbin', num2cell(sMz));
for k = 1:n
    entries(k).name = string(entries(k).name);
    entries(k).group = string(entries(k).group);
end
[~, iMx] = max(abs(sMx));
[~, iMz] = max(abs(sMz));
ub = struct();
ub.Mx_lbin = 1.0;
ub.Mz_lbin = 1.0;
ub.entries = entries;
ub.denominator_in8 = detI;
ub.maxAbsSigY_Mx_per_lbin = abs(sMx(iMx));
ub.maxAbsSigY_Mx_location = entries(iMx);
ub.maxAbsSigY_Mz_per_lbin = abs(sMz(iMz));
ub.maxAbsSigY_Mz_location = entries(iMz);
ub.method = "Unit bending coefficients from the generalised flexure formula (unsymmetricBendingStress).";
ub.units.coefficient = "psi per lb-in";
end

%% ------------------------------------------------------ wall properties
function walls = loopWalls(walls, group, isUpper, c, xf, B)
% Adds the properties a wall needs inside a cell loop: developed length,
% ds/t, area of the boom at the start of the wall, and swept area about
% the reference point (xf, 0).
out = cell(numel(walls), 1);
for k = 1:numel(walls)
    w = walls(k);
    w.group = string(group);
    w.dx_in = w.xEnd_in - w.xStart_in;
    w.dz_in = w.zEnd_in - w.zStart_in;
    if isnan(isUpper)
        w.surfaceFlag = NaN;
        w.ds_in = abs(w.dz_in);
        contourArea = 0.0;
    else
        w.surfaceFlag = double(isUpper);
        w.ds_in = surfaceArcLength(w.xStart_in / c, w.xEnd_in / c, double(isUpper)) * c;
        contourArea = signedContourArea(w.xStart_in / c, w.xEnd_in / c, isUpper) * c^2;
    end
    w.ds_over_t = w.ds_in / w.thickness_in;
    w.dp_area_in2 = startBoomArea(w.startBoomName, B);
    % Swept area of the wall seen from (xf, 0): area under the contour
    % minus the triangles to the reference point at both wall ends.
    w.sweptArea_in2 = contourArea - 0.5 * ((w.xEnd_in - xf) * w.zEnd_in - (w.xStart_in - xf) * w.zStart_in);
    if group == "frontSpar" || group == "frontSparCell2"
        w.sweptArea_in2 = 0.0;   % wall lies on the reference line
    end
    w.sweptAreaSource = "Swept area between the chord line and the contour (surfaceAreaIntegral) plus the triangular terms to the moment reference.";
    out{k} = w;
end
walls = vertcat(out{:});
end

function a = signedContourArea(xa, xb, isUpper)
% Integral of z/c from x/c = xa to xb along the chosen surface.
a = surfaceAreaIntegral(min(xa, xb), max(xa, xb), double(isUpper));
if ~isUpper
    a = -a;      % lower surface lies below the chord line
end
if xb < xa
    a = -a;
end
end

function area = startBoomArea(name, B)
idx = find(B.name == name, 1);
if isempty(idx)
    area = 0.0;  % leading-edge contour points carry no boom area
else
    area = B.area(idx);
end
end

%% -------------------------------------------------- open-section flows
function walls = openSectionFlows(walls, centroid, I)
% Cut the cell at the start of its first wall (q = 0) and accumulate the
% flow change produced by each boom met along the loop.
n = numel(walls);
dqX = zeros(n, 1);
dqZ = zeros(n, 1);
for k = 2:n
    x = walls(k).xStart_in - centroid.x_in;
    z = walls(k).zStart_in - centroid.z_in;
    dqX(k) = boomShearFlowIncrement(x, z, 1.0, 0.0, I(1), I(2), I(3), walls(k).dp_area_in2);
    dqZ(k) = boomShearFlowIncrement(x, z, 0.0, 1.0, I(1), I(2), I(3), walls(k).dp_area_in2);
end
qX = -cumsum(dqX);
qZ = -cumsum(dqZ);
out = cell(n, 1);
for k = 1:n
    w = walls(k);
    w.dP_X_per_lb = dqX(k);
    w.dP_Z_per_lb = dqZ(k);
    w.qPrime_X_per_lb = qX(k);
    w.qPrime_Z_per_lb = qZ(k);
    w.q_dS_over_t_X = qX(k) * w.ds_over_t;
    w.q_dS_over_t_Z = qZ(k) * w.ds_over_t;
    w.two_A_qprime_X = 2 * w.sweptArea_in2 * qX(k);
    w.two_A_qprime_Z = 2 * w.sweptArea_in2 * qZ(k);
    w.qp_dx_X = qX(k) * w.dx_in;
    w.qp_dx_Z = qZ(k) * w.dx_in;
    w.qp_dz_X = qX(k) * w.dz_in;
    w.qp_dz_Z = qZ(k) * w.dz_in;
    out{k} = w;
end
walls = vertcat(out{:});
end

function s = openFlowSummary(loops, section)
all = vertcat(loops{:});
[~, iX] = max(abs([all.qPrime_X_per_lb]));
[~, iZ] = max(abs([all.qPrime_Z_per_lb]));
s = struct();
s.Vx_unit_lb = 1.0;
s.Vz_unit_lb = 1.0;
s.panels = all;
s.mainBoxPanels = loops{1};
s.noseCellPanels = loops{2};
s.panelCount = numel(all);
webs = struct2cell(section.webs);
s.expectedPanelCount = sum(cellfun(@numel, webs));
s.mainBoxLoop.first = loops{1}(1);
s.mainBoxLoop.last = loops{1}(end);
s.noseCellLoop.first = loops{2}(1);
s.noseCellLoop.last = loops{2}(end);
s.maxAbsQPrime_X_per_lb = abs(all(iX).qPrime_X_per_lb);
s.maxAbsQPrime_X_location = all(iX);
s.maxAbsQPrime_Z_per_lb = abs(all(iZ).qPrime_Z_per_lb);
s.maxAbsQPrime_Z_location = all(iZ);
s.method = "Stage 3A open-section qPrime values from boom shear-flow increments (boomShearFlowIncrement).";
s.units.qPrime = "lb/in per lb";
s.units.ds_over_t = "dimensionless";
s.units.sweptArea = "in^2";
s.exclusions = ["corrected physical shear flows"; "shear center"; "torsional shear flow"; ...
    "qtot_X, qtot_Z, qtot_Y"; "load-case shear-flow scaling"; "shear stress"; "von Mises"; ...
    "buckling"; "yield"; "divergence"; "weight"; "optimisation"];
end

%% -------------------------------------------- multicell compatibility
function [qs, compat] = cellConstants(loops, sharedGroup, cellNames)
% Zero-twist condition for every cell: sum over the cell walls of
% (q_open + q_cell - q_neighbour) ds/t = 0, solved for the cell constants.
nCells = numel(loops);
K = zeros(nCells);
rX = zeros(nCells, 1);
rZ = zeros(nCells, 1);
sharedWall = cell(nCells, 1);
for i = 1:nCells
    w = loops{i};
    K(i, i) = sum([w.ds_over_t]);
    sharedWall{i} = w([w.group] == sharedGroup(i));
    rX(i) = -sum([w.q_dS_over_t_X]);
    rZ(i) = -sum([w.q_dS_over_t_Z]);
end
K(1, 2) = -sharedWall{1}.ds_over_t;
K(2, 1) = -sharedWall{2}.ds_over_t;
qs.X = K \ rX;
qs.Z = K \ rZ;

compat = struct();
compat.A11 = K(1, 1);
compat.A22 = K(2, 2);
compat.A12 = K(1, 2);
compat.A21 = K(2, 1);
compat.Amat = K;
compat.Bmat_X = rX;
compat.Bmat_Z = rZ;
compat.qs_X_per_lb = qs.X;
compat.qs_Z_per_lb = qs.Z;
compat.conditionNumber = cond(K);
compat.cellLoopNames = cellNames;
compat.sharedWall.mainBoxPanelName = sharedWall{1}.name;
compat.sharedWall.noseCellPanelName = sharedWall{2}.name;
compat.sharedWall.mainBox_ds_over_t = sharedWall{1}.ds_over_t;
compat.sharedWall.noseCell_ds_over_t = sharedWall{2}.ds_over_t;
compat.sharedWall.absDifference = abs(sharedWall{1}.ds_over_t - sharedWall{2}.ds_over_t);
compat.sharedWall.relativeDifference = compat.sharedWall.absDifference / ...
    max(abs([sharedWall{1}.ds_over_t, sharedWall{2}.ds_over_t]));
compat.method = "Stage 3B-1 compatibility from the two-cell zero-twist equations.";
compat.units.Amat = "dimensionless";
compat.units.Bmat = "lb/in per lb times ds/t";
compat.units.qs = "lb/in per lb";
compat.exclusions = ["torsional shear flow"; "total shear flow with torsion"; ...
    "qtot_X, qtot_Z, qtot_Y"; "load-case shear-flow scaling"; "shear stress"; "von Mises"; ...
    "buckling"; "yield"; "divergence"; "weight"; "optimisation"];
end

%% ------------------------------------------------ closed-section flows
function [pre, sharedOut, physical, cellArea] = closedSectionFlows(loops, qs, sharedGroup, panelCount)
nCells = numel(loops);
closed = loops;
cellArea = zeros(nCells, 1);
for i = 1:nCells
    for k = 1:numel(closed{i})
        closed{i}(k).qcorr_X_per_lb = closed{i}(k).qPrime_X_per_lb + qs.X(i);
        closed{i}(k).qcorr_Z_per_lb = closed{i}(k).qPrime_Z_per_lb + qs.Z(i);
        closed{i}(k).qcorrSource = "pre_torsion_qPrime_plus_cell_qs";
    end
    own = closed{i}([closed{i}.group] ~= sharedGroup(i));
    cellArea(i) = sum([own.sweptArea_in2]);
end

% The shared front spar exists once physically: main-cell open flow plus
% the main-cell constant minus the nose-cell constant.
mainCopy = closed{1}([closed{1}.group] == sharedGroup(1));
noseCopy = closed{2}([closed{2}.group] == sharedGroup(2));
sharedOut = mainCopy;
sharedOut.name = "frontSpar_physical_shared";
sharedOut.group = "frontSparSharedPhysical";
sharedOut.sourceMainCellPanelName = mainCopy.name;
sharedOut.sourceNoseCellPanelName = noseCopy.name;
sharedOut.qPrime_X_main_per_lb = mainCopy.qPrime_X_per_lb;
sharedOut.qPrime_Z_main_per_lb = mainCopy.qPrime_Z_per_lb;
sharedOut.qPrime_X_nose_per_lb = noseCopy.qPrime_X_per_lb;
sharedOut.qPrime_Z_nose_per_lb = noseCopy.qPrime_Z_per_lb;
sharedOut.qcorr_X_per_lb = mainCopy.qPrime_X_per_lb + qs.X(1) - qs.X(2);
sharedOut.qcorr_Z_per_lb = mainCopy.qPrime_Z_per_lb + qs.Z(1) - qs.Z(2);
sharedOut.signConventionNote = "Physical shared-wall flow: qPrime_main + qs_cell1 - qs_cell2; no torsional qt included.";

sharedPanel = mainCopy;
sharedPanel.name = sharedOut.name;
sharedPanel.group = sharedOut.group;
sharedPanel.qcorr_X_per_lb = sharedOut.qcorr_X_per_lb;
sharedPanel.qcorr_Z_per_lb = sharedOut.qcorr_Z_per_lb;
sharedPanel.qcorrSource = "physical_shared_wall_qPrime_main_plus_qs_cell1_minus_qs_cell2";

physical = [closed{1}([closed{1}.group] ~= sharedGroup(1))
    closed{2}([closed{2}.group] ~= sharedGroup(2))
    sharedPanel];

pre = struct();
pre.Vx_unit_lb = 1.0;
pre.Vz_unit_lb = 1.0;
pre.mainBoxPanels = closed{1};
pre.noseCellPanels = closed{2};
pre.physicalPanels = physical;
pre.physicalPanelCount = numel(physical);
pre.expectedPhysicalPanelCount = panelCount - 1;
pre.area.Area1_mainCell_in2 = cellArea(1);
pre.area.Area2_noseCell_in2 = cellArea(2);
pre.area.source = "Area1/Area2 are the enclosed cell areas; shared front spar has zero swept area and is not double-counted.";
pre.momentSums = struct();
pre.method = "Stage 3B-2 corrected pre-torsion shear flows and shear center from 2Aq moment equilibrium.";
pre.units.qcorr = "lb/in per lb";
pre.units.area = "in^2";
pre.units.momentSum = "lb-in per lb";
pre.exclusions = ["My = 1 torsional shear flow"; "qt_Y"; "qtot_X, qtot_Z, qtot_Y"; ...
    "load-case shear-flow scaling"; "shear stress"; "von Mises"; "buckling"; "yield"; ...
    "divergence"; "weight"; "optimisation"];
end

function [sc, sums] = shearCentre(loops, qs, cellArea, sharedGroup, xf)
% Moment of the unit shear flows about (xf, 0); the applied unit shear
% force must produce the same moment when it acts at the shear centre.
mX = 0;
mZ = 0;
for i = 1:numel(loops)
    own = loops{i}([loops{i}.group] ~= sharedGroup(i));
    mX = mX + sum([own.two_A_qprime_X]);
    mZ = mZ + sum([own.two_A_qprime_Z]);
end
mX = mX + 2 * sum(qs.X .* cellArea);
mZ = mZ + 2 * sum(qs.Z .* cellArea);

sums = struct();
sums.sum_2_a_q_X = mX;
sums.sum_2_a_q_Z = mZ;
sums.reference = "2Aq moment sums are taken about the front spar reference line.";

sc = struct();
sc.x_sc_in = mZ / 1.0 + xf;
sc.z_sc_in = -mX / 1.0;
sc.xReference = "x_sc_in measured from leading edge.";
sc.zReference = "z_sc_in measured from airfoil chord line.";
sc.momentReference = "The x moment arm is computed about the front spar line, then frontSpar_x_in is added.";
sc.signConvention = "Sign convention: x_sc = sum_2_a_q_Z/Vz + frontSpar_x; z_sc = -sum_2_a_q_X/Vx.";
sc.frontSpar_x_in = xf;
sc.source = "Stage 3B-2 shear center from pre-torsion unit Vx/Vz flows; no torsional qt terms included.";
end

%% --------------------------------------------------------------- torsion
function tor = unitTorsion(physical, compat, cellArea, pre)
% Bredt-Batho two-cell torsion for My = 1 lb-in. Unknowns [q1; q2; Gtheta]:
%   (1/(2 A_i)) * (sum_i ds/t * q_i - shared ds/t * q_j) - G*theta = 0
%   2 A_1 q_1 + 2 A_2 q_2 = My
T = 1.0;
s = compat.sharedWall.mainBox_ds_over_t;
M = [compat.A11 / cellArea(1), -s / cellArea(1), -1
     -s / cellArea(2), compat.A22 / cellArea(2), -1
     2 * cellArea(1), 2 * cellArea(2), 0];
sol = M \ [0; 0; T];
qt = sol(1:2);
qtShared = qt(1) - qt(2);
recovered = 2 * cellArea(1) * qt(1) + 2 * cellArea(2) * qt(2);

physicalOut = physical;
for k = 1:numel(physicalOut)
    switch physicalOut(k).group
        case "frontSparSharedPhysical"
            physicalOut(k).qt_Y_per_lbin = qtShared;
            physicalOut(k).qtSource = "physical_shared_wall_qt_cell1_minus_qt_cell2";
        otherwise
            physicalOut(k).qt_Y_per_lbin = qt(physicalOut(k).cell);
            physicalOut(k).qtSource = "pure_My_unit_torsion_only";
    end
end

tor = struct();
tor.My_unit_lbin = T;
tor.qt_Y_per_lbin = qt;
tor.qt_Y_mainCell_per_lbin = qt(1);
tor.qt_Y_noseCell_per_lbin = qt(2);
tor.qt_Y_shared_per_lbin = qtShared;
tor.q1t_over_q2t = qt(1) / qt(2);
tor.recoveredTorque_lbin = recovered;
tor.torqueError_lbin = recovered - T;
tor.mainBoxPanels = setFieldPair(pre.mainBoxPanels, qt(1));
tor.noseCellPanels = setFieldPair(pre.noseCellPanels, qt(2));
tor.physicalPanels = physicalOut;
tor.physicalPanelCount = numel(physicalOut);
tor.expectedPhysicalPanelCount = pre.expectedPhysicalPanelCount;
tor.sharedWallPhysical.panelName = "frontSpar_physical_shared";
tor.sharedWallPhysical.qt_Y_per_lbin = qtShared;
tor.sharedWallPhysical.signConventionNote = "Physical shared-wall torsion: qt_cell1 - qt_cell2.";
tor.method = "Pure My=1 torsion uses the two-cell torsion solution; no Vx/Vz moment-arm corrections are included.";
tor.units.qt = "lb/in per lb-in";
tor.units.torque = "lb-in";
tor.exclusions = ["qtot_X_per_lb"; "qtot_Z_per_lb"; "full load-case shear-flow scaling"; ...
    "moment-arm torque corrections from Vx/Vz"; "shear stress"; "von Mises"; "buckling"; ...
    "yield"; "divergence"; "weight"; "optimisation"];
end

function walls = setFieldPair(walls, q)
for k = 1:numel(walls)
    walls(k).qt_Y_per_lbin = q;
    walls(k).qtSource = "pure_My_unit_torsion_only";
end
end

%% ------------------------------------------------- combined unit flows
function ut = unitTotals(section, c, pre)
% Vz acts at the quarter chord on the chord line, Vx along the chord line
% and My is the aerodynamic pitching moment. Their torque about the shear
% centre is carried by the unit torsion flows.
sc = section.shearCenter;
tor = section.unitShear.stage3C.torsion;
xRef = 0.25 * c;
zRef = 0.0;
armZ = sc.x_sc_in - xRef;
armX = -sc.z_sc_in;

walls = tor.physicalPanels;
for k = 1:numel(walls)
    walls(k).qt_X_per_lb = armX * walls(k).qt_Y_per_lbin;
    walls(k).qt_Z_per_lb = armZ * walls(k).qt_Y_per_lbin;
    walls(k).qtot_X_per_lb = walls(k).qcorr_X_per_lb + walls(k).qt_X_per_lb;
    walls(k).qtot_Z_per_lb = walls(k).qcorr_Z_per_lb + walls(k).qt_Z_per_lb;
    walls(k).qtot_Y_per_lbin = walls(k).qt_Y_per_lbin;
    walls(k).qtotSource = "unit_coefficients_only_quarter_chord_reference_no_load_case_scaling";
end
[~, iX] = max(abs([walls.qtot_X_per_lb]));
[~, iZ] = max(abs([walls.qtot_Z_per_lb]));
[~, iY] = max(abs([walls.qtot_Y_per_lbin]));

ut = struct();
ut.loadReference.x_load_ref_in = xRef;
ut.loadReference.z_load_ref_in = zRef;
ut.loadReference.source = "Load reference convention: aerodynamic load reference at quarter chord, z = chord line.";
ut.torqueArms.torque_X_per_lb = armX;
ut.torqueArms.torque_Z_per_lb = armZ;
ut.torqueArms.torque_Y_per_lbin = 1.0;
ut.torqueArms.source = "Load reference convention: torque_Z = Vz*(x_sc - 0.25*c), torque_X = -Vx*z_sc, torque_Y = My.";
ut.physicalPanels = walls;
ut.physicalPanelCount = numel(walls);
ut.expectedPhysicalPanelCount = pre.expectedPhysicalPanelCount;
ut.maxAbsQtot_X_per_lb = abs(walls(iX).qtot_X_per_lb);
ut.maxAbsQtot_X_location = walls(iX);
ut.maxAbsQtot_Z_per_lb = abs(walls(iZ).qtot_Z_per_lb);
ut.maxAbsQtot_Z_location = walls(iZ);
ut.maxAbsQtot_Y_per_lbin = abs(walls(iY).qtot_Y_per_lbin);
ut.maxAbsQtot_Y_location = walls(iY);
ut.method = "Combined unit coefficients only: qtot_X/Z = qcorr_X/Z + torque_X/Z*qt_Y; qtot_Y = qt_Y.";
ut.units.qtot_X = "lb/in per lb";
ut.units.qtot_Z = "lb/in per lb";
ut.units.qtot_Y = "lb/in per lb-in";
ut.units.torqueArm = "in";
ut.exclusions = ["nine-case shear-flow scaling"; "shear stress"; "von Mises"; "buckling"; ...
    "yield"; "divergence"; "weight"; "optimisation"];
end

%% ------------------------------------------------------------- checks
function checkInputs(wing, opts)
if ~isstruct(wing)
    error('computeSectionUnitLoads:MissingWing', 'Input wing must be a struct from setGeometry.');
end
need = ["planform.chord_in", "box.frontSpar_xc", "box.rearSpar_xc", "sparCaps.areas_in2", ...
    "stringers.numTop", "stringers.numBottom", "stringers.numNoseTop", "stringers.numNoseBottom", ...
    "stringers.areaTop_in2", "stringers.areaBottom_in2", "stringers.areaNoseTop_in2", ...
    "stringers.areaNoseBottom_in2", "thickness.t1_frontSpar_in", "thickness.t2_rearSpar_in", ...
    "thickness.t3_upperMainSkin_in", "thickness.t4_lowerMainSkin_in", ...
    "thickness.upperNoseSkin_in", "thickness.lowerNoseSkin_in"];
for f = need
    if ~hasPath(wing, f)
        error('computeSectionUnitLoads:MissingField', ...
            'wing.%s is required for section unit-load geometry setup.', f);
    end
end
c = wing.planform.chord_in;
xs = [wing.box.frontSpar_xc, wing.box.rearSpar_xc];
counts = [wing.stringers.numTop, wing.stringers.numBottom, wing.stringers.numNoseTop, wing.stringers.numNoseBottom];
areas = [wing.stringers.areaTop_in2, wing.stringers.areaBottom_in2, ...
    wing.stringers.areaNoseTop_in2, wing.stringers.areaNoseBottom_in2];
t = wing.thickness;
thick = [t.t1_frontSpar_in, t.t2_rearSpar_in, t.t3_upperMainSkin_in, ...
    t.t4_lowerMainSkin_in, t.upperNoseSkin_in, t.lowerNoseSkin_in];
assertTrue(isfinite(c) && c > 0, 'InvalidChord', 'wing.planform.chord_in must be positive and finite.');
assertTrue(all(isfinite(xs)) && xs(1) > 0 && xs(1) < xs(2) && xs(2) < 1, ...
    'InvalidSparLocations', 'Expected 0 < frontSpar_xc < rearSpar_xc < 1.');
assertTrue(numel(wing.sparCaps.areas_in2) == 4 && all(isfinite(wing.sparCaps.areas_in2)) && ...
    all(wing.sparCaps.areas_in2 > 0), 'InvalidSparCapAreas', ...
    'wing.sparCaps.areas_in2 must contain four positive finite areas.');
assertTrue(all(isfinite(counts)) && all(counts >= 0) && all(mod(counts, 1) == 0), ...
    'InvalidStringerCounts', 'Stringer counts must be nonnegative finite integers.');
assertTrue(all(isfinite(areas)) && all(areas > 0), 'InvalidStringerAreas', ...
    'Stringer areas must be positive and finite.');
assertTrue(all(isfinite(thick)) && all(thick > 0), 'InvalidThicknesses', ...
    'Required web/panel thicknesses must be positive and finite.');

assertTrue(isstruct(opts), 'InvalidAnalysisOptions', 'analysisOptions must be a struct.');
assertTrue(isfield(opts, 'usePanelIdealisation'), 'MissingPanelIdealisationOption', ...
    'analysisOptions.usePanelIdealisation is required.');
assertTrue(islogical(opts.usePanelIdealisation) || (isnumeric(opts.usePanelIdealisation) && ...
    isscalar(opts.usePanelIdealisation)), 'InvalidPanelIdealisationOption', ...
    'analysisOptions.usePanelIdealisation must be a scalar logical value.');
assertTrue(isfield(opts, 'panelIdealisationMethod') && ...
    strlength(string(opts.panelIdealisationMethod)) > 0, 'MissingPanelIdealisationMethod', ...
    'analysisOptions.panelIdealisationMethod is required.');
end

function checkSection(section, wing)
B = boomTable(section);
assertTrue(all(isfinite([B.x; B.z])), 'InvalidCoordinates', 'All generated boom coordinates must be finite.');
assertTrue(all(isfinite(B.area)) && all(B.area > 0), 'InvalidGeneratedAreas', ...
    'All spar cap and stringer areas must be positive and finite.');
assertTrue(section.centroid.x_in >= 0 && section.centroid.x_in <= wing.planform.chord_in, ...
    'CentroidOutOfRange', 'Boom-only centroid x-coordinate must lie within the chord range.');
assertTrue(section.inertia.Ix_in4 > 0 && isfinite(section.inertia.Ix_in4), 'InvalidIx', 'Ix must be positive and finite.');
assertTrue(section.inertia.Iz_in4 > 0 && isfinite(section.inertia.Iz_in4), 'InvalidIz', 'Iz must be positive and finite.');
assertTrue(numel(section.unitBending.entries) == numel(B.x) && ...
    all(isfinite([section.unitBending.entries.sigY_Mx_per_lbin, section.unitBending.entries.sigY_Mz_per_lbin])), ...
    'InvalidBendingCoefficients', 'Every boom needs finite unit bending coefficients.');

s3A = section.unitShear.stage3A;
p = s3A.panels;
assertTrue(s3A.panelCount == s3A.expectedPanelCount, 'Stage3APanelCountChanged', ...
    'Every wall must appear in exactly one cell loop (shared wall once per loop).');
assertTrue(all([p.ds_in] > 0) && all(isfinite([p.ds_in, p.ds_over_t, p.qPrime_X_per_lb, p.qPrime_Z_per_lb])), ...
    'InvalidOpenFlows', 'Wall lengths and open-section flows must be finite (lengths positive).');

cmp = section.unitShear.stage3B.compatibility;
assertTrue(all(isfinite(cmp.Amat), 'all') && abs(det(cmp.Amat)) >= 1e-12, ...
    'SingularStage3BAmat', 'Cell compatibility matrix must be finite and non-singular.');
assertTrue(all(isfinite([cmp.qs_X_per_lb; cmp.qs_Z_per_lb])), 'InvalidStage3BQs', ...
    'Cell constant flows must be finite.');
tol = 1e-10 * max(1.0, max(abs([cmp.sharedWall.mainBox_ds_over_t, cmp.sharedWall.noseCell_ds_over_t])));
assertTrue(cmp.sharedWall.absDifference <= tol, 'Stage3BSharedWallMismatch', ...
    'Both loop copies of the shared wall must have the same ds/t.');

pre = section.unitShear.stage3B.preTorsion;
assertTrue(pre.area.Area1_mainCell_in2 > 0 && pre.area.Area2_noseCell_in2 > 0, ...
    'InvalidPreTorsionAreas', 'Cell areas must be positive.');
assertTrue(isfinite(section.shearCenter.x_sc_in) && isfinite(section.shearCenter.z_sc_in), ...
    'InvalidShearCenter', 'Shear centre coordinates must be finite.');

tor = section.unitShear.stage3C.torsion;
assertTrue(all(isfinite(tor.qt_Y_per_lbin)) && tor.q1t_over_q2t > 0 && ...
    abs(tor.recoveredTorque_lbin - tor.My_unit_lbin) < 1e-9, 'Stage3CTorqueMismatch', ...
    'Unit torsion flows must be finite and reproduce the applied torque.');

ut = section.unitShear.stage3C.unitTotals;
w = ut.physicalPanels;
assertTrue(ut.physicalPanelCount == ut.expectedPhysicalPanelCount && ...
    nnz([w.group] == "frontSparSharedPhysical") == 1 && ...
    ~any([w.group] == "frontSpar" | [w.group] == "frontSparCell2"), ...
    'UnitTotalsSharedWallDoubleCounted', 'The shared front spar must appear exactly once.');
assertTrue(all(isfinite([w.qt_X_per_lb, w.qt_Z_per_lb, w.qtot_X_per_lb, w.qtot_Z_per_lb, w.qtot_Y_per_lbin])), ...
    'InvalidUnitTotalPanelValues', 'All combined unit shear-flow coefficients must be finite.');
end

function tf = hasPath(s, dotted)
tf = true;
for part = split(string(dotted), ".")'
    if ~isstruct(s) || ~isfield(s, part)
        tf = false;
        return;
    end
    s = s.(part);
end
end

function assertTrue(cond, id, msg)
if ~cond
    error(['computeSectionUnitLoads:' id], '%s', msg);
end
end
