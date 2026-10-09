function design = baselineDesign()
% baselineDesign Return the deterministic baseline wing-box design.
%
% Outputs:
%   design - Struct containing the current baseline design choices in an
%            optimisation-ready form.
%
% Scope:
%   This function only packages the baseline assumptions. It does not run
%   analysis and does not modify geometry generation.

design = struct();

design.label = "default_final_design_template";
design.status = "baseline_design_variables_packaged_for_optimisation";

design.box.frontSpar_xc = 0.15;
design.box.rearSpar_xc = 0.75;

design.thickness.t1_frontSpar_in = 0.04;
design.thickness.t2_rearSpar_in = 0.04;
design.thickness.t3_upperMainSkin_in = 0.02;
design.thickness.t4_lowerMainSkin_in = 0.02;
design.thickness.t5_noseSkin_in = 0.02;
design.thickness.upperNoseSkin_in = 0.02;
design.thickness.lowerNoseSkin_in = 0.02;

design.stringers.section_key = "A3_10_Z_40";
design.stringers.section_type = "Z-section";
design.stringers.sec_no = 40;
design.stringers.numTop = 8;
design.stringers.numBottom = 8;
design.stringers.numNoseTop = 3;
design.stringers.numNoseBottom = 3;

design.sparCaps.section_type = "T-section";
design.sparCaps.sec_no = 4;

design.material.skinWeb.alloy = "2024";
design.material.skinWeb.temper = "T3";
design.material.skinWeb.product_form = "Sheet_0.010-0.128in";

design.material.sparCapStringer.alloy = "7075";
design.material.sparCapStringer.temper = "T6";
design.material.sparCapStringer.product_form = "Extrusion";

design.ribs.ribSpacing_in = 10;
design.ribs.weightModel = "geometric_airfoil_area";
design.ribs.weightPerRib_lb = 1.0;
design.ribs.material.alloy = "2024";
design.ribs.material.temper = "T3";
design.ribs.material.product_form = "Sheet_0.010-0.128in";
design.ribs.airfoil = "NACA2412";
design.ribs.normalizedArea_A_over_c2 = 0.082169;
design.ribs.chord_in = 60.0;
design.ribs.grossArea_in2 = 295.8093;
design.ribs.solidAreaFactor = 0.7;
design.ribs.thickness_in = 0.040;
design.ribs.countingConvention = ...
    "semispan_includes_root_and_tip_fullspan_counts_root_once";
design.ribs.notes = ...
    "Geometric rib weight estimate using NACA 2412 enclosed area, 2024-T3 sheet density, 0.040 in rib thickness, and 0.7 solid-area factor.";
% Reference analyses used K_column = 1.0. The final-design Euler
% column-buckling default uses K_column = 0.7.
design.columnBuckling.K_column = 0.7;

design.warnings = [
    "ribSpacing_in = 10 in is the default baseline rib spacing."
    "Rib weight uses a geometric NACA 2412 enclosed-area estimate; constant_per_rib remains available as an alternate model."
    "K_column = 0.7 is the final-design default; K_column = 1.0 remains available as an alternate reference value."
    "setGeometry uses the deterministic baseline unless optimisation variables are supplied."
    "7075-T6 spar cap/stringer material uses the active Extrusion row from material_properties.csv."
    "2024-T3 skin/web values use the active Sheet_0.010-0.128in row from material_properties.csv."
    ];

design.exclusions = [
    "panel buckling"
    "crippling"
    "combined R-ratios"
    "final margins"
    "Monte Carlo"
    "optimisation"
    ];
end

