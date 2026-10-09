function sigma = unsymmetricBendingStress(x, z, Mx, Mz, Ixx, Izz, Ixz)
% unsymmetricBendingStress Direct stress from bending about two axes.
%
%   sigma = unsymmetricBendingStress(x, z, Mx, Mz, Ixx, Izz, Ixz) returns
%   the axial stress at a point located (x, z) from the section centroid,
%   for bending moments Mx and Mz and centroidal second moments Ixx, Izz
%   and product Ixz.
%
% Method: generalised flexure formula for a section with a non-zero
% product of inertia (Megson, "Aircraft Structures for Engineering
% Students"; Bruhn, "Analysis and Design of Flight Vehicle
% Structures"), written in the sign convention used by this
% project's section model.

detI = Ixx * Izz - Ixz^2;
sigma = -(x * (Ixx * Mz - Ixz * Mx) + z * (Izz * Mx - Ixz * Mz)) / detI;
end
