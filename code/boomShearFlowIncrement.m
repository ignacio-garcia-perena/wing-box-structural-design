function dq = boomShearFlowIncrement(x, z, Vx, Vz, Ixx, Izz, Ixz, area)
% boomShearFlowIncrement Change in open-section shear flow across a boom.
%
%   dq = boomShearFlowIncrement(x, z, Vx, Vz, Ixx, Izz, Ixz, area) returns
%   the jump in shear flow produced by a boom of the given area located
%   (x, z) from the centroid, for shear forces Vx and Vz.
%
% Method: open-section shear flow in an idealised (boom + panel) section
% with unsymmetric bending (Megson, "Aircraft Structures for Engineering
% Students", shear of idealised open and closed sections): each boom changes the shear flow by its area
% times the rate of change of its direct stress, written in the sign
% convention used by this project's section model.

detI = Ixx * Izz - Ixz^2;
dq = -area * (x * (Ixx * Vx - Ixz * Vz) + z * (Izz * Vz - Ixz * Vx)) / detI;
end
