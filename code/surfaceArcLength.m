function s = surfaceArcLength(xStart, xEnd, isUpper)
% surfaceArcLength Arc length of the NACA 2412 surface between two stations.
%
%   s = surfaceArcLength(xStart, xEnd, isUpper) returns the length (in
%   chord units) of the upper (isUpper = 1) or lower surface between the
%   chordwise stations xStart and xEnd (both x/c).
%
% Method: the surface is replaced by a polyline through 11 equally spaced
% stations (10 straight segments) and the segment lengths are summed.
% Stations that fall at or ahead of the leading edge (x/c <= 0) are placed
% on the chord line, z = 0.

nSegments = 10;
h = (xEnd - xStart) / nSegments;

s = 0;
xPrev = xStart;
zPrev = nacaSurfaceOrdinate(xStart, isUpper);
for k = 1:nSegments
    xNext = xPrev + h;
    if xNext > 0
        zNext = nacaSurfaceOrdinate(xNext, isUpper);
    else
        zNext = 0;
    end
    s = s + hypot(h, zNext - zPrev);
    xPrev = xNext;
    zPrev = zNext;
end
end
