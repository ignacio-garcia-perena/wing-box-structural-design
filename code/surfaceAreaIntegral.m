function a = surfaceAreaIntegral(xStart, xEnd, isUpper)
% surfaceAreaIntegral Area between the chord line and the NACA 2412 surface.
%
%   a = surfaceAreaIntegral(xStart, xEnd, isUpper) integrates z/c of the
%   upper (isUpper = 1) or lower surface from xStart to xEnd (x/c), in
%   chord^2 units. For the lower surface the magnitude of the integral is
%   returned, so the result is always a positive enclosed area.
%
% Method: exact integration of the NACA four-digit mean line and thickness
% polynomials (see nacaSurfaceOrdinate for the definition). The mean-line
% antiderivative is evaluated piecewise about the maximum-camber station p.

m = 0.02;
p = 0.4;
t = 0.12;

% Antiderivatives of the two mean-line arcs.
fwdArc = @(x) m / p^2 * (p * x.^2 - x.^3 / 3);
aftArc = @(x) m / (1 - p)^2 * ((1 - 2 * p) * x + p * x.^2 - x.^3 / 3);

if xEnd < p
    meanLineArea = fwdArc(xEnd) - fwdArc(xStart);
elseif xStart < p
    meanLineArea = (fwdArc(p) - fwdArc(xStart)) + (aftArc(xEnd) - aftArc(p));
else
    meanLineArea = aftArc(xEnd) - aftArc(xStart);
end

% Antiderivative of the half-thickness distribution.
thicknessArea = @(x) t / 0.2 * (0.2969 * x.^1.5 / 1.5 - 0.1260 * x.^2 / 2 ...
    - 0.3516 * x.^3 / 3 + 0.2843 * x.^4 / 4 - 0.1015 * x.^5 / 5);
halfThicknessArea = thicknessArea(xEnd) - thicknessArea(xStart);

if isUpper == 1
    a = meanLineArea + halfThicknessArea;
else
    a = abs(meanLineArea - halfThicknessArea);
end
end
