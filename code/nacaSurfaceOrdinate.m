function zc = nacaSurfaceOrdinate(xc, isUpper)
% nacaSurfaceOrdinate Height of the NACA 2412 upper or lower surface.
%
%   zc = nacaSurfaceOrdinate(xc, isUpper) returns z/c at the chordwise
%   station x/c. isUpper = 1 selects the upper surface, any other value the
%   lower surface.
%
% Method: standard NACA four-digit definition (Abbott & von Doenhoff,
% "Theory of Wing Sections", 1959, Sec. 6.4; NASA TM-4741). For 2412 the
% maximum camber is m = 0.02 at p = 0.4 and the thickness is t = 0.12.
% The half-thickness is added vertically to the mean line (the common
% small-camber simplification), and the open trailing-edge coefficient
% -0.1015 is used.

if xc < 0
    warning('nacaSurfaceOrdinate:NegativeStation', ...
        'Chordwise station x/c = %g is outside the airfoil.', xc);
end

m = 0.02;      % maximum camber (fraction of chord)
p = 0.4;       % chordwise position of maximum camber
t = 0.12;      % maximum thickness (fraction of chord)

% Mean line: two parabolic arcs joined at x/c = p.
if xc < p
    meanLine = m / p^2 * (2 * p * xc - xc^2);
else
    meanLine = m / (1 - p)^2 * ((1 - 2 * p) + 2 * p * xc - xc^2);
end

% Thickness distribution, written as t/0.2 times the 20 %-thick polynomial.
thicknessPoly = [-0.1015, 0.2843, -0.3516, -0.1260];   % x^4, x^3, x^2, x^1
halfThickness = t / 0.2 * (0.2969 * sqrt(xc) + polyval([thicknessPoly 0], xc));

if isUpper == 1
    zc = meanLine + halfThickness;
else
    zc = meanLine - halfThickness;
end
end
