%% naca2412_area.m
% Computes the cross-sectional area of the NACA 2412 airfoil
% using numerical integration (trapezoidal rule) on the official
% UIUC coordinate database geometry.
%
% Coordinates source: UIUC Airfoil Coordinates Database
%   https://m-selig.ae.illinois.edu/ads/coord/naca2412.dat
%   (same data served by airfoiltools.com)
%
% Convention:
%   - Coordinates are normalized: chord c = 1
%   - Upper surface listed first (LE to TE), then lower surface (LE to TE)
%   - x runs 0 (leading edge) to 1 (trailing edge)
%   - y_upper > 0, y_lower < 0
%
% Units:
%   - Normalized area:  [chord^2], dimensionless ratio A/c^2
%   - Scaled to c = 5 ft  --> area in ft^2
%   - Scaled to c = 60 in --> area in in^2

clear; clc;

%% --- Airfoil coordinates (normalized, c = 1) --------------------------
% Format: [x, y]  -- upper surface, leading edge to trailing edge
upper = [
    0.0000    0.0000
    0.0125    0.0215
    0.0250    0.0299
    0.0500    0.0413
    0.0750    0.0496
    0.1000    0.0563
    0.1500    0.0661
    0.2000    0.0726
    0.2500    0.0767
    0.3000    0.0788
    0.4000    0.0780
    0.5000    0.0724
    0.6000    0.0636
    0.7000    0.0518
    0.8000    0.0375
    0.9000    0.0208
    0.9500    0.0114
    1.0000    0.0013
];

% Lower surface, leading edge to trailing edge
lower = [
    0.0000    0.0000
    0.0125   -0.0165
    0.0250   -0.0227
    0.0500   -0.0301
    0.0750   -0.0346
    0.1000   -0.0375
    0.1500   -0.0410
    0.2000   -0.0423
    0.2500   -0.0422
    0.3000   -0.0412
    0.4000   -0.0380
    0.5000   -0.0334
    0.6000   -0.0276
    0.7000   -0.0214
    0.8000   -0.0150
    0.9000   -0.0082
    0.9500   -0.0048
    1.0000   -0.0013
];

%% --- Interpolate onto a fine common x grid ----------------------------
% Use 10000 points for high-precision integration
x_fine = linspace(0, 1, 10000)';

y_upper_fine = interp1(upper(:,1), upper(:,2), x_fine, 'pchip');
y_lower_fine = interp1(lower(:,1), lower(:,2), x_fine, 'pchip');

%% --- Compute local thickness distribution ----------------------------
thickness = y_upper_fine - y_lower_fine;   % local height at each x

%% --- Integrate area using trapezoidal rule ----------------------------
% A = integral from 0 to c of [y_upper(x) - y_lower(x)] dx
% For normalized chord (c=1): A_norm = integral of thickness dx

A_norm = trapz(x_fine, thickness);   % dimensionless: A / c^2

%% --- Scale to actual chord -------------------------------------------
c_ft = 5;          % chord in feet
c_in = c_ft * 12;  % chord in inches = 60 in

A_ft2 = A_norm * c_ft^2;
A_in2 = A_norm * c_in^2;

%% --- Display results --------------------------------------------------
fprintf('============================================\n');
fprintf('  NACA 2412 Cross-Section Area\n');
fprintf('  Source: UIUC Coordinate Database\n');
fprintf('  Integration: pchip interp + trapz, N=10000\n');
fprintf('============================================\n');
fprintf('  Normalized area  A/c^2 = %.6f  [-]\n',   A_norm);
fprintf('  Chord            c     = %.4f  ft  = %.4f  in\n', c_ft, c_in);
fprintf('  Area             A     = %.4f  ft^2\n',  A_ft2);
fprintf('  Area             A     = %.4f  in^2\n',  A_in2);
fprintf('============================================\n');

%% --- Verification checks ---------------------------------------------
fprintf('\n--- Verification ---\n');

% 1. Max thickness should be ~12% at x ~ 0.30c
[t_max, idx_max] = max(thickness);
fprintf('  Max thickness / c = %.4f  (expected ~0.12)\n', t_max);
fprintf('  Location x/c      = %.4f  (expected ~0.30)\n', x_fine(idx_max));

% 2. Upper surface at x=0.3 should be ~+0.0788c
y_up_30 = interp1(upper(:,1), upper(:,2), 0.3, 'pchip');
fprintf('  y_upper at x/c=0.3 = %.4f  (expected ~+0.0788)\n', y_up_30);

% 3. Area should be positive
assert(A_norm > 0, 'ERROR: computed area is non-positive');
fprintf('  Area > 0 check: PASSED\n');

% 4. Rough sanity: NACA 4-digit ~0.06 to 0.09 * c^2
assert(A_norm > 0.06 && A_norm < 0.10, ...
    'ERROR: area outside expected range [0.06, 0.10]');
fprintf('  Area in expected range [0.06, 0.10]: PASSED\n');

%% --- Plot -------------------------------------------------------------
figure;
fill([upper(:,1); flipud(lower(:,1))], ...
     [upper(:,2); flipud(lower(:,2))], ...
     [0.7 0.85 1.0], 'EdgeColor', 'b', 'LineWidth', 1.5);
hold on;
plot(upper(:,1), upper(:,2), 'b.-', 'MarkerSize', 10);
plot(lower(:,1), lower(:,2), 'r.-', 'MarkerSize', 10);
axis equal;
grid on;
xlabel('x/c');
ylabel('y/c');
title(sprintf('NACA 2412  |  A = %.4f in^2  (c = %d in)', A_in2, c_in));
legend('Cross-section area', 'Upper surface', 'Lower surface', ...
       'Location', 'northeast');
