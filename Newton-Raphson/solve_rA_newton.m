% Newton-Raphson solver validation
% Using the appropriate minimization method didn't provide a significant
% improvement compared to a simpler search grid method

clc, clear
[rA_out, Va_out, rP_out, a_out, e_out, mn_out] = solve_rA_newton([1 1.5 2], 63.4, 63.4, 0,'nscan',600, 'verbose', true);

%% FUNCTION
function [rA_out, Va_out, rP_out, a_out, e_out, mn_out] = ...
    solve_rA_newton(mn_targets, i_deg, phi_deg, use_special, varargin)
%SOLVE_RA_NEWTON  Newton-Raphson solver — LGO apogee radius at target m/n
%
%  Finds rA [km] such that  m_n(rA) = mn_target  for each entry in
%  mn_targets, using Newton-Raphson with bisection fallback.
%
%  Algorithm
%  ─────────
%  1. Coarse log-spaced scan (N_SCAN=30 pts) to bracket the root.
%  2. Newton-Raphson step with central-difference derivative.
%  3. If the Newton step leaves the bracket, bisect instead (Brent-style).
%  4. Repeat until |m_n(rA) - target| < tol.
%
%  ── Mandatory inputs ────────────────────────────────────────────────────
%  mn_targets  [-]    target m/n values               (scalar or row vector)
%  i_deg       [deg]  orbital inclination              (scalar)
%  phi_deg     [deg]  apogee sub-satellite latitude    (scalar)
%                       set equal to i_deg for the special LGO case
%  use_special [0/1]  1 = special LGO (phi = i),  0 = general LGO
%
%  ── Optional name-value pairs ───────────────────────────────────────────
%  'tol'      convergence tolerance on |m_n - target|   (default 1e-9)
%  'maxiter'  max Newton/bisection iterations per target (default 80)
%  'nscan'    number of points in the coarse bracket scan (default 30)
%  'verbose'  print iteration log per target (default false)
%
%  ── Outputs ─────────────────────────────────────────────────────────────
%  rA_out   [km]    apogee geocentric radius    (NaN if no convergence)
%  Va_out   [km/s]  apogee velocity
%  rP_out   [km]    perigee geocentric radius
%  a_out    [km]    semimajor axis
%  e_out    [-]     eccentricity
%  mn_out   [-]     m/n at solution             (verify ≈ mn_target)
%
%  ── Example ─────────────────────────────────────────────────────────────
%  % Special LGO, i = 63.4 deg, find orbits for m/n = 1, 1.5, 2
%  [rA, Va, rP, a, e, mn] = solve_rA_newton([1, 1.5, 2], 63.4, 63.4, 1);
%  H = rA - 6371   % [km] apogee altitude
%
%  % General LGO, i = 63.4 deg, phi = 53 deg
%  [rA, Va, rP, a, e, mn] = solve_rA_newton([1, 2], 63.4, 53.0, 0);

% ── Physical constants ────────────────────────────────────────────────────
mu    = 398600.4;       % [km^3/s^2] gravitational parameter
omE   = 7.2921159e-5;   % [s^-1]     Earth rotation rate
RE    = 6378.0;         % [km]       Earth equatorial radius
J2    = 0.0010827;      % [-]        second zonal harmonic
Rmean = 6371.0;         % [km]       Earth mean radius

rA_lo_bound = RE    + 300.0;       % [km] min valid apogee (above atmosphere)
rA_hi_bound = Rmean + 60000.0;     % [km] max valid apogee (grid limit)

% ── Parse optional arguments ──────────────────────────────────────────────
p = inputParser;
addParameter(p, 'tol',     1e-9,  @(x) isscalar(x) && x > 0);
addParameter(p, 'maxiter', 80,    @(x) isscalar(x) && x > 0);
addParameter(p, 'nscan',   30,    @(x) isscalar(x) && x >= 2);
addParameter(p, 'verbose', false, @islogical);
parse(p, varargin{:});
tol     = p.Results.tol;
maxiter = p.Results.maxiter;
nscan   = p.Results.nscan;
verbose = p.Results.verbose;

% ── Convert to radians ────────────────────────────────────────────────────
i_rad   = i_deg   * pi / 180;   % [rad] inclination
phi_rad = phi_deg * pi / 180;   % [rad] sub-satellite latitude

% gamma [rad]: angle between Va and frame velocity Vom
% cosg = cos(i)/cos(phi)  (special case: phi=i → cosg=1)
if use_special || abs(cos(phi_rad)) < 1e-12
    cosg = 1.0;
else
    cosg = cos(i_rad) / cos(phi_rad);
end

% ── Pre-allocate outputs ──────────────────────────────────────────────────
rA_out = nan(size(mn_targets));
Va_out = nan(size(mn_targets));
rP_out = nan(size(mn_targets));
a_out  = nan(size(mn_targets));
e_out  = nan(size(mn_targets));
mn_out = nan(size(mn_targets));

% ── Coarse bracket scan (shared across all targets) ──────────────────────
% Build log-spaced rA grid and evaluate m_n once
log_rA_vals = linspace(log(rA_lo_bound), log(rA_hi_bound), nscan);
rA_scan     = exp(log_rA_vals);   % [km]
mn_scan     = zeros(1, nscan);

for s = 1:nscan
    mn_scan(s) = eval_mn_lgo(rA_scan(s), i_rad, phi_rad, cosg, ...
                              use_special, mu, omE, RE, J2);
end

% ── Solve for each target ─────────────────────────────────────────────────
for t = 1:numel(mn_targets)
    mnt = mn_targets(t);
    if mnt <= 0 || ~isfinite(mnt), continue; end

    % ── Step 1: find bracket from coarse scan ─────────────────────────
    rA_lo = NaN;  rA_hi = NaN;
    for s = 1:nscan-1
        f_lo = mn_scan(s)   - mnt;
        f_hi = mn_scan(s+1) - mnt;
        if isfinite(f_lo) && isfinite(f_hi) && f_lo*f_hi <= 0
            rA_lo = rA_scan(s);
            rA_hi = rA_scan(s+1);
            break
        end
    end

    if isnan(rA_lo)
        if verbose
            fprintf('  Target m/n=%.4f: no bracket found (out of range)\n', mnt);
        end
        continue
    end

    % ── Step 2: Newton-Raphson with bisection fallback ────────────────
    rA = (rA_lo + rA_hi) / 2.0;   % start at bracket midpoint
    converged = false;

    for iter = 1:maxiter
        mn_val = eval_mn_lgo(rA, i_rad, phi_rad, cosg, ...
                              use_special, mu, omE, RE, J2);

        if ~isfinite(mn_val)
            rA = (rA_lo + rA_hi) / 2.0;   % bisect if invalid
            continue
        end

        f = mn_val - mnt;

        if abs(f) < tol
            converged = true;
            if verbose
                fprintf('  iter %2d  rA=%9.2f  m_n=%.9f  f=%+.3e\n', ...
                    iter, rA, mn_val, f);
            end
            break
        end

        % Tighten bracket for bisection fallback
        if f > 0,  rA_lo = rA;
        else,      rA_hi = rA;
        end

        % Central-difference derivative
        h    = rA * 1e-5;                         % step [km]
        mn_p = eval_mn_lgo(rA+h, i_rad, phi_rad, cosg, use_special, mu, omE, RE, J2);
        mn_m = eval_mn_lgo(rA-h, i_rad, phi_rad, cosg, use_special, mu, omE, RE, J2);

        if     isfinite(mn_p) && isfinite(mn_m),  df = (mn_p - mn_m) / (2*h);
        elseif isfinite(mn_p),                     df = (mn_p - mn_val) / h;
        elseif isfinite(mn_m),                     df = (mn_val - mn_m) / h;
        else
            rA = (rA_lo + rA_hi) / 2.0;           % fallback: bisect
            continue
        end

        if abs(df) < 1e-20
            rA = (rA_lo + rA_hi) / 2.0;           % flat region → bisect
            continue
        end

        rA_newton = rA - f / df;

        % Accept Newton step only if it stays inside the bracket
        if rA_newton > rA_lo && rA_newton < rA_hi
            rA = rA_newton;
        else
            rA = (rA_lo + rA_hi) / 2.0;           % bisect fallback
        end

        if verbose
            fprintf('  iter %2d  rA=%9.2f  m_n=%.9f  f=%+.3e\n', ...
                    iter, rA, mn_val, f);
        end
    end

    if ~converged, continue; end

    % ── Step 3: extract orbit parameters ──────────────────────────────
    [rp, a, e] = orbit_params_lgo(rA, i_rad, phi_rad, cosg, ...
                                   use_special, mu, omE);
    if isnan(rp), continue; end

    Va  = sqrt(2*mu*rp / (rA*(rA+rp)));           % [km/s] apogee velocity
    mn_sol = eval_mn_lgo(rA, i_rad, phi_rad, cosg, ...
                          use_special, mu, omE, RE, J2);

    if ~isfinite(Va) || ~isfinite(mn_sol), continue; end

    rA_out(t) = rA;
    Va_out(t) = Va;
    rP_out(t) = rp;
    a_out(t)  = a;
    e_out(t)  = e;
    mn_out(t) = mn_sol;
end

end   % ── end main function ──────────────────────────────────────────────────


%% ── eval_mn_lgo ─────────────────────────────────────────────────────────────
% Evaluates the m/n repetition factor at a given apogee radius rA [km].
% Returns NaN for physically invalid configurations.
function mn = eval_mn_lgo(rA, i_rad, phi_rad, cosg, use_special, mu, omE, RE, J2)

ci  = cos(i_rad);   ci2 = ci^2;
cp  = cos(phi_rad); cp2 = cp^2;
rA3 = rA^3;

if use_special
    den = 2*mu - omE^2 * rA3 * ci2;          % denominator: special LGO
else
    den = 2*mu*cosg - omE^2 * rA3 * cp2;     % denominator: general LGO
end

if den <= 0 || ~isfinite(den), mn = NaN; return; end

rp = omE^2 * rA^4 * ci2 / den;               % [km] LGO perigee radius

if rp <= 0 || rp >= rA, mn = NaN; return; end

a        = (rA + rp) / 2;                    % [km] semimajor axis
e        = (rA - rp) / (rA + rp);            % [-]  eccentricity
one_m_e2 = 1 - e^2;

if a <= 0 || one_m_e2 <= 0, mn = NaN; return; end

a3h  = a^1.5;
a7h  = a^3.5;
e2sq = one_m_e2^2;
sqe  = sqrt(one_m_e2);

K = 1 - (1.5*J2*RE^2 / (a^2*e2sq)) * ...
    (  (10*ci2 - 2)/4 ...
     - (3*ci2  - 1)/4 * sqe ...
     + sqrt(mu)*ci / (omE*a7h*e2sq) );

mn_val = 1 / (omE * a3h/sqrt(mu) * K);       % [-] m/n repetition factor

if isfinite(mn_val) && mn_val > 0
    mn = mn_val;
else
    mn = NaN;
end

end  % eval_mn_lgo


%% ── orbit_params_lgo ────────────────────────────────────────────────────────
% Returns rP [km], a [km], e [-] for the LGO orbit at given rA [km].
function [rp, a, e] = orbit_params_lgo(rA, i_rad, phi_rad, cosg, ...
                                        use_special, mu, omE)
ci2 = cos(i_rad)^2;
cp2 = cos(phi_rad)^2;

if use_special
    den = 2*mu - omE^2 * rA^3 * ci2;
else
    den = 2*mu*cosg - omE^2 * rA^3 * cp2;
end

if den <= 0
    rp = NaN; a = NaN; e = NaN; return
end

rp = omE^2 * rA^4 * ci2 / den;

if rp <= 0 || rp >= rA
    rp = NaN; a = NaN; e = NaN; return
end

a = (rA + rp) / 2;
e = (rA - rp) / (rA + rp);

end  % orbit_params_lgo
