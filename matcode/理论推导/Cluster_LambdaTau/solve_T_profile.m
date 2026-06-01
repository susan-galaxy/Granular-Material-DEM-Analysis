function [z, Tnorm, info] = solve_T_profile(phi, en, d, L, varargin)
% Steady-state granular-gas temperature profile T(z)/T_wall in a 1D slab
% under microgravity (n = P/T = const along z), from the dimensionless
% equation derived in the PPT:
%
%       d^2 u / d xi^2  =  1 / (2 u)        with u = sqrt(T/T_wall),
%       u(0) = u(Lbar) = 1,   xi = z / lambda(phi, en, d).
%
% This is solved here by the analytic first integral
%
%       u'^2  =  s^2 + ln(u)     ;   u_min = exp(-s^2)
%       Lbar/2 = int_{u_min}^1  du / sqrt(s^2 + ln u)
%
% so we 1-D root-find for s (i.e. the slope at the wall) such that the
% integral equals Lbar/2. There is a maximum Lbar_max ~ 1.78 above which
% no smooth stationary solution exists - that is the clustering collapse;
% in that case we fall back to the linearized cosh profile so the plotted
% curve stays smooth and monotone.
%
% Inputs:
%   phi  : volume fraction (scalar)
%   en   : restitution coefficient
%   d    : particle diameter
%   L    : system size
%   Name-Value:
%     'N'      : number of grid points  (default 121)
%     'Cscale' : prefactor passed to theory_lambda  (default 1)
%     'mode'   : 'nonlinear' (default) or 'linear' (force cosh form)
%
% Outputs:
%   z      : grid (length N)
%   Tnorm  : T(z)/T_wall on grid
%   info   : struct with fields lambda, Lbar, u_min, mode_used, collapsed

    p = inputParser;
    p.addParameter('N', 121);
    p.addParameter('Cscale', 1.0);
    p.addParameter('mode', 'nonlinear');
    p.parse(varargin{:});

    lambda = theory_lambda(phi, en, d, 'Cscale', p.Results.Cscale);
    Lbar   = L / lambda;
    halfLbar = Lbar / 2;

    info = struct('lambda', lambda, 'Lbar', Lbar, ...
                  'u_min', NaN, 'mode_used', '', 'collapsed', false);

    z = linspace(0, L, p.Results.N);

    if strcmpi(p.Results.mode, 'linear')
        Tnorm = cosh( (z - L/2)/lambda ) / cosh(L/(2*lambda));
        info.mode_used = 'linear-cosh';
        info.u_min = sqrt(min(Tnorm));
        return;
    end

    % ----- find Lbar_max numerically (it's a universal constant) ------
    % half-domain length as function of s :
    %   half(s) = int_{u_min(s)}^1 du / sqrt(s^2 + ln u),   u_min = exp(-s^2)
    half_of_s = @(s) integral( @(u) 1./sqrt( max(s.^2 + log(u), 1e-12) ), ...
                               exp(-s.^2), 1, 'AbsTol',1e-7,'RelTol',1e-5 );

    % bracket the maximum by sampling
    s_grid = linspace(0.02, 4.0, 400);
    h_grid = zeros(size(s_grid));
    for i = 1:numel(s_grid), h_grid(i) = half_of_s(s_grid(i)); end
    [hMax, iMax] = max(h_grid);
    halfMax = hMax;

    if halfLbar > halfMax * 0.999
        % collapsed - no smooth solution. Fall back to cosh form.
        Tnorm = cosh( (z - L/2)/lambda ) / cosh(L/(2*lambda));
        info.mode_used = 'collapsed-cosh-fallback';
        info.collapsed = true;
        info.u_min = sqrt(min(Tnorm));
        return;
    end

    % halfLbar lies on the increasing branch (s small) OR decreasing branch
    % (s large). Physically the increasing branch is the "smooth" gas
    % solution; the decreasing branch is the "compacted" branch. Pick
    % increasing branch.
    s_lo = s_grid(1);     h_lo = h_grid(1);
    s_hi = s_grid(iMax);  h_hi = h_grid(iMax);
    if halfLbar < h_lo
        % very short box -- u_min very close to 1, use linear cosh
        Tnorm = cosh( (z - L/2)/lambda ) / cosh(L/(2*lambda));
        info.mode_used = 'small-box-cosh';
        info.u_min = sqrt(min(Tnorm));
        return;
    end

    % bisection on the increasing branch
    f = @(s) half_of_s(s) - halfLbar;
    try
        s_star = fzero(f, [s_lo, s_hi]);
    catch
        Tnorm = cosh( (z - L/2)/lambda ) / cosh(L/(2*lambda));
        info.mode_used = 'fzero-failed-cosh';
        info.u_min = sqrt(min(Tnorm));
        return;
    end

    u_min = exp(-s_star^2);
    info.u_min = u_min;

    % --- reconstruct u(xi) on half-domain by inverse-integrating ----
    u_grid = linspace(u_min*1.001, 1.0, p.Results.N);
    xi_grid = zeros(size(u_grid));
    for i = 2:numel(u_grid)
        xi_grid(i) = integral( @(u) 1./sqrt(s_star^2 + log(u)), ...
                               u_grid(i), 1, 'AbsTol',1e-7, 'RelTol',1e-5 );
    end
    % xi=0 at wall (u=1), xi=halfLbar at center (u=u_min).
    % Sort by xi:
    [xi_sorted, idx] = sort(xi_grid);
    u_sorted = u_grid(idx);

    % map xi to physical z on left half : z = xi * lambda
    z_left  = xi_sorted * lambda;
    Tn_left = u_sorted.^2;
    % mirror to right half by symmetry
    z_right = L - z_left;
    [z_full, ord] = sort([z_left, z_right]);
    Tn_full = [Tn_left, Tn_left];
    Tn_full = Tn_full(ord);
    % remove dupes
    [z_full, idu] = unique(z_full);
    Tn_full = Tn_full(idu);

    z = linspace(0, L, p.Results.N);
    Tnorm = interp1(z_full, Tn_full, z, 'linear', 'extrap');
    Tnorm = max(Tnorm, 1e-6);
    info.mode_used = 'nonlinear-FirstIntegral';
end
