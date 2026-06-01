function out = dem_extract_run(matfile, varargin)
% Extract granular-gas observables from one DEM .mat file produced by
% Pos_Vel_Extrace_v3.m  (the file contains struct ExtractedData with a
% sub-field like Tag_2.Pos / .Vel / .Radius / .TimeStep).
%
% Outputs the following fields in struct 'out':
%
%   phi          : volume fraction (read from file name; NaN if absent)
%   L            : box length (cm) - assumed cubic, taken from option
%   d            : effective particle diameter (cm) - from data or option
%   m            : particle mass (g) - from option (rho, d)
%   n_mean       : time-averaged number density (1/cm^3)
%   T_mean       : time-averaged granular temperature
%                  T = m * <(v - <v>)^2> / 3      [g*cm^2/s^2]
%   T_wall       : granular temperature in a thin layer adjacent to the
%                  walls (the "boundary temperature" used in theory)
%   T_center     : granular temperature inside a central sphere
%   r_edges      : bin edges for r_wall  = min distance to any wall
%   T_of_rwall   : T(r_wall) profile, time- and shell-averaged
%   n_of_rwall   : n(r_wall) profile
%   lam_dem      : characteristic length fitted from T(r_wall) ~ cosh
%   tau_dem_mfp  : DEM mean free time from Enskog formula at (n_mean,T_mean)
%   tau_dem_H    : DEM Haff time at (n_mean,T_mean,en)
%   sigma2_n     : variance of local density across bins (cluster indicator)
%   cluster_flag : true if the system is judged "clustered"
%
%   The thoroughly time-averaged quantities use the last 'frac_late'
%   fraction of frames (default 0.5) to skip transient.

    p = inputParser;
    p.addParameter('L',       4.0);          % cm (4cm cube)
    p.addParameter('box_min', NaN);          % cm, lower-left corner; NaN = auto
    p.addParameter('d_force', 0.08);         % override diameter if data has none
    p.addParameter('rho',     6.0);          % g/cm^3 (ZrO2)
    p.addParameter('en',      0.94);         % wall-particle restitution
    p.addParameter('NBin',    25);           % # bins for r_wall profile
    p.addParameter('frac_late', 0.5);        % use last fraction of frames
    p.addParameter('TagField', '');          % e.g. 'Tag_2'; auto-detect if empty
    p.addParameter('phi_override', NaN);     % manual phi (else parsed)
    p.parse(varargin{:});

    %% 1. Load file ---------------------------------------------------------
    S = load(matfile);
    if isfield(S, 'ExtractedData')
        ED = S.ExtractedData;
    else
        error('dem_extract_run:NoExtractedData', ...
              'File %s has no "ExtractedData" struct.', matfile);
    end

    % Auto-detect Tag field
    tagField = p.Results.TagField;
    if isempty(tagField)
        fn = fieldnames(ED);
        tagField = fn{1};
    end
    DD = ED.(tagField);

    nSteps = numel(DD.Pos);
    if nSteps == 0
        error('dem_extract_run:NoFrames', 'No frames in %s.', matfile);
    end

    %% 2. Box / particle geometry ------------------------------------------
    L = p.Results.L;
    if iscell(DD.Radius) && ~isempty(DD.Radius{1})
        d = 2 * median(DD.Radius{1});
    else
        d = p.Results.d_force;
    end
    rho = p.Results.rho;
    m   = rho * (4/3) * pi * (d/2)^3;        % g

    % Auto-detect box origin if not provided ---------------------------
    box_min = p.Results.box_min;
    if isnan(box_min)
        sampleP = DD.Pos{ceil(numel(DD.Pos)/2)};
        if isempty(sampleP), sampleP = DD.Pos{1}; end
        if isempty(sampleP)
            box_min = 0;
        else
            xmin = min(min(sampleP));
            box_min = floor(xmin * 10) / 10;     % round down to 0.1cm
            if box_min < -L/4, box_min = -L/2; end
            if box_min > L/4,  box_min = 0;    end
        end
    end

    %% 3. Phi ---------------------------------------------------------------
    [~, fname, ~] = fileparts(matfile);
    if ~isnan(p.Results.phi_override)
        phi = p.Results.phi_override;
    else
        phi = phi_from_name(fname);
    end

    %% 4. Pre-allocate accumulators ----------------------------------------
    nBin = p.Results.NBin;
    r_edges = linspace(0, L/2, nBin+1);
    r_mid   = 0.5 * (r_edges(1:end-1) + r_edges(2:end));

    T_acc = zeros(1, nBin);
    n_acc = zeros(1, nBin);
    cnt   = zeros(1, nBin);

    % shell volumes (rough): cube with walls at 0 & L. The set { r_wall in
    % [a,b] } has volume = L^3 - (L-2b)^3 - [L^3 - (L-2a)^3] = (L-2a)^3 - (L-2b)^3
    shell_vol = (L - 2*r_edges(1:end-1)).^3 - (L - 2*r_edges(2:end)).^3;
    shell_vol(shell_vol <= 0) = NaN;         % central bin can blow up

    % Frame selection (transient skip)
    iStart = max(1, round((1 - p.Results.frac_late) * nSteps));
    framesUsed = iStart:nSteps;
    nFrames = numel(framesUsed);

    T_glob = zeros(1, nFrames);
    n_glob = zeros(1, nFrames);
    T_wallLayer = zeros(1, nFrames);
    T_centerLayer = zeros(1, nFrames);

    wallSlab = 0.20;                          % cm , layer thickness near wall
    centerR  = 0.40;                          % cm , radius for central T

    sigma2_n_frame = zeros(1, nFrames);

    for k = 1:nFrames
        i = framesUsed(k);
        P = DD.Pos{i};
        V = DD.Vel{i};
        if isempty(P), continue; end

        % Global granular temperature
        Vmean = mean(V, 1);
        dv2   = sum( (V - Vmean).^2, 2);      % per particle
        T_glob(k) = m * mean(dv2) / 3;        % T = m <dv^2>/3

        % Number density
        N = size(P, 1);
        n_glob(k) = N / L^3;

        % r_wall : distance to nearest wall (box [box_min, box_min + L]^3)
        Prel = P - box_min;
        rW = min( [Prel, L - Prel], [], 2 );
        rW = max(rW, 0);

        % Wall slab / center T
        maskWall   = rW < wallSlab;
        maskCenter = sqrt( sum((Prel - L/2).^2, 2) ) < centerR;

        if any(maskWall)
            vw = V(maskWall, :); vwm = mean(vw,1);
            T_wallLayer(k) = m * mean(sum( (vw - vwm).^2, 2)) / 3;
        else
            T_wallLayer(k) = NaN;
        end
        if any(maskCenter)
            vc = V(maskCenter, :); vcm = mean(vc,1);
            T_centerLayer(k) = m * mean(sum( (vc - vcm).^2, 2)) / 3;
        else
            T_centerLayer(k) = NaN;
        end

        % Binning along r_wall
        [~, ~, idx] = histcounts(rW, r_edges);
        for b = 1:nBin
            sel = (idx == b);
            if any(sel)
                vb = V(sel, :);
                vbm = mean(vb, 1);
                T_b = m * mean(sum( (vb - vbm).^2, 2)) / 3;
                T_acc(b) = T_acc(b) + T_b;
                n_acc(b) = n_acc(b) + sum(sel);
                cnt(b)   = cnt(b) + 1;
            end
        end

        % local density variance : cubic sub-cells
        nGrid = 5;
        edges = linspace(0, L, nGrid+1);
        count = local_cubic_hist(Prel, edges);
        cellVol = (L/nGrid)^3;
        nLoc = count(:) / cellVol;
        if mean(nLoc) > 0
            sigma2_n_frame(k) = var(nLoc) / mean(nLoc)^2;
        else
            sigma2_n_frame(k) = NaN;
        end
    end

    %% 5. Averages ---------------------------------------------------------
    T_of_rwall = T_acc ./ max(cnt, 1);
    n_of_rwall = n_acc ./ max(cnt, 1) ./ shell_vol;     % per-frame mean N / shell_vol

    T_mean    = mean(T_glob, 'omitnan');
    n_mean    = mean(n_glob, 'omitnan');
    T_wall    = mean(T_wallLayer,   'omitnan');
    T_center  = mean(T_centerLayer, 'omitnan');
    sigma2_n  = mean(sigma2_n_frame, 'omitnan');

    %% 6. Fit DEM characteristic length from T(r_wall) ----------------------
    %     T(r_wall) = T_wall * cosh( (L/2 - r_wall) / lam ) / cosh( L/(2*lam) )
    %  Rearranged: ln( T(r_wall)/T_wall ) ~ ln(cosh(...)) ...
    %  We do a robust 1-parameter fit for lam.
    valid = T_of_rwall > 0 & isfinite(T_of_rwall);
    lam_dem = NaN;
    if sum(valid) >= 4
        Tfit  = T_of_rwall(valid) / T_wall;
        rfit  = r_mid(valid);
        % model: M(lam) = cosh((L/2 - r)/lam) / cosh(L/(2*lam))
        modelfun = @(lam, r) cosh( (L/2 - r) ./ lam ) ./ cosh( L/(2*lam) );
        cost = @(lam) sum( ( log( max(Tfit, 1e-6) ) - log( max(modelfun(lam, rfit), 1e-6) ) ).^2 );
        try
            lam0 = max(L/8, theory_lambda(phi, p.Results.en, d));
            opts = optimset('Display','off','TolX',1e-5,'TolFun',1e-6);
            lam_dem = fminbnd(cost, d, 20*L, opts);
        catch
            lam_dem = NaN;
        end
    end

    %% 7. Time scales (DEM) -------------------------------------------------
    [tau_dem_H, tau_dem_mfp] = theory_tau(phi, p.Results.en, d, T_mean, m);

    %% 8. Cluster judgement -------------------------------------------------
    cluster_flag = (sigma2_n > 0.05) || ( (T_center / T_wall) < 0.3 );

    %% 9. Pack out ---------------------------------------------------------
    out = struct();
    out.matfile     = matfile;
    out.phi         = phi;
    out.L           = L;
    out.box_min     = box_min;
    out.d           = d;
    out.m           = m;
    out.en          = p.Results.en;
    out.nFrames     = nFrames;
    out.n_mean      = n_mean;
    out.T_mean      = T_mean;
    out.T_wall      = T_wall;
    out.T_center    = T_center;
    out.r_edges     = r_edges;
    out.r_mid       = r_mid;
    out.T_of_rwall  = T_of_rwall;
    out.n_of_rwall  = n_of_rwall;
    out.lam_dem     = lam_dem;
    out.tau_dem_H   = tau_dem_H;
    out.tau_dem_mfp = tau_dem_mfp;
    out.sigma2_n    = sigma2_n;
    out.cluster_flag = cluster_flag;
end

% ----------------------------------------------------------------------
function C = local_cubic_hist(P, edges)
% Minimal replacement for histcn : counts particles in a regular 3D grid.
    nx = numel(edges) - 1;
    C  = zeros(nx, nx, nx);
    ix = discretize(P(:,1), edges);
    iy = discretize(P(:,2), edges);
    iz = discretize(P(:,3), edges);
    ok = ~isnan(ix) & ~isnan(iy) & ~isnan(iz);
    for k = find(ok)'
        C(ix(k), iy(k), iz(k)) = C(ix(k), iy(k), iz(k)) + 1;
    end
end
