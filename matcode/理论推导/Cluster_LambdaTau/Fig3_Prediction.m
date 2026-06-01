function Fig3_Prediction(results, opt)
% Figure 3 :  Prediction figure
%
%   Panel A : T_center / T_wall  vs phi  (theory + DEM markers + critical phi)
%   Panel B : full T(z)/T_wall profiles for several phi values
%             (theory solid line + DEM dots ; demonstrates that knowing
%             lambda and tau lets us predict the spatial state and the
%             onset of clustering.)

    if nargin < 2, opt = struct(); end
    if ~isfield(opt, 'L'),        opt.L = 4.0; end
    if ~isfield(opt, 'en'),       opt.en = 0.94; end
    if ~isfield(opt, 'Cscale'),   opt.Cscale = 1.0; end
    if ~isfield(opt, 'savepath'), opt.savepath = pwd; end
    if ~isfield(opt, 'd')
        if ~isempty(results), opt.d = median([results.d]); else, opt.d = 0.08; end
    end

    %% Panel A:  Tc/Tw vs phi ---------------------------------------------
    % Plot BOTH the linearized cosh curve (gentle) and the full nonlinear
    % ODE prediction (sharp drop at collapse). The latter is the actual
    % theoretical prediction; the former is shown for reference.
    phi_th = logspace(-3, log10(0.30), 400);
    lam_th = theory_lambda(phi_th, opt.en, opt.d, 'Cscale', opt.Cscale);
    Tratio_lin = 1 ./ cosh( opt.L ./ (2*lam_th) );
    Tratio_lin = max(Tratio_lin, 1e-6);

    % Nonlinear prediction : solve the granular ODE explicitly
    Tratio_nl = NaN(size(phi_th));
    Lbar_max_half = nonlinear_max_halfLbar();    % universal constant
    for k = 1:numel(phi_th)
        halfLbar = opt.L / (2 * lam_th(k));
        if halfLbar < Lbar_max_half
            s = nonlinear_find_s(halfLbar);
            u_min = exp(-s^2);
            Tratio_nl(k) = u_min^2;
        else
            Tratio_nl(k) = 1e-6;
        end
    end
    Tratio_nl = max(Tratio_nl, 1e-6);

    Tratio_th = Tratio_nl;   % use nonlinear as the headline curve

    % critical phi : Tc/Tw < 0.3 say (heuristic clustering threshold)
    Tcrit = 0.3;
    iC = find(Tratio_th < Tcrit, 1, 'first');
    phi_c_pred = NaN;
    if ~isempty(iC) && iC > 1
        phi_c_pred = interp1(Tratio_th(iC-1:iC), phi_th(iC-1:iC), Tcrit);
    end

    % DEM
    phi_d  = [results.phi];
    Tw_d   = [results.T_wall];
    Tc_d   = [results.T_center];
    cluster = [results.cluster_flag];
    valid = ~isnan(phi_d) & ~isnan(Tw_d) & ~isnan(Tc_d) & Tw_d > 0;
    phi_d = phi_d(valid);
    Tratio_d = Tc_d(valid) ./ Tw_d(valid);
    cluster = cluster(valid);

    %% Plot ----------------------------------------------------------------
    f = figure('Color','w','Position', [120 80 1500 700]);

    axA = subplot(1, 2, 1, 'Parent', f); hold(axA,'on'); box(axA,'on');
    plot(axA, phi_th, Tratio_lin, '--', 'LineWidth', 1.8, ...
        'Color',[0.55 0.65 0.85], 'DisplayName','Linearised (cosh)');
    plot(axA, phi_th, Tratio_th, '-', 'LineWidth', 2.8, ...
        'Color',[0.10 0.30 0.85], 'DisplayName','Nonlinear ODE  T_c/T_w(\phi)');

    yline(axA, Tcrit, ':', 'Color', [0.6 0.6 0.6], 'LineWidth',1.5, ...
        'Label', sprintf('T_c/T_w = %.2f (cluster threshold)', Tcrit), ...
        'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    if ~isnan(phi_c_pred)
        xline(axA, phi_c_pred, ':', 'Color',[0.20 0.55 0.20], 'LineWidth',1.8, ...
            'Label', sprintf('\\phi_c^{predict} = %.4f', phi_c_pred), ...
            'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    end

    if ~isempty(phi_d)
        scatter(axA, phi_d(~cluster), Tratio_d(~cluster), 110, ...
            'MarkerEdgeColor',[0 0.45 0.74], 'MarkerFaceColor',[0.65 0.85 1.0], ...
            'LineWidth',1.2, 'DisplayName','DEM uniform');
        scatter(axA, phi_d( cluster), Tratio_d( cluster), 130, 'd', ...
            'MarkerEdgeColor',[0.85 0.1 0.1], 'MarkerFaceColor',[1 0.75 0.75], ...
            'LineWidth',1.2, 'DisplayName','DEM cluster');
    end

    set(axA,'XScale','log');
    xlabel(axA,'\phi','FontSize',14);
    ylabel(axA,'T_{center} / T_{wall}','FontSize',14);
    title (axA,'(A)  Predicted central temperature ratio','FontSize',14);
    legend(axA,'Location','southwest','FontSize',12);
    axA.FontSize = 12;
    ylim(axA,[0, 1.05]);
    xlim(axA,[1e-3, 0.3]);
    grid(axA,'on');

    %% Panel B :  T(z)/T_wall profiles for selected phi ------------------
    axB = subplot(1, 2, 2, 'Parent', f); hold(axB,'on'); box(axB,'on');

    % pick a few representative phi from DEM data (or fall back to grid)
    if numel(phi_d) >= 3
        sel_phi = unique(round(phi_d * 1e4) / 1e4);
        if numel(sel_phi) > 6
            ix = round(linspace(1, numel(sel_phi), 6));
            sel_phi = sel_phi(ix);
        end
    else
        sel_phi = [0.005 0.01 0.02 0.04 0.08];
    end
    cmap = parula(max(numel(sel_phi)+1, 6));

    handles = gobjects(0);
    legs    = strings(0);
    for k = 1:numel(sel_phi)
        phi_k = sel_phi(k);
        lam_k = theory_lambda(phi_k, opt.en, opt.d, 'Cscale', opt.Cscale);
        z_grid = linspace(0, opt.L, 121);
        Tn = cosh( (z_grid - opt.L/2)/lam_k ) / cosh(opt.L/(2*lam_k));
        h = plot(axB, z_grid/(opt.L/2), Tn, '-', 'LineWidth', 2.2, 'Color', cmap(k,:));
        handles(end+1) = h; %#ok<AGROW>
        legs(end+1)    = sprintf('\\phi = %.3g (theory)', phi_k); %#ok<AGROW>
    end

    % overlay DEM profiles
    for kk = 1:numel(results)
        r = results(kk);
        if isnan(r.phi) || isnan(r.T_wall) || r.T_wall <= 0, continue; end
        rmid = r.r_mid;
        z_eq = opt.L/2 - rmid;
        z_full = [z_eq, opt.L - z_eq];
        T_full = [r.T_of_rwall, r.T_of_rwall];
        Tnorm  = T_full / r.T_wall;
        [~, kClosest] = min(abs(sel_phi - r.phi));
        plot(axB, z_full/(opt.L/2), Tnorm, ':o', ...
            'Color', cmap(kClosest,:), 'MarkerFaceColor', cmap(kClosest,:), ...
            'MarkerSize', 5, 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    xlabel(axB, 'z / (L/2)','FontSize',14);
    ylabel(axB, 'T(z) / T_{wall}','FontSize',14);
    title (axB, '(B)  T(z) profiles : line = theory  |  dotted = DEM','FontSize',14);
    legend(axB, handles, legs, 'Location','south','FontSize', 10, 'NumColumns', 2);
    axB.FontSize = 12;
    ylim(axB,[0, 1.1]);
    xlim(axB,[0, 2]);
    grid(axB,'on');

    %% save
    if ~exist(opt.savepath,'dir'), mkdir(opt.savepath); end
    base = fullfile(opt.savepath, 'Fig3_prediction');
    saveas(f, [base '.png']);
    saveas(f, [base '.fig']);
    fprintf('\n[Fig3]  saved : %s.png / .fig\n', base);
    if ~isnan(phi_c_pred)
        fprintf('[Fig3]  predicted phi_c (T_c/T_w < %.2f) = %.4f\n', Tcrit, phi_c_pred);
    end
end

% =======================================================================
function Lh = nonlinear_max_halfLbar()
% Universal maximum half-Lbar at which the smooth solution of
% u'' = 1/(2u),  u(0)=u(Lbar)=1  still exists.
% Returns ~1.7841 (computed once and cached).

    persistent val
    if isempty(val)
        half_of_s = @(s) integral( @(u) 1./sqrt( max(s.^2 + log(u), 1e-12) ), ...
                                   exp(-s.^2), 1, 'AbsTol',1e-7,'RelTol',1e-5 );
        [~, val] = fminbnd(@(s) -half_of_s(s), 0.05, 4, ...
                           optimset('TolX',1e-6));
        val = -val;
    end
    Lh = val;
end

function s = nonlinear_find_s(halfLbar)
% Solve  half_of_s(s) = halfLbar for s on the SMOOTH branch (s small).
    half_of_s = @(s) integral( @(u) 1./sqrt( max(s.^2 + log(u), 1e-12) ), ...
                               exp(-s.^2), 1, 'AbsTol',1e-7,'RelTol',1e-5 );
    s_max_val = nonlinear_max_halfLbar();
    [~, smax] = fminbnd(@(s) -half_of_s(s), 0.05, 4, optimset('TolX',1e-6));
    smax = -smax; %#ok<NASGU>
    % Bisection on (0.001, s_at_max)
    try
        % find s value at peak first
        sgrid = linspace(0.01, 4, 200);
        hg = arrayfun(half_of_s, sgrid);
        [~, im] = max(hg);
        s_peak = sgrid(im);
        s = fzero(@(s) half_of_s(s) - halfLbar, [0.001, s_peak]);
    catch
        s = NaN;
    end
end
