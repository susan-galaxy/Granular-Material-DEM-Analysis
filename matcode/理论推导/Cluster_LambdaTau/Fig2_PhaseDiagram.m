function Fig2_PhaseDiagram(results, opt)
% Figure 2 :  Phase diagram in the dimensionless (lambda*, tau*) plane.
%
%   x-axis : Lambda_hat = lambda / (L/2)         (energy reach ratio)
%   y-axis : Tau_hat    = tau_H  / tau_mfp       (collisions per Haff time)
%                       =  tau_H_per_dt = how many mean free times per Haff
%                       (low value -> very dissipative)
%
% Each DEM run is one point ; theory curve = locus of states with phi
% running on [phi_min, phi_max]. Background heat-map encodes the cluster
% indicator (sigma^2_n / <n>^2).
%
% A vertical dashed line at lambda* = 1 separates uniform-gas regime
% (right) from clustering regime (left), based on the theoretical
% criterion lambda < L/2.

    if nargin < 2, opt = struct(); end
    if ~isfield(opt, 'L'),        opt.L = 4.0; end
    if ~isfield(opt, 'en'),       opt.en = 0.94; end
    if ~isfield(opt, 'Cscale'),   opt.Cscale = 1.0; end
    if ~isfield(opt, 'savepath'), opt.savepath = pwd; end
    if ~isfield(opt, 'd')
        if ~isempty(results), opt.d = median([results.d]); else, opt.d = 0.08; end
    end
    if ~isfield(opt, 'm')
        if ~isempty(results), opt.m = median([results.m]); else, opt.m = 6.0*(4/3)*pi*(0.04)^3; end
    end
    if ~isfield(opt, 'T_ref')
        if ~isempty(results)
            opt.T_ref = median([results.T_mean]);
        else
            opt.T_ref = 1.0;
        end
    end

    %% Theory phi-sweep (parametric curve in (lambda*, tau*) plane) -------
    phi_th = logspace(-3.2, log10(0.30), 300);
    lam_th = theory_lambda(phi_th, opt.en, opt.d, 'Cscale', opt.Cscale);
    [tauH_th, tauMfp_th] = theory_tau(phi_th, opt.en, opt.d, opt.T_ref, opt.m);

    Lhat_th = lam_th / (opt.L/2);
    That_th = tauH_th ./ tauMfp_th;   % ~ 1/(1-en^2) , nearly constant in phi

    %% DEM points ---------------------------------------------------------
    phi_d   = [results.phi];
    lam_d   = [results.lam_dem];
    tauH_d  = [results.tau_dem_H];
    tauM_d  = [results.tau_dem_mfp];
    sig_d   = [results.sigma2_n];
    cluster = [results.cluster_flag];

    valid = ~isnan(phi_d) & ~isnan(lam_d) & ~isnan(tauH_d) & ~isnan(tauM_d) ...
            & isfinite(lam_d) & isfinite(tauH_d) & isfinite(tauM_d);
    phi_d   = phi_d(valid);
    lam_d   = lam_d(valid);
    tauH_d  = tauH_d(valid);
    tauM_d  = tauM_d(valid);
    sig_d   = sig_d(valid);
    cluster = cluster(valid);

    Lhat_d = lam_d / (opt.L/2);
    That_d = tauH_d ./ tauM_d;

    %% Build heat-map background : cluster indicator vs (phi, e_n) -------
    phi_grid = logspace(-3.2, log10(0.20), 120);
    en_grid  = linspace(0.80, 0.999, 100);
    [PG, EG] = meshgrid(phi_grid, en_grid);
    LamG = theory_lambda(PG, EG, opt.d, 'Cscale', opt.Cscale);
    LhatG = LamG / (opt.L/2);
    [tH, tM] = theory_tau(PG, EG, opt.d, opt.T_ref, opt.m);
    ThatG = tH ./ tM;

    %% Plot ---------------------------------------------------------------
    f = figure('Color','w','Position',[120 80 1500 700]);

    % ---- panel A : (phi , en) phase plane with critical curve ----------
    axA = subplot(1, 2, 1, 'Parent', f); hold(axA,'on');
    [~, hContF] = contourf(axA, PG, EG, log10(LhatG), 30, 'LineStyle','none');
    set(hContF, 'HandleVisibility', 'off');
    cb = colorbar(axA);
    cb.Label.String = 'log_{10}( \lambda / (L/2) )';
    cb.Label.FontSize = 12;
    colormap(axA, parula);

    % critical curve lambda = L/2
    [~, hCrit] = contour(axA, PG, EG, LhatG, [1 1], 'k-', 'LineWidth', 2.5);
    set(hCrit, 'DisplayName', '\lambda = L/2 (critical curve)');
    text(axA, 5e-3, 0.86, '\lambda = L/2', 'FontSize', 14, 'Color','k', ...
         'BackgroundColor', [1 1 1 0.7]);

    % DEM operating point
    if ~isempty(phi_d)
        scatter(axA, phi_d(~cluster), 0*phi_d(~cluster)+opt.en, 90, ...
            'MarkerEdgeColor','w', 'MarkerFaceColor',[0.10 0.40 0.85], ...
            'LineWidth',1.1, 'DisplayName','DEM uniform');
        scatter(axA, phi_d( cluster), 0*phi_d( cluster)+opt.en, 130, 'd', ...
            'MarkerEdgeColor','w', 'MarkerFaceColor',[0.85 0.10 0.10], ...
            'LineWidth',1.1, 'DisplayName','DEM cluster');
    end

    set(axA, 'XScale','log');
    xlabel(axA, '\phi  (volume fraction)','FontSize',14);
    ylabel(axA, 'e_n  (restitution)','FontSize',14);
    title (axA, '(A)  Phase plane (\phi , e_n) — color: log_{10}(\lambda/(L/2))','FontSize',14);
    axA.FontSize = 12;
    xlim(axA, [1e-3, 0.2]); ylim(axA, [0.80, 0.999]);
    legend(axA, 'Location','southwest');
    grid(axA,'on'); box(axA,'on');

    % ---- panel B : dimensionless (Lambda_hat , Tau_hat) trajectory ----
    axB = subplot(1, 2, 2, 'Parent', f); hold(axB,'on'); box(axB,'on');

    % shade regimes
    yl = [1, max([That_th That_d 100])*1.4];
    yl(1) = max(1, 0.5*min([That_th That_d 10]));
    patch(axB, [1e-3 1 1 1e-3], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], ...
          'EdgeColor','none','FaceAlpha',0.55,'DisplayName','Clustering');
    patch(axB, [1 1e3 1e3 1], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.95 1.0], ...
          'EdgeColor','none','FaceAlpha',0.55,'DisplayName','Uniform gas');

    % theory locus
    plot(axB, Lhat_th, That_th, '-', 'Color',[0.10 0.30 0.85], ...
        'LineWidth', 2.5, 'DisplayName','Theory locus (sweep \phi)');

    % phi tick labels along theory curve
    phi_marks = [0.005 0.01 0.02 0.04 0.08];
    for q = phi_marks
        if q >= phi_th(1) && q <= phi_th(end)
            xs = interp1(phi_th, Lhat_th, q);
            ys = interp1(phi_th, That_th, q);
            plot(axB, xs, ys, 'k.', 'MarkerSize', 14, 'HandleVisibility','off');
            text(axB, xs*1.05, ys*1.30, sprintf('\\phi=%.3g', q), ...
                'FontSize',10, 'HorizontalAlignment','left', ...
                'BackgroundColor',[1 1 1 0.6]);
        end
    end

    % DEM points colored by sigma^2_n
    if ~isempty(phi_d)
        scatter(axB, Lhat_d, That_d, 150, sig_d, 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8, ...
            'DisplayName','DEM runs');
        cb2 = colorbar(axB);
        cb2.Label.String = '\sigma^2_n / \langle n\rangle^2 (density fluct.)';
        colormap(axB, hot);
    end

    % critical line lambda = L/2
    xline(axB, 1, '--k', 'LineWidth', 2, ...
        'Label','\lambda = L/2','LabelOrientation','horizontal', ...
        'HandleVisibility','off');

    set(axB, 'XScale','log','YScale','log');
    xlabel(axB,'\Lambda^* = \lambda / (L/2)','FontSize',14);
    ylabel(axB,'\tau^*    = \tau_H / \tau_{mfp}','FontSize',14);
    title (axB,'(B)  Dimensionless phase plane (\Lambda^* , \tau^*)','FontSize',14);
    axB.FontSize = 12;
    grid(axB,'on');
    legend(axB,'Location','southwest');

    % limits
    if ~isempty(Lhat_d)
        xlim(axB, [0.5*min([Lhat_d Lhat_th]) , 2*max([Lhat_d Lhat_th])]);
    else
        xlim(axB, [min(Lhat_th)*0.5, max(Lhat_th)*2]);
    end
    ylim(axB, yl);

    %% save
    if ~exist(opt.savepath, 'dir'), mkdir(opt.savepath); end
    base = fullfile(opt.savepath, 'Fig2_phase_diagram');
    saveas(f, [base '.png']);
    saveas(f, [base '.fig']);
    fprintf('\n[Fig2]  saved : %s.png / .fig\n', base);
end
