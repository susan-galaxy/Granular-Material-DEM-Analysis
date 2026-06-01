function [C_fit, phi_c_theory] = Fig1_LengthVsPhi(results, opt)
% Figure 1 :  characteristic length lambda(phi) — theory vs DEM
%
% Inputs:
%   results : struct array from dem_load_all
%   opt     : struct with fields (all optional)
%       .L           - box size (cm)            default 4
%       .d           - particle diameter (cm)   default median of data
%       .en          - restitution coeff        default 0.94
%       .Cscale      - prefactor in theory_lambda  default 1
%       .savepath    - directory to save .png/.fig (default current)
%       .fit_C       - if true, fit Cscale so theory matches DEM in
%                      log–log sense (default true)

    if nargin < 2, opt = struct(); end
    if ~isfield(opt, 'L'),        opt.L = 4.0; end
    if ~isfield(opt, 'en'),       opt.en = 0.94; end
    if ~isfield(opt, 'Cscale'),   opt.Cscale = 1.0; end
    if ~isfield(opt, 'savepath'), opt.savepath = pwd; end
    if ~isfield(opt, 'fit_C'),    opt.fit_C = true; end
    if ~isfield(opt, 'd')
        if ~isempty(results)
            opt.d = median([results.d]);
        else
            opt.d = 0.08;
        end
    end

    %% Theory curve --------------------------------------------------------
    phi_th  = logspace(-3.5, log10(0.6), 400);
    lam_th0 = theory_lambda(phi_th, opt.en, opt.d, 'Cscale', 1.0);

    %% DEM points ----------------------------------------------------------
    phi_d   = [results.phi];
    lam_d   = [results.lam_dem];
    cluster = [results.cluster_flag];

    valid = ~isnan(phi_d) & ~isnan(lam_d) & isfinite(lam_d);
    phi_d   = phi_d(valid);
    lam_d   = lam_d(valid);
    cluster = cluster(valid);

    %% Optional prefactor fit ---------------------------------------------
    C_fit = opt.Cscale;
    if opt.fit_C && ~isempty(phi_d)
        lam_th_at = theory_lambda(phi_d, opt.en, opt.d, 'Cscale', 1.0);
        good = isfinite(lam_th_at) & lam_th_at > 0;
        if any(good)
            C_fit = exp( mean( log(lam_d(good)) - log(lam_th_at(good)) ) );
        end
    end
    lam_th = lam_th0 * C_fit;

    %% Critical phi (theory) : lambda = L/2 -------------------------------
    phi_c_theory = NaN;
    diff_curve = lam_th - opt.L/2;
    sgn = sign(diff_curve);
    iCross = find(sgn(1:end-1) .* sgn(2:end) < 0, 1, 'first');
    if ~isempty(iCross)
        phi_c_theory = interp1(diff_curve(iCross:iCross+1), ...
                               phi_th  (iCross:iCross+1), 0, 'linear');
    end

    %% Plot ---------------------------------------------------------------
    f = figure('Color','w','Position',[100 100 900 700]);
    ax = axes('Parent', f); hold(ax,'on'); box(ax,'on');

    % theory pre-factor in "absolute" form:  lambda = (Cabs * d) / sqrt((1-en^2) phi chi)
    C0 = sqrt( (sqrt(pi)/6) * pi / (6 * (12/sqrt(pi))) );
    C_abs = C_fit * C0;

    % Theory line
    plot(ax, phi_th, lam_th, '-', 'LineWidth', 2.8, ...
         'Color', [0.10 0.30 0.85], ...
         'DisplayName', sprintf('Theory  \\lambda = %.3f\\cdot d / [(1-e_n^2)\\phi\\chi(\\phi)]^{1/2}', C_abs) );

    % L/2 line
    yline(ax, opt.L/2, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, ...
          'Label','\lambda = L/2  (critical)', 'LabelHorizontalAlignment','left', ...
          'HandleVisibility','off');

    % DEM points -- split by cluster flag
    if ~isempty(phi_d)
        if any(~cluster)
            scatter(ax, phi_d(~cluster), lam_d(~cluster), 110, ...
                'MarkerEdgeColor', [0 0.45 0.74], ...
                'MarkerFaceColor', [0.65 0.85 1.0], ...
                'LineWidth', 1.2, ...
                'DisplayName', 'DEM (uniform gas)');
        end
        if any(cluster)
            scatter(ax, phi_d(cluster), lam_d(cluster), 130, ...
                'Marker','d', ...
                'MarkerEdgeColor', [0.85 0.10 0.10], ...
                'MarkerFaceColor', [1.00 0.75 0.75], ...
                'LineWidth', 1.2, ...
                'DisplayName', 'DEM (clustered)');
        end
    end

    % phi_c theory mark
    if ~isnan(phi_c_theory)
        xline(ax, phi_c_theory, ':', 'Color', [0.20 0.55 0.20], ...
              'LineWidth', 1.6, ...
              'Label', sprintf('\\phi_c^{th} = %.4f', phi_c_theory), ...
              'LabelVerticalAlignment','bottom', ...
              'HandleVisibility','off');
    end

    set(ax, 'XScale','log', 'YScale','log');
    grid(ax, 'on');
    xlabel(ax, '\phi  (volume fraction)', 'FontSize', 16);
    ylabel(ax, '\lambda  (cm) — energy penetration depth', 'FontSize', 16);
    title (ax, sprintf('Characteristic length  \\lambda(\\phi)   |   d=%.2fcm, e_n=%.2f, L=%.1fcm', ...
        opt.d, opt.en, opt.L), 'FontSize', 15);
    legend(ax, 'Location', 'southwest', 'FontSize', 12);
    ax.FontSize = 13;
    xlim(ax, [3e-4, 0.6]);
    ylim(ax, [opt.d/10, max(opt.L*4, opt.d*100)]);

    %% Save ----------------------------------------------------------------
    if ~exist(opt.savepath, 'dir'), mkdir(opt.savepath); end
    base = fullfile(opt.savepath, 'Fig1_lambda_vs_phi');
    saveas(f, [base '.png']);
    saveas(f, [base '.fig']);
    fprintf('\n[Fig1]  saved : %s.png / .fig\n', base);
    fprintf('[Fig1]  fitted prefactor C_scale = %.3f\n', C_fit);
    if ~isnan(phi_c_theory)
        fprintf('[Fig1]  theoretical phi_c (lambda = L/2) = %.4f\n', phi_c_theory);
    end
end
