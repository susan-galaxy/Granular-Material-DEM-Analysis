%% =====================================================================
%  Microgravity Granular Gas: Characteristic Length & Time Scales,
%  Phase Diagram, and Predictive Clustering Map  (REVISION 3)
%  ---------------------------------------------------------------------
%  Changes vs v2:
%   • Empirical "Model A" effective-dissipation fit added: γ_eff captures
%     the residual dissipation that pure smooth-sphere theory misses.
%     One-parameter fit, plotted as golden curve in Fig 1.
%   • Fixed row/column vector mismatch bug in Model A fit (Bug-fix).
%   • Fig 1 legend now correctly includes Model A entry.
%  ---------------------------------------------------------------------
%  Theory:
%     λ_Noir = d / [ 6 φ χ(φ) √(1-e_n²) ]              (smooth sphere)
%     λ_full = d / [ 6 φ χ(φ) √ξ_full ]                (with rotation)
%     λ_modA = d / [ 6 φ χ(φ) √ξ_eff  ],
%             ξ_eff = (1-e_n²)·(1 + γ_eff), γ_eff fit from DEM
% =====================================================================
clear all; close all; clc;

%% ========== 1. System parameters =====================================
d       = 0.08;          % particle diameter [cm]
L       = 4.0;           % box side length [cm]
half_L  = L/2;
rho_s   = 1.08;
m_part  = rho_s*(pi/6)*d^3;

e_n     = 0.94;
beta_t  = 0.10;
mu_fric = 0.68;
q2      = 2/5;
kappa_T = 0.16;

xi_smooth = 1 - e_n^2;
xi_rot    = (8/7)*(1-beta_t^2)*(1 - kappa_T/q2);
xi_full   = xi_smooth + xi_rot;

fprintf('========= Dissipation parameters =========\n');
fprintf('  ξ_smooth (Noirhomme)        = %.4f\n', xi_smooth);
fprintf('  ξ_rot                       = %.4f\n', xi_rot);
fprintf('  ξ_full = ξ_smooth + ξ_rot   = %.4f\n', xi_full);
fprintf('==========================================\n\n');

%% ========== 2. Data path & file discovery ============================
pos_vel_dir = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Pos_Vel_data';
pv_pattern  = 'Data_*.mat';
tag_name    = 'Tag_2';

n_shells       = 20;
shell_edges    = linspace(0, half_L, n_shells+1);
shell_centers  = (shell_edges(1:end-1) + shell_edges(2:end))/2;
dist_from_wall = half_L - shell_centers;

delta_step = 2;
grid_bins  = 10;
grid_edges = linspace(-half_L, half_L, grid_bins+1);

% Fit configuration
saturation_handling = 'lower_bound';

lambda_upper_bound   = 500;
lambda_upper_extended= 5000;
lambda_max_phys      = 20;
sigma_thresh         = 2.0;
Tratio_thresh        = 0.3;
outer_T_frac         = 0.10;

R2_thresh            = 0.70;
Kn_thresh            = 1.0;

switch saturation_handling
    case 'lower_bound'
        lam_fit_upper = lambda_upper_bound;
        use_lower_bound_markers = true;
    case 'extend_bound'
        lam_fit_upper = lambda_upper_extended;
        use_lower_bound_markers = false;
    otherwise
        error('saturation_handling must be ''lower_bound'' or ''extend_bound''');
end
fprintf('Saturation handling: %s (λ fit upper bound = %g)\n\n', ...
        saturation_handling, lam_fit_upper);

pv_files = dir(fullfile(pos_vel_dir, pv_pattern));
if isempty(pv_files)
    error('No files matching "%s" in %s', pv_pattern, pos_vel_dir);
end
extract_phi = @(fn) str2double(regexp(fn,'(?<=phi)\d+\.?\d*','match','once'));
pv_map = containers.Map('KeyType','double','ValueType','char');
for k = 1:length(pv_files)
    p = extract_phi(pv_files(k).name);
    if ~isnan(p), pv_map(p) = pv_files(k).name; end
end
phi_data = sort(cell2mat(pv_map.keys));
phi_data = phi_data(:);   % force column to avoid row/col mismatch later
n_data   = length(phi_data);
fprintf('Found %d volume-fraction datasets:\n  φ = %s\n\n', ...
        n_data, mat2str(phi_data',4));

%% ========== 3. Per-φ DEM analysis ====================================
R                  = struct();
R.phi              = phi_data;
R.T_global         = zeros(n_data,1);
R.T_wall_fit       = zeros(n_data,1);
R.T_center_fit     = zeros(n_data,1);
R.lambda_DEM       = zeros(n_data,1);
R.lambda_DEM_CI    = zeros(n_data,2);
R.lambda_is_lb     = false(n_data,1);
R.fit_mode         = zeros(n_data,1);
R.R2               = nan(n_data,1);
R.lambda_theory    = zeros(n_data,1);
R.Kn_th            = zeros(n_data,1);
R.trust_fit        = false(n_data,1);
R.tau_lambda       = zeros(n_data,1);
R.tau_E            = zeros(n_data,1);
R.N_coll           = zeros(n_data,1);
R.sigma2_n         = zeros(n_data,1);
R.T_shell          = zeros(n_data, n_shells);
R.Tc_over_Tw       = zeros(n_data,1);

chi_f = @(phi) (1 - phi/2)./(1 - phi).^3;

for i = 1:n_data
    phi_val = phi_data(i);
    fprintf('[%2d/%2d] processing φ = %.4f ...', i, n_data, phi_val);
    
    ld = load(fullfile(pos_vel_dir, pv_map(phi_val)));
    PD = ld.ExtractedData.(tag_name);
    valid_frames = 1:delta_step:length(PD.Pos);
    
    sh_v2 = zeros(n_shells,1); sh_cnt = zeros(n_shells,1);
    g_v2 = 0; g_cnt = 0;
    fluc_list = [];
    
    for jf = 1:length(valid_frames)
        j = valid_frames(jf);
        p = double(PD.Pos{j});  v = double(PD.Vel{j});
        if isempty(p), continue; end
        bad = any(isnan(p),2) | any(isnan(v),2);
        p(bad,:) = []; v(bad,:) = [];
        if isempty(p), continue; end
        
        vf = v - mean(v,1);
        g_v2  = g_v2  + sum(sum(vf.^2));
        g_cnt = g_cnt + size(p,1);
        
        chebd = max(abs(p),[],2);
        s_idx = discretize(chebd, shell_edges);
        for s = 1:n_shells
            mk = (s_idx == s);
            if any(mk)
                sh_v2(s)  = sh_v2(s)  + sum(sum(vf(mk,:).^2));
                sh_cnt(s) = sh_cnt(s) + sum(mk);
            end
        end
        
        ix = discretize(p(:,1), grid_edges);
        iy = discretize(p(:,2), grid_edges);
        iz = discretize(p(:,3), grid_edges);
        good = ~isnan(ix) & ~isnan(iy) & ~isnan(iz);
        if any(good)
            lin = sub2ind([grid_bins grid_bins grid_bins], ...
                          ix(good), iy(good), iz(good));
            counts = histcounts(lin, 0.5:1:(grid_bins^3+0.5));
            mn = mean(counts);
            if mn > 0
                fluc_list(end+1) = var(counts)/mn; %#ok<SAGROW>
            end
        end
    end
    
    R.T_global(i) = g_v2 / (3*g_cnt);
    T_prof = zeros(n_shells,1);
    for s = 1:n_shells
        if sh_cnt(s) > 0, T_prof(s) = sh_v2(s)/(3*sh_cnt(s)); end
    end
    R.T_shell(i,:) = T_prof';
    
    if ~isempty(fluc_list)
        R.sigma2_n(i) = mean(fluc_list);
    end
    
    valid_s = T_prof > 0;
    idxs    = find(valid_s);
    if length(idxs) < 4
        fprintf(' [SKIP: insufficient shells]\n');
        continue;
    end
    
    T_wall_outer = T_prof(idxs(end));
    T_cent_inner = T_prof(idxs(1));
    R.Tc_over_Tw(i) = T_cent_inner / T_wall_outer;
    
    % ===== STAGE 1: full fit =====
    z_full = dist_from_wall(idxs)';
    T_full = T_prof(idxs);
    [lam1, Tw1, ci1, ok1, R2_1] = do_exp_fit(z_full, T_full, ...
                                       lam_fit_upper, 50*max(T_prof));
    
    % ===== STAGE 2: refit if clustering indicators trigger =====
    is_clustered = (R.sigma2_n(i) > sigma_thresh) || ...
                   (R.Tc_over_Tw(i) < Tratio_thresh);
    
    used_outer = false;
    if is_clustered && ok1
        T_thresh = outer_T_frac * T_wall_outer;
        keep_mask = T_prof(idxs) > T_thresh;
        idxs_o = idxs(keep_mask);
        if length(idxs_o) >= 4
            [lam2, Tw2, ci2, ok2, R2_2] = do_exp_fit( ...
                dist_from_wall(idxs_o)', T_prof(idxs_o), ...
                lam_fit_upper, 50*max(T_prof));
            if ok2 && lam2 > 0 && isfinite(lam2)
                R.lambda_DEM(i)      = lam2;
                R.lambda_DEM_CI(i,:) = ci2;
                R.T_wall_fit(i)      = Tw2;
                R.fit_mode(i)        = 2;
                R.R2(i)              = R2_2;
                used_outer = true;
            end
        end
    end
    
    if ~used_outer
        if ok1
            R.lambda_DEM(i)      = lam1;
            R.lambda_DEM_CI(i,:) = ci1;
            R.T_wall_fit(i)      = Tw1;
            R.fit_mode(i)        = 1;
            R.R2(i)              = R2_1;
        else
            pf = polyfit(z_full, log(T_full), 1);
            R.lambda_DEM(i)      = -1/pf(1);
            R.T_wall_fit(i)      =  exp(pf(2));
            R.lambda_DEM_CI(i,:) = [R.lambda_DEM(i) R.lambda_DEM(i)];
            R.fit_mode(i)        = 3;
            y_pred = R.T_wall_fit(i) * exp(-z_full / R.lambda_DEM(i));
            SS_res = sum((T_full - y_pred).^2);
            SS_tot = sum((T_full - mean(T_full)).^2);
            if SS_tot > 0
                R.R2(i) = 1 - SS_res/SS_tot;
            end
        end
    end
    R.T_center_fit(i) = R.T_wall_fit(i) * exp(-half_L/R.lambda_DEM(i));
    
    R.lambda_theory(i) = d / (6 * phi_val * chi_f(phi_val) * sqrt(xi_smooth));
    R.Kn_th(i)         = R.lambda_theory(i) / L;
    R.trust_fit(i) = (R.R2(i) > R2_thresh) && (R.Kn_th(i) < Kn_thresh);
    
    if use_lower_bound_markers
        R.lambda_is_lb(i) = R.lambda_DEM(i) > lambda_max_phys;
    else
        R.lambda_is_lb(i) = false;
    end
    
    if R.lambda_DEM(i) > 0 && R.T_wall_fit(i) > 0
        v_T              = sqrt(R.T_wall_fit(i));
        R.tau_lambda(i)  = R.lambda_DEM(i) / v_T;
        chi_i            = chi_f(phi_val);
        R.tau_E(i)       = sqrt(pi)*d / (24 * phi_val * chi_i * v_T);
        R.N_coll(i)      = R.tau_lambda(i) / R.tau_E(i);
    end
    
    mode_tag = {'full','outer','poly'}; 
    sat_tag = ''; if R.lambda_is_lb(i), sat_tag = ' [SAT]'; end
    trust_tag = ''; if ~R.trust_fit(i), trust_tag = '  [UNTRUSTED]'; end
    fprintf(' λ=%6.2f%s  R²=%.2f  Kn=%5.2f  T_c/T_w=%.3f  σ²=%5.2f  mode:%s%s\n', ...
        R.lambda_DEM(i), sat_tag, R.R2(i), R.Kn_th(i), ...
        R.Tc_over_Tw(i), R.sigma2_n(i), ...
        mode_tag{max(R.fit_mode(i),1)}, trust_tag);
end

%% ========== 4. Theory on dense φ-grid ================================
phi_th  = logspace(log10(0.001), log10(0.5), 600)';
chi_th  = chi_f(phi_th);

lambda_Noir = d ./ (6*phi_th.*chi_th*sqrt(xi_smooth));
lambda_full = d ./ (6*phi_th.*chi_th*sqrt(xi_full));

ic1 = find(lambda_Noir <= half_L, 1, 'first');
ic2 = find(lambda_full <= half_L, 1, 'first');
phi_c_Noir = phi_th(ic1);
phi_c_full = phi_th(ic2);
fprintf('\nCritical phi: Noirhomme=%.4f , full=%.4f\n', phi_c_Noir, phi_c_full);

lambda_Noir_at_DEM = d ./ (6*phi_data .* chi_f(phi_data) * sqrt(xi_smooth));

T_wall_ref = mean(R.T_wall_fit(R.T_wall_fit > 0 & ~R.lambda_is_lb));
v_T_ref    = sqrt(T_wall_ref);
tau_lambda_th = lambda_Noir / v_T_ref;
tau_E_th      = sqrt(pi)*d ./ (24 * phi_th .* chi_th * v_T_ref);

%% ========== 4b. Model A: effective dissipation fit ====================
fit_mask = R.trust_fit & R.lambda_DEM > 0 & ~R.lambda_is_lb;
if sum(fit_mask) >= 3
    % Force column vectors to avoid implicit row/col expansion
    phi_fit = phi_data(fit_mask); phi_fit = phi_fit(:);
    lam_fit = R.lambda_DEM(fit_mask); lam_fit = lam_fit(:);
    chi_fit = chi_f(phi_fit); chi_fit = chi_fit(:);
    
    lam_Noir_fit = d ./ (6 * phi_fit .* chi_fit * sqrt(xi_smooth));
    log_ratio = log(lam_Noir_fit ./ lam_fit);
    gamma_eff = exp(2 * mean(log_ratio)) - 1;
    
    % Bootstrap 95% CI on γ_eff
    nboot = 1000;
    boot_g = zeros(nboot,1);
    Nf = length(lam_fit);
    for b = 1:nboot
        idx_b = randi(Nf, Nf, 1);
        boot_g(b) = exp(2 * mean(log(lam_Noir_fit(idx_b) ./ lam_fit(idx_b)))) - 1;
    end
    g_ci = prctile(boot_g, [2.5 97.5]);
    
    xi_eff = (1 - e_n^2) * (1 + gamma_eff);
    lambda_modA = d ./ (6 * phi_th .* chi_th * sqrt(xi_eff));
    
    lam_modA_at_fit = d ./ (6 * phi_fit .* chi_fit * sqrt(xi_eff));
    rms_log = sqrt(mean((log(lam_fit) - log(lam_modA_at_fit)).^2));
    
    fprintf('\n=========== Model A (effective dissipation) ===========\n');
    fprintf('  γ_eff     = %.3f  (95%% CI: [%.3f, %.3f])\n', ...
            gamma_eff, g_ci(1), g_ci(2));
    fprintf('  ξ_eff     = %.3f  (smooth: %.3f, full: %.3f)\n', ...
            xi_eff, xi_smooth, xi_full);
    fprintf('  ξ_eff/ξ_smooth = %.2f  (interpretation: %.0f%% of γ_max activated)\n', ...
            xi_eff/xi_smooth, 100*gamma_eff/(xi_full/xi_smooth - 1));
    fprintf('  RMS log-residual = %.3f\n', rms_log);
    fprintf('=======================================================\n\n');
else
    warning('Insufficient trusted points to fit γ_eff. Using γ_eff=0.');
    gamma_eff = 0;
    g_ci = [0 0];
    xi_eff = xi_smooth;
    lambda_modA = lambda_Noir;
    rms_log = NaN;
end

%% ========== 5. PLOTTING ==============================================
col_th     = [0.85 0.10 0.10];
col_th2    = [0.20 0.20 0.85];
col_modA   = [0.85 0.55 0.10];
col_dem    = [0.10 0.55 0.20];
col_dem2   = [0.85 0.35 0.10];
col_grid   = [0.45 0.45 0.45];
col_untrust= [0.65 0.65 0.65];
fnt = 'Times New Roman';

valid_m_full_tr  = R.fit_mode == 1 & ~R.lambda_is_lb &  R.trust_fit;
valid_m_full_un  = R.fit_mode == 1 & ~R.lambda_is_lb & ~R.trust_fit;
valid_m_outer_tr = R.fit_mode == 2 &  R.trust_fit;
valid_m_outer_un = R.fit_mode == 2 & ~R.trust_fit;
valid_m_lb       = R.lambda_is_lb;
valid_t          = R.tau_lambda > 0 & isfinite(R.tau_lambda) & ...
                   ~R.lambda_is_lb & R.trust_fit;
valid_s_arr      = R.sigma2_n > 0;

phi_Kn1 = phi_th(find(lambda_Noir <= L*Kn_thresh, 1, 'first'));
if isempty(phi_Kn1), phi_Kn1 = NaN; end
fprintf('Fit-validity boundary: Kn=%.1f at φ ≈ %.4f\n\n', Kn_thresh, phi_Kn1);

%% ---------- Figure 1: λ(φ) -------------------------------------------
fig1 = figure('Name','Fig1_Lambda_vs_phi','Position',[60 60 760 580],'Color','w');

% Three theory curves
loglog(phi_th, lambda_Noir, '-' , 'Color', col_th  , 'LineWidth', 3.2); hold on;
loglog(phi_th, lambda_full, '--', 'Color', col_th2 , 'LineWidth', 2.0);
loglog(phi_th, lambda_modA, '-' , 'Color', col_modA, 'LineWidth', 2.8);

% Trusted DEM points (filled)
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_full_tr, ...
            'o', col_dem);
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_outer_tr, ...
            's', col_dem2);

% Untrusted DEM points (hollow gray)
if any(valid_m_full_un)
    plot(phi_data(valid_m_full_un), R.lambda_DEM(valid_m_full_un), 'o', ...
        'MarkerSize', 10, 'LineWidth', 1.4, ...
        'MarkerEdgeColor', col_untrust, 'MarkerFaceColor', 'none');
end
if any(valid_m_outer_un)
    plot(phi_data(valid_m_outer_un), R.lambda_DEM(valid_m_outer_un), 's', ...
        'MarkerSize', 10, 'LineWidth', 1.4, ...
        'MarkerEdgeColor', col_untrust, 'MarkerFaceColor', 'none');
end

% Saturated lower-bound points
if any(valid_m_lb)
    phi_lb = phi_data(valid_m_lb);
    lam_lb = lambda_max_phys * ones(size(phi_lb));
    plot(phi_lb, lam_lb, '^', 'MarkerSize', 12, 'LineWidth', 1.6, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.85 0.85 0.85]);
    for k = 1:length(phi_lb)
        text(phi_lb(k), lam_lb(k)*1.5, '\uparrow', ...
            'HorizontalAlignment','center', 'FontSize', 18, 'Color','k', ...
            'FontWeight', 'bold');
    end
end

yline(half_L, ':', 'Color', col_grid, 'LineWidth', 2.2);
text(0.18, half_L*1.45, '$\lambda = L/2$ (cluster onset)', ...
    'Interpreter','latex','FontSize',15,'Color',col_grid);

xline(phi_c_Noir, '-.', 'Color', col_th, 'LineWidth', 1.6);
text(phi_c_Noir*1.12, 0.025, sprintf('$\\phi_c \\approx %.3f$', phi_c_Noir), ...
    'Interpreter','latex','FontSize',15,'Color',col_th);

if isfinite(phi_Kn1)
    xline(phi_Kn1, ':', 'Color', col_untrust, 'LineWidth', 2.0, ...
        'HandleVisibility','off');
    text(phi_Kn1*0.55, 200, ...
        sprintf('$\\leftarrow$ Knudsen regime (Kn$>%g$)', Kn_thresh), ...
        'Interpreter','latex','FontSize',13,'Color',col_untrust);
end

if use_lower_bound_markers
    yl = [1e-2 1e2];
else
    yl = [1e-2 1e4];
end
patch([phi_c_Noir 0.5 0.5 phi_c_Noir], [yl(1) yl(1) half_L half_L], ...
      [1 0.92 0.92], 'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Energy penetration depth $\lambda$ [cm]', 'Interpreter','latex','FontSize',22);

% Build legend (single definition, includes Model A)
leg_entries = { ...
    'Theory: $\xi_{\rm smooth}=1-e_n^2$ (Noirhomme)', ...
    sprintf('Theory: $\\xi_{\\rm full}=%.2f$ (with rotation)', xi_full), ...
    sprintf('Model A: $\\xi_{\\rm eff}=%.3f$, $\\gamma_{\\rm eff}=%.2f\\pm%.2f$', ...
            xi_eff, gamma_eff, max(gamma_eff-g_ci(1), g_ci(2)-gamma_eff)), ...
    'DEM (trusted, full fit)', ...
    'DEM (trusted, outer-only refit)'};
if any(valid_m_full_un) || any(valid_m_outer_un)
    leg_entries{end+1} = sprintf('DEM (untrusted: R$^2<%.2f$ or Kn$>%g$)', ...
                                  R2_thresh, Kn_thresh);
end
if use_lower_bound_markers && any(valid_m_lb)
    leg_entries{end+1} = sprintf('DEM lower bound ($\\lambda > %d$ cm)', lambda_max_phys);
end
legend(leg_entries, 'Interpreter','latex','FontSize',12,'Location','southwest');
legend boxoff;

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on');
xlim([min(phi_th) max(phi_th)]); ylim(yl); box on;

%% ---------- Figure 2: τ(φ) -------------------------------------------
fig2 = figure('Name','Fig2_Tau_vs_phi','Position',[110 110 760 580],'Color','w');

loglog(phi_th, tau_lambda_th, '-' , 'Color', col_th , 'LineWidth', 3.2); hold on;
loglog(phi_th, tau_E_th     , '--', 'Color', col_th2, 'LineWidth', 2.4);

if any(valid_t)
    loglog(phi_data(valid_t), R.tau_lambda(valid_t), 'o', ...
        'MarkerSize', 12, 'LineWidth', 1.6, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem);
    loglog(phi_data(valid_t), R.tau_E(valid_t), 's', ...
        'MarkerSize', 11, 'LineWidth', 1.6, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', [0.95 0.78 0.30]);
end

tau_c = half_L / v_T_ref;
yline(tau_c, ':', 'Color', col_grid, 'LineWidth', 2.2);
text(0.20, tau_c*1.6, '$\tau_c = (L/2)/v_T$', ...
    'Interpreter','latex','FontSize',15,'Color',col_grid);
xline(phi_c_Noir, '-.', 'Color', col_th, 'LineWidth', 1.6);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Characteristic time [s]', 'Interpreter','latex','FontSize',22);
legend({ ...
    '$\tau_\lambda = \lambda/v_T$ (theory)', ...
    '$\tau_E = \sqrt{\pi}d/(24\phi\chi v_T)$ (theory)', ...
    '$\tau_\lambda$ from DEM', ...
    '$\tau_E$ from DEM'}, ...
    'Interpreter','latex','FontSize',13,'Location','southwest');
legend boxoff;
set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on');
xlim([min(phi_th) max(phi_th)]); box on;

%% ---------- Figure 3: Coincident theory breakdown ===================
fig3 = figure('Name','Fig3_Breakdown','Position',[160 60 820 600],'Color','w');

lam_ratio = R.lambda_DEM ./ lambda_Noir_at_DEM;
ok_ratio  = R.lambda_DEM > 0 & ~R.lambda_is_lb;

yyaxis left
plot(phi_data(ok_ratio), lam_ratio(ok_ratio), 'o-', ...
    'LineWidth', 2.5, 'MarkerSize', 11, ...
    'Color', col_dem, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_dem, 'LineWidth', 1.8, 'HandleVisibility','off');
ylabel('$\lambda_{\rm DEM}/\lambda_{\rm theory}$', ...
       'Interpreter','latex','FontSize',22,'Color',col_dem);
ax = gca; ax.YColor = col_dem;
ylim([0 1.8]);
text(min(phi_data)*1.5, 1.08, 'theory holds', ...
    'Interpreter','latex','FontSize',13,'Color',col_dem);

yyaxis right
semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 's-', ...
    'LineWidth', 2.5, 'MarkerSize', 11, ...
    'Color', col_th, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_th);
yline(1, '--', 'Color', col_th, 'LineWidth', 1.8, 'HandleVisibility','off');
ylabel('$\sigma_n^2/\langle n\rangle$', ...
       'Interpreter','latex','FontSize',22,'Color',col_th);
ax.YColor = col_th;
set(gca, 'YScale', 'log');

xline(phi_c_Noir, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8, ...
    'HandleVisibility','off');
text(phi_c_Noir*1.10, 30, sprintf('$\\phi_c \\approx %.3f$', phi_c_Noir), ...
    'Interpreter','latex','FontSize',15,'Color',[0.4 0.4 0.4]);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
xlim([min(phi_data)*0.7, max(phi_data)*1.4]);

legend({'$\lambda_{\rm DEM}/\lambda_{\rm theory}$', ...
        '$\sigma_n^2/\langle n\rangle$'}, ...
    'Interpreter','latex','FontSize',14,'Location','northwest');
legend boxoff;
box on;

%% ---------- Figure 4: Predictive phase diagram (φ, e_n) =============
fig4 = figure('Name','Fig4_Predictive_Phase_Map','Position',[210 110 840 620],'Color','w');

phi_g = logspace(log10(0.001), log10(0.3), 360);
en_g  = linspace(0.55, 0.999, 240);
[PHI, EN] = meshgrid(phi_g, en_g);
CHI   = chi_f(PHI);
XI_g  = 1 - EN.^2;
LAMR  = d ./ (6 * PHI .* CHI .* sqrt(XI_g)) / L;

contourf(PHI, EN, log10(LAMR), 24, 'LineColor','none'); hold on;
colormap(gca, parula);
cb = colorbar;
cb.Label.String       = '$\log_{10}(\lambda/L)$';
cb.Label.Interpreter  = 'latex';
cb.Label.FontSize     = 18;
cb.FontSize           = 14;
cb.FontName           = fnt;

[C_h, h_h] = contour(PHI, EN, LAMR, [0.5 0.5], 'r-', 'LineWidth', 3.5);
clabel(C_h, h_h, 'Color', 'r', 'FontSize', 13, 'LabelSpacing', 500, ...
       'FontWeight','bold');

[C2, h2] = contour(PHI, EN, LAMR, [0.1 1 5], 'w--', 'LineWidth', 1.3);
clabel(C2, h2, 'Color', 'w', 'FontSize', 11);

sigma_valid = R.sigma2_n(R.sigma2_n > 0);
if isempty(sigma_valid), clust_thr = inf;
else, clust_thr = max(2.0, 2*median(sigma_valid)); end
is_clust = R.sigma2_n > clust_thr;

plot(phi_data(~is_clust), e_n*ones(sum(~is_clust),1), 'o', ...
    'MarkerSize', 11, 'LineWidth', 2.2, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','w');
plot(phi_data(is_clust),  e_n*ones(sum(is_clust),1),  's', ...
    'MarkerSize', 13, 'LineWidth', 2.5, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor', col_th);

text(0.10, 0.70, '\textbf{Clustering}', 'Interpreter','latex', ...
    'FontSize', 19, 'Color', 'w');
text(0.0035, 0.96, '\textbf{Homogeneous Gas}', 'Interpreter','latex', ...
    'FontSize', 19, 'Color', 'w');
text(0.012, 0.91, sprintf('Our DEM: $e_n=%.2f$', e_n), ...
    'Interpreter','latex','FontSize',14,'Color','k', ...
    'BackgroundColor','w','EdgeColor','k','Margin',2);

set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Normal restitution $e_n$', 'Interpreter','latex','FontSize',22);
xlim([min(phi_g) max(phi_g)]); ylim([min(en_g) max(en_g)]);
box on;

datacursormode(fig4, 'off');
delete(findall(fig4, 'Tag', 'DataTipMarker'));
delete(findall(fig4, 'Type', 'datatip'));

legend({'$\lambda/L = 1/2$ (boundary)', '$\lambda/L = 0.1,\ 1,\ 5$', ...
        'DEM: gas', 'DEM: clustered'}, ...
    'Interpreter','latex','FontSize',12,'Location','southwest', ...
    'Color',[1 1 1 0.85]);

%% ---------- Figure 5: Order parameter ================================
fig5 = figure('Name','Fig5_Order_Param','Position',[260 60 720 540],'Color','w');

semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.6, 'MarkerSize', 12, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;

yline(1, '--', 'Color', col_grid, 'LineWidth', 2);
text(min(phi_data)*1.5, 1.3, 'Poisson (ideal gas)', ...
    'FontSize', 14, 'Color', col_grid);

xline(phi_c_Noir, '-.', 'Color', col_th, 'LineWidth', 2);
text(phi_c_Noir*1.10, max(R.sigma2_n)*0.7, ...
    sprintf('$\\phi_c \\approx %.3f$', phi_c_Noir), ...
    'Interpreter','latex','FontSize',16,'Color',col_th);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter','latex','FontSize',22);

set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
box on;

%% ---------- Figure 6: T(z) profiles + fits (sparse legend) ============
fig6 = figure('Name','Fig6_Tz_Profiles','Position',[310 110 800 580],'Color','w');

cmap6 = jet(n_data);
legend_step = max(1, round(n_data/8));

for i = 1:n_data
    T_prof = R.T_shell(i,:);
    vs = T_prof > 0;
    if ~any(vs), continue; end
    
    Tw_norm = T_prof(find(vs,1,'last'));
    show_legend = mod(i-1, legend_step) == 0 || i == n_data;
    
    plot(dist_from_wall(vs), T_prof(vs)/Tw_norm, 'o', ...
        'MarkerSize', 5, 'Color', cmap6(i,:), ...
        'MarkerFaceColor', cmap6(i,:), 'HandleVisibility','off'); hold on;
    
    if R.lambda_DEM(i) > 0
        zd = linspace(0, half_L, 200);
        if show_legend
            leg_str = sprintf('$\\phi=%.3f,\\ \\lambda=%.2f$', ...
                              phi_data(i), R.lambda_DEM(i));
            plot(zd, exp(-zd/R.lambda_DEM(i)), '-', ...
                'Color', cmap6(i,:), 'LineWidth', 1.4, ...
                'DisplayName', leg_str);
        else
            plot(zd, exp(-zd/R.lambda_DEM(i)), '-', ...
                'Color', cmap6(i,:), 'LineWidth', 1.4, ...
                'HandleVisibility','off');
        end
    end
end

xlabel('Distance from wall $z$ [cm]', 'Interpreter','latex','FontSize',22);
ylabel('$T(z) / T_{\rm wall}$', 'Interpreter','latex','FontSize',22);
yline(1, ':', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
legend('Interpreter','latex','FontSize',12,'Location','northeast');
legend boxoff;
set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
ylim([0 1.2]); xlim([0 half_L]); box on;

%% ---------- Figure 7: Four-witness panel =============================
fig7 = figure('Name','Fig7_Four_Witnesses','Position',[60 60 1100 800],'Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

x_lo = min(phi_data)*0.7;
x_hi = max(phi_data)*1.4;

% --- (a) λ(φ) ---
nexttile;
loglog(phi_th, lambda_Noir, '-', 'Color', col_th, 'LineWidth', 2.5); hold on;
loglog(phi_th, lambda_modA, '-', 'Color', col_modA, 'LineWidth', 2.2);
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_full_tr, 'o', col_dem);
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_outer_tr, 's', col_dem2);
if any(valid_m_full_un)
    plot(phi_data(valid_m_full_un), R.lambda_DEM(valid_m_full_un), 'o', ...
        'MarkerSize', 7, 'LineWidth', 1.2, ...
        'MarkerEdgeColor', col_untrust, 'MarkerFaceColor', 'none');
end
if any(valid_m_outer_un)
    plot(phi_data(valid_m_outer_un), R.lambda_DEM(valid_m_outer_un), 's', ...
        'MarkerSize', 7, 'LineWidth', 1.2, ...
        'MarkerEdgeColor', col_untrust, 'MarkerFaceColor', 'none');
end
if any(valid_m_lb)
    plot(phi_data(valid_m_lb), lambda_max_phys*ones(sum(valid_m_lb),1), ...
        '^', 'MarkerSize',9, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.85 .85 .85]);
end
yline(half_L, ':', 'Color', col_grid, 'LineWidth', 1.8);
xline(phi_c_Noir, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6);
if isfinite(phi_Kn1)
    xline(phi_Kn1, ':', 'Color', col_untrust, 'LineWidth', 1.5);
end
ylabel('$\lambda$ [cm]', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(a)', 'Units','normalized','FontSize',18,'FontWeight','bold');
text(0.02, 0.78, '$\lambda = L/2$', 'Units','normalized', ...
    'Interpreter','latex','FontSize',13,'Color',col_grid);
set(gca,'FontSize',15,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
if use_lower_bound_markers
    ylim_a = [1e-1 1e2];
else
    ylim_a = [1e-1 1e4];
end
xlim([x_lo x_hi]); ylim(ylim_a);
box on;

% --- (b) T_c/T_w(φ) ---
nexttile;
ok_tc = R.Tc_over_Tw > 0;
semilogx(phi_data(ok_tc), R.Tc_over_Tw(ok_tc), 'o-', ...
    'LineWidth', 2.2, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5);
xline(phi_c_Noir, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6);
ylabel('$T_{\rm center}/T_{\rm wall}$', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(b)', 'Units','normalized','FontSize',18,'FontWeight','bold');
text(0.02, 0.78, 'equilibrium', 'Units','normalized', ...
    'FontSize',12,'Color',col_grid);
set(gca,'FontSize',15,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
xlim([x_lo x_hi]); ylim([0 1.2]);
box on;

% --- (c) σ²_n/⟨n⟩(φ) ---
nexttile;
loglog(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.2, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5);
xline(phi_c_Noir, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6);
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\sigma_n^2/\langle n\rangle$', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(c)', 'Units','normalized','FontSize',18,'FontWeight','bold');
text(0.02, 0.78, 'Poisson', 'Units','normalized', ...
    'FontSize',12,'Color',col_grid);
set(gca,'FontSize',15,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
xlim([x_lo x_hi]);
box on;

% --- (d) τ_E(φ) ---
nexttile;
loglog(phi_th, tau_E_th, '--', 'Color', col_th2, 'LineWidth', 2.2); hold on;
if any(valid_t)
    loglog(phi_data(valid_t), R.tau_E(valid_t), 's', ...
        'MarkerSize', 9, 'LineWidth', 1.4, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', [0.95 0.78 0.30]);
end
xline(phi_c_Noir, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6);
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\tau_E$ [s]', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(d)', 'Units','normalized','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',15,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
xlim([x_lo x_hi]);
box on;

title(tlo, 'Four independent witnesses of the gas-to-clustering transition', ...
      'FontSize', 16, 'FontName', fnt);

%% ========== 6. Summary table & save ==================================
fprintf('\n=========== Summary ===========\n');
fprintf('  φ      T_global   T_wall    λ_DEM   λ_Noir  λ_modA  Kn   R²   fit    T_c/T_w  σ²/⟨n⟩  trust\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
mode_str = {'full ','outer','poly '};
% Interpolate Model A theory at each DEM φ for comparison
lambda_modA_at_DEM = d ./ (6 * phi_data .* chi_f(phi_data) * sqrt(xi_eff));
for i = 1:n_data
    lam_th_i = lambda_Noir_at_DEM(i);
    lam_mA_i = lambda_modA_at_DEM(i);
    sat = ' '; if R.lambda_is_lb(i), sat = '*'; end
    fm = max(R.fit_mode(i),1);
    trust_str = '  Y'; if ~R.trust_fit(i), trust_str = '  N'; end
    fprintf('%.4f  %.3e  %.3e  %6.2f%s  %5.2f   %5.2f  %5.2f  %.2f  %s  %6.3f   %6.2f%s\n', ...
        phi_data(i), R.T_global(i), R.T_wall_fit(i), R.lambda_DEM(i), sat, ...
        lam_th_i, lam_mA_i, R.Kn_th(i), R.R2(i), mode_str{fm}, ...
        R.Tc_over_Tw(i), R.sigma2_n(i), trust_str);
end
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('* = λ saturated (lower bound only)\n');
fprintf('Y/N: trusted (R² > %.2f AND Kn < %g) / untrusted\n\n', ...
        R2_thresh, Kn_thresh);

save_path = fullfile(pos_vel_dir, 'Phase_Diagram_Analysis_v3.mat');
save(save_path, 'R', 'phi_data', 'phi_th', 'lambda_Noir', 'lambda_full', ...
     'lambda_modA', 'lambda_Noir_at_DEM', 'lambda_modA_at_DEM', ...
     'tau_lambda_th', 'tau_E_th', ...
     'phi_c_Noir', 'phi_c_full', 'xi_smooth', 'xi_full', 'xi_eff', ...
     'gamma_eff', 'g_ci', 'rms_log', ...
     'lambda_max_phys', 'sigma_thresh', 'Tratio_thresh', ...
     'saturation_handling', 'lam_fit_upper', ...
     'R2_thresh', 'Kn_thresh', 'phi_Kn1', '-v7.3');
fprintf('Results saved to: %s\n', save_path);

fprintf('\nDone! 7 figures produced.\n');


%% ========== LOCAL FUNCTIONS (must be at end of script) ===============
function [lam, T_w, ci, ok, R2] = do_exp_fit(z_, T_, lam_upper, T_upper)
    try
        ft  = fittype('a*exp(-x/b)','independent','x');
        opt = fitoptions(ft);
        opt.StartPoint = [T_(end), 0.5];
        opt.Lower      = [0, 1e-3];
        opt.Upper      = [T_upper, lam_upper];
        [fr, ~, ~] = fit(z_, T_, ft, opt);
        ci_all = confint(fr, 0.95);
        lam = fr.b;
        T_w = fr.a;
        ci  = ci_all(:,2)';
        ok  = isfinite(lam) && lam > 0;
        y_pred = fr.a * exp(-z_ / fr.b);
        SS_res = sum((T_ - y_pred).^2);
        SS_tot = sum((T_ - mean(T_)).^2);
        if SS_tot > 0
            R2 = 1 - SS_res / SS_tot;
        else
            R2 = NaN;
        end
    catch
        lam = NaN; T_w = NaN; ci = [NaN NaN]; ok = false; R2 = NaN;
    end
end

function plot_helper(phi_v, lam_v, lci_v, mask, marker, faceColor)
    if ~any(mask), return; end
    phi_s = phi_v(mask);
    lam_s = lam_v(mask);
    lci_s = lci_v(mask,:);
    
    lneg = lam_s - lci_s(:,1);
    lpos = lci_s(:,2) - lam_s;
    bad_n = lneg <= 0 | ~isfinite(lneg) | lneg > 0.5*lam_s;
    bad_p = lpos <= 0 | ~isfinite(lpos) | lpos > 0.5*lam_s;
    lneg(bad_n) = 0;
    lpos(bad_p) = 0;
    
    errorbar(phi_s, lam_s, lneg, lpos, ...
        marker, 'MarkerSize', 11, 'LineWidth', 1.5, ...
        'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColor, ...
        'CapSize', 6, 'LineStyle', 'none');
end