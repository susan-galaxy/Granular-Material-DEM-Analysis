%% =====================================================================
%  Microgravity Granular Gas: TWO-LENGTH-SCALE Phase Diagram (REV 4)
%  ---------------------------------------------------------------------
%  KEY THEORETICAL CHANGE vs v3:
%
%      v3 criterion (single length):    λ_T(φ) = L/2
%      v4 criterion (two lengths):      λ_b(φ) + λ_T(φ) = L/2
%
%  Reason: the v3 scaling balances diffusion (∂²u/∂ξ²) against
%  dissipation (χ/u) and gives the thermal decay length λ_T. It drops
%  the source S(ξ) because S → 0 far from walls. But ballistic particles
%  carry energy a distance λ_b INTO the bulk before thermalizing, so the
%  total reach of injected energy is λ_b + λ_T, not λ_T alone.
%
%  Length-scale definitions (all in cm):
%      λ_b(φ)     = d / [6 φ χ(φ)]               (ballistic m.f.p., e-indep)
%      λ_T(φ,ξ)   = d / [6 φ χ(φ) √ξ]             (thermal decay)
%      λ_total    = λ_b + λ_T                     (energy penetration)
%
%  Three flavors of ξ (the dissipation parameter):
%      ξ_smooth = 1 - e_n²                        (smooth Noirhomme form)
%      ξ_full   = ξ_smooth + ξ_rot                (with rotational d.o.f.)
%      ξ_eff    = (1-e_n²)(1+γ_eff)               (Model A empirical fit)
%
%  Comparison strategy:
%      • DEM measures the bulk T(z) decay → this is λ_T (NOT λ_b)
%      • So DEM points are plotted against λ_T curves (as in v3).
%      • The phase boundary is a separate horizontal line at λ_total=L/2,
%        which sits BELOW (i.e., at higher φ_c than) the v3 line.
%
%  Sanity check (analytic):  setting λ_b(1 + 1/√ξ) = L/2 in the dilute
%  limit χ≈1 gives  φ_c = d(1 + 1/√ξ) / (3L).
%  For d=0.08, L=4:
%      ξ=ξ_smooth=0.116 → 1+1/√ξ=3.93 → φ_c ≈ 0.0262
%      ξ=ξ_full ≈ 0.79  → 1+1/√ξ=2.12 → φ_c ≈ 0.0141
%  v3 (old) gave φ_c ≈ 0.0195 (smooth) and 0.0075 (full). The new
%  numbers are HIGHER because the ballistic contribution buys extra
%  reach. The DEM σ²-onset at ≈0.02 sits between smooth/full predictions.
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
mu_fric = 0.5;
q2      = 2/5;
kappa_T = 0.16;

xi_smooth = 1 - e_n^2;
xi_rot    = (8/7)*(1-beta_t^2)*(1 - kappa_T/q2);
xi_full   = xi_smooth + xi_rot;

fprintf('========= Dissipation parameters =========\n');
fprintf('  ξ_smooth (Noirhomme)        = %.4f\n', xi_smooth);
fprintf('  ξ_rot                       = %.4f\n', xi_rot);
fprintf('  ξ_full = ξ_smooth + ξ_rot   = %.4f\n', xi_full);
fprintf('  λ_T/λ_b ratios:\n');
fprintf('    smooth: 1/√ξ_smooth = %.2f  (thermal length is %.1f× ballistic)\n', ...
        1/sqrt(xi_smooth), 1/sqrt(xi_smooth));
fprintf('    full  : 1/√ξ_full   = %.2f\n', 1/sqrt(xi_full));
fprintf('==========================================\n\n');

%% ========== 2. Data path & file discovery ============================
pos_vel_dir = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Pos_Vel_data/En0_94';
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

%% ========== 3. Per-φ DEM analysis (UNCHANGED from v3) ================
% This block measures the empirical thermal decay length λ_T_DEM from
% exponential fits to T(z). All comparisons remain valid: λ_DEM ↔ λ_T.
R                  = struct();
R.phi              = phi_data;
R.T_global         = zeros(n_data,1);
R.T_wall_fit       = zeros(n_data,1);
R.T_center_fit     = zeros(n_data,1);
R.lambda_DEM       = zeros(n_data,1);   % this is λ_T from DEM, NOT λ_b
R.lambda_DEM_CI    = zeros(n_data,2);
R.lambda_is_lb     = false(n_data,1);
R.fit_mode         = zeros(n_data,1);
R.R2               = nan(n_data,1);
R.lambda_T_theory  = zeros(n_data,1);   % renamed from lambda_theory
R.lambda_b_theory  = zeros(n_data,1);   % NEW
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
    
    % Theory values at this φ (using smooth ξ as the reference here)
    R.lambda_T_theory(i) = d / (6 * phi_val * chi_f(phi_val) * sqrt(xi_smooth));
    R.lambda_b_theory(i) = d / (6 * phi_val * chi_f(phi_val));
    R.Kn_th(i)           = R.lambda_T_theory(i) / L;
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
    fprintf(' λ_T=%6.2f%s  λ_b=%5.2f  R²=%.2f  Kn=%5.2f  T_c/T_w=%.3f  σ²=%5.2f  mode:%s%s\n', ...
        R.lambda_DEM(i), sat_tag, R.lambda_b_theory(i), R.R2(i), R.Kn_th(i), ...
        R.Tc_over_Tw(i), R.sigma2_n(i), ...
        mode_tag{max(R.fit_mode(i),1)}, trust_tag);
end

%% ========== 4. Theory on dense φ-grid ================================
phi_th  = logspace(log10(0.001), log10(0.5), 600)';
chi_th  = chi_f(phi_th);

% Ballistic mean free path (e-independent)
lambda_b   = d ./ (6*phi_th.*chi_th);

% Thermal decay length, three flavors of ξ
lambda_T_Noir = d ./ (6*phi_th.*chi_th*sqrt(xi_smooth));
lambda_T_full = d ./ (6*phi_th.*chi_th*sqrt(xi_full));

% Total energy penetration (NEW key quantity)
lambda_tot_Noir = lambda_b + lambda_T_Noir;
lambda_tot_full = lambda_b + lambda_T_full;

% === Old (v3) criterion: where λ_T = L/2 ============================
ic1_old = find(lambda_T_Noir <= half_L, 1, 'first');
ic2_old = find(lambda_T_full <= half_L, 1, 'first');
phi_c_Noir_old = phi_th(ic1_old);
phi_c_full_old = phi_th(ic2_old);

% === New (v4) criterion: where λ_b + λ_T = L/2 =======================
ic1 = find(lambda_tot_Noir <= half_L, 1, 'first');
ic2 = find(lambda_tot_full <= half_L, 1, 'first');
phi_c_Noir = phi_th(ic1);
phi_c_full = phi_th(ic2);

fprintf('\n=========== Critical volume fractions ===========\n');
fprintf('  Old criterion (λ_T = L/2):\n');
fprintf('     φ_c_smooth = %.4f      φ_c_full = %.4f\n', ...
        phi_c_Noir_old, phi_c_full_old);
fprintf('  New criterion (λ_b + λ_T = L/2):\n');
fprintf('     φ_c_smooth = %.4f      φ_c_full = %.4f\n', ...
        phi_c_Noir, phi_c_full);
fprintf('  DEM observation: σ²/⟨n⟩ jumps at φ ≈ 0.02\n');
fprintf('=================================================\n\n');

% Theory values at each DEM phi (for Model A residuals and plots)
lambda_b_at_DEM    = d ./ (6*phi_data .* chi_f(phi_data));
lambda_T_Noir_at_DEM = d ./ (6*phi_data .* chi_f(phi_data) * sqrt(xi_smooth));

T_wall_ref = mean(R.T_wall_fit(R.T_wall_fit > 0 & ~R.lambda_is_lb));
v_T_ref    = sqrt(T_wall_ref);
tau_lambda_th = lambda_T_Noir / v_T_ref;
tau_E_th      = sqrt(pi)*d ./ (24 * phi_th .* chi_th * v_T_ref);

%% ========== 4b. Model A: effective dissipation fit ====================
% Same fit as v3: fit ξ_eff so λ_T_modA matches λ_DEM.
fit_mask = R.trust_fit & R.lambda_DEM > 0 & ~R.lambda_is_lb;
if sum(fit_mask) >= 3
    phi_fit = phi_data(fit_mask); phi_fit = phi_fit(:);
    lam_fit = R.lambda_DEM(fit_mask); lam_fit = lam_fit(:);
    chi_fit = chi_f(phi_fit); chi_fit = chi_fit(:);
    
    lam_Noir_fit = d ./ (6 * phi_fit .* chi_fit * sqrt(xi_smooth));
    log_ratio = log(lam_Noir_fit ./ lam_fit);
    gamma_eff = exp(2 * mean(log_ratio)) - 1;
    
    nboot = 1000;
    boot_g = zeros(nboot,1);
    Nf = length(lam_fit);
    for b = 1:nboot
        idx_b = randi(Nf, Nf, 1);
        boot_g(b) = exp(2 * mean(log(lam_Noir_fit(idx_b) ./ lam_fit(idx_b)))) - 1;
    end
    g_ci = prctile(boot_g, [2.5 97.5]);
    
    xi_eff = (1 - e_n^2) * (1 + gamma_eff);
    lambda_T_modA   = d ./ (6 * phi_th .* chi_th * sqrt(xi_eff));
    lambda_tot_modA = lambda_b + lambda_T_modA;
    
    lam_modA_at_fit = d ./ (6 * phi_fit .* chi_fit * sqrt(xi_eff));
    rms_log = sqrt(mean((log(lam_fit) - log(lam_modA_at_fit)).^2));
    
    % Critical phi for Model A under NEW criterion
    ic3 = find(lambda_tot_modA <= half_L, 1, 'first');
    if isempty(ic3), phi_c_modA = NaN; else, phi_c_modA = phi_th(ic3); end
    
    % And under OLD criterion for comparison
    ic3_old = find(lambda_T_modA <= half_L, 1, 'first');
    if isempty(ic3_old), phi_c_modA_old = NaN; else, phi_c_modA_old = phi_th(ic3_old); end
    
    fprintf('=========== Model A (effective dissipation) ===========\n');
    fprintf('  γ_eff     = %.3f  (95%% CI: [%.3f, %.3f])\n', ...
            gamma_eff, g_ci(1), g_ci(2));
    fprintf('  ξ_eff     = %.3f  (smooth: %.3f, full: %.3f)\n', ...
            xi_eff, xi_smooth, xi_full);
    fprintf('  ξ_eff/ξ_smooth = %.2f  (%.0f%% of γ_max activated)\n', ...
            xi_eff/xi_smooth, 100*gamma_eff/(xi_full/xi_smooth - 1));
    fprintf('  RMS log-residual = %.3f\n', rms_log);
    fprintf('  φ_c (Model A, OLD criterion λ_T=L/2)       = %.4f\n', phi_c_modA_old);
    fprintf('  φ_c (Model A, NEW criterion λ_b+λ_T=L/2)  = %.4f\n', phi_c_modA);
    fprintf('=======================================================\n\n');
else
    warning('Insufficient trusted points to fit γ_eff. Using γ_eff=0.');
    gamma_eff = 0;
    g_ci = [0 0];
    xi_eff = xi_smooth;
    lambda_T_modA   = lambda_T_Noir;
    lambda_tot_modA = lambda_tot_Noir;
    phi_c_modA      = phi_c_Noir;
    phi_c_modA_old  = phi_c_Noir_old;
    rms_log = NaN;
end

lambda_T_modA_at_DEM = d ./ (6 * phi_data .* chi_f(phi_data) * sqrt(xi_eff));

%% ========== 5. PLOTTING ==============================================
col_th     = [0.85 0.10 0.10];   % smooth theory (red)
col_th2    = [0.20 0.20 0.85];   % full theory (blue)
col_modA   = [0.85 0.55 0.10];   % Model A (orange)
col_dem    = [0.10 0.55 0.20];   % DEM trusted (green)
col_dem2   = [0.85 0.35 0.10];   % DEM outer-refit (red-orange)
col_grid   = [0.45 0.45 0.45];
col_untrust= [0.65 0.65 0.65];
col_ball   = [0.40 0.20 0.55];   % λ_b (purple) — NEW
col_tot    = [0.00 0.00 0.00];   % λ_total (black bold) — NEW
fnt = 'Times New Roman';

valid_m_full_tr  = R.fit_mode == 1 & ~R.lambda_is_lb &  R.trust_fit;
valid_m_full_un  = R.fit_mode == 1 & ~R.lambda_is_lb & ~R.trust_fit;
valid_m_outer_tr = R.fit_mode == 2 &  R.trust_fit;
valid_m_outer_un = R.fit_mode == 2 & ~R.trust_fit;
valid_m_lb       = R.lambda_is_lb;
valid_t          = R.tau_lambda > 0 & isfinite(R.tau_lambda) & ...
                   ~R.lambda_is_lb & R.trust_fit;
valid_s_arr      = R.sigma2_n > 0;

phi_Kn1 = phi_th(find(lambda_T_Noir <= L*Kn_thresh, 1, 'first'));
if isempty(phi_Kn1), phi_Kn1 = NaN; end
fprintf('Fit-validity boundary: Kn=%.1f at φ ≈ %.4f\n\n', Kn_thresh, phi_Kn1);

%% ---------- Figure 1: λ_T, λ_b and λ_total vs φ ----------------------
fig1 = figure('Name','Fig1_Lambda_vs_phi','Position',[60 60 820 620],'Color','w');

% --- Thermal decay length λ_T (three flavors) — solid colored ---
loglog(phi_th, lambda_T_Noir, '-' , 'Color', col_th  , 'LineWidth', 2.5); hold on;
loglog(phi_th, lambda_T_full, '--', 'Color', col_th2 , 'LineWidth', 2.0);
loglog(phi_th, lambda_T_modA, '-' , 'Color', col_modA, 'LineWidth', 2.5);

% --- Ballistic mean free path λ_b — dot-dashed purple ---
loglog(phi_th, lambda_b, '-.', 'Color', col_ball, 'LineWidth', 2.4);

% --- Total energy penetration λ_b + λ_T_modA — bold black ---
loglog(phi_th, lambda_tot_modA, '-', 'Color', col_tot, 'LineWidth', 3.2);

% --- DEM points (these are empirical λ_T, compare to λ_T curves) ---
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_full_tr, ...
            'o', col_dem);
plot_helper(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_outer_tr, ...
            's', col_dem2);

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

% Phase boundary line (L/2) — now interpreted against λ_total
yline(half_L, ':', 'Color', col_grid, 'LineWidth', 2.2);
text(0.18, half_L*1.45, ...
    '$\lambda_b+\lambda_T = L/2$ (cluster onset, NEW)', ...
    'Interpreter','latex','FontSize',13,'Color',col_grid);

% Critical phi vertical line (NEW criterion, Model A)
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.8);
text(phi_c_modA*1.10, 0.025, ...
    sprintf('$\\phi_c \\approx %.3f$', phi_c_modA), ...
    'Interpreter','latex','FontSize',14,'Color',col_tot);

% Old criterion as faint dashed reference
xline(phi_c_Noir_old, ':', 'Color', col_th, 'LineWidth', 1.2);
text(phi_c_Noir_old*0.55, 0.012, ...
    sprintf('$\\phi_c^{\\rm old} = %.3f$', phi_c_Noir_old), ...
    'Interpreter','latex','FontSize',11,'Color',col_th);

if isfinite(phi_Kn1)
    xline(phi_Kn1, ':', 'Color', col_untrust, 'LineWidth', 2.0, ...
        'HandleVisibility','off');
    text(phi_Kn1*0.55, 200, ...
        sprintf('$\\leftarrow$ Knudsen regime (Kn$>%g$)', Kn_thresh), ...
        'Interpreter','latex','FontSize',12,'Color',col_untrust);
end

if use_lower_bound_markers
    yl = [1e-2 1e2];
else
    yl = [1e-2 1e4];
end
patch([phi_c_modA 0.5 0.5 phi_c_modA], [yl(1) yl(1) half_L half_L], ...
      [1 0.92 0.92], 'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Length scales [cm]', 'Interpreter','latex','FontSize',22);

leg_entries = { ...
    'Thermal $\lambda_T$: $\xi_{\rm smooth}$', ...
    sprintf('Thermal $\\lambda_T$: $\\xi_{\\rm full}=%.2f$', xi_full), ...
    sprintf('Thermal $\\lambda_T$: Model A ($\\xi_{\\rm eff}=%.2f$)', xi_eff), ...
    'Ballistic m.f.p. $\lambda_b = d/(6\phi\chi)$', ...
    '$\lambda_{\rm tot}=\lambda_b+\lambda_T$ (Model A)', ...
    'DEM $\lambda_T$ (trusted, full fit)', ...
    'DEM $\lambda_T$ (outer-only refit)'};
if any(valid_m_full_un) || any(valid_m_outer_un)
    leg_entries{end+1} = 'DEM (untrusted)';
end
if use_lower_bound_markers && any(valid_m_lb)
    leg_entries{end+1} = sprintf('DEM lower bound ($>%d$ cm)', lambda_max_phys);
end
legend(leg_entries, 'Interpreter','latex','FontSize',11,'Location','southwest');
legend boxoff;

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on');
xlim([min(phi_th) max(phi_th)]); ylim(yl); box on;

%% ---------- Figure 2: τ(φ) (UNCHANGED) -------------------------------
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
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Characteristic time [s]', 'Interpreter','latex','FontSize',22);
legend({ ...
    '$\tau_\lambda = \lambda_T/v_T$ (theory)', ...
    '$\tau_E = \sqrt{\pi}d/(24\phi\chi v_T)$ (theory)', ...
    '$\tau_\lambda$ from DEM', ...
    '$\tau_E$ from DEM'}, ...
    'Interpreter','latex','FontSize',13,'Location','southwest');
legend boxoff;
set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on');
xlim([min(phi_th) max(phi_th)]); box on;

%% ---------- Figure 3: λ_DEM/λ_T_theory + σ² (small label update) ====
fig3 = figure('Name','Fig3_Breakdown','Position',[160 60 820 600],'Color','w');

lam_ratio = R.lambda_DEM ./ lambda_T_Noir_at_DEM;
ok_ratio  = R.lambda_DEM > 0 & ~R.lambda_is_lb;

yyaxis left
plot(phi_data(ok_ratio), lam_ratio(ok_ratio), 'o-', ...
    'LineWidth', 2.5, 'MarkerSize', 11, ...
    'Color', col_dem, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_dem, 'LineWidth', 1.8, 'HandleVisibility','off');
ylabel('$\lambda_{T,\rm DEM}/\lambda_{T,\rm theory}^{\rm smooth}$', ...
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

xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.8, ...
    'HandleVisibility','off');
text(phi_c_modA*1.10, 30, ...
    sprintf('$\\phi_c \\approx %.3f$ (new)', phi_c_modA), ...
    'Interpreter','latex','FontSize',15,'Color',col_tot);

xlabel('$\phi$', 'Interpreter','latex','FontSize',22);
set(gca,'XScale','log','FontSize',26,'FontName',fnt,'LineWidth',2.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
xlim([min(phi_data)*0.7, max(phi_data)*1.4]);

legend({'$\lambda_{T,\rm DEM}/\lambda_{T,\rm theory}$', ...
        '$\sigma_n^2/\langle n\rangle$'}, ...
    'Interpreter','latex','FontSize',24,'Location','northwest');
legend boxoff;
box on;

%% ---------- Figure 4: Predictive phase diagram (φ, e_n) — NEW =========
% Coloring uses log10(λ_total / L) = log10[(λ_b + λ_T)/L].
% Phase boundary is the contour λ_total/L = 1/2.
fig4 = figure('Name','Fig4_Predictive_Phase_Map','Position',[210 110 860 640],'Color','w');

phi_g = logspace(log10(0.001), log10(0.3), 360);
en_g  = linspace(0.55, 0.999, 240);
[PHI, EN] = meshgrid(phi_g, en_g);
CHI   = chi_f(PHI);
XI_g  = 1 - EN.^2;
LB    = d ./ (6 * PHI .* CHI);                 % ballistic
LT    = d ./ (6 * PHI .* CHI .* sqrt(XI_g));   % thermal
LTOT  = LB + LT;                                % total (NEW)
LTOT_norm = LTOT / L;

contourf(PHI, EN, log10(LTOT_norm), 24, 'LineColor','none'); hold on;
colormap(gca, parula);
cb = colorbar;
cb.Label.String       = '$\log_{10}[(\lambda_b+\lambda_T)/L]$';
cb.Label.Interpreter  = 'latex';
cb.Label.FontSize     = 17;
cb.FontSize           = 14;
cb.FontName           = fnt;

% NEW phase boundary: λ_total/L = 1/2
[C_h, h_h] = contour(PHI, EN, LTOT_norm, [0.5 0.5], 'r-', 'LineWidth', 3.5);
clabel(C_h, h_h, 'Color', 'r', 'FontSize', 13, 'LabelSpacing', 500, ...
       'FontWeight','bold');

% OLD phase boundary (λ_T/L = 1/2) shown for comparison as dashed white
LT_norm = LT/L;
[Co, ho] = contour(PHI, EN, LT_norm, [0.5 0.5], 'w--', 'LineWidth', 2.0);
clabel(Co, ho, 'Color', 'w', 'FontSize', 11);

[C2, h2] = contour(PHI, EN, LTOT_norm, [0.1 1 5], 'w:', 'LineWidth', 1.0);
clabel(C2, h2, 'Color', 'w', 'FontSize', 10);

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

legend({'$(\lambda_b+\lambda_T)/L = 1/2$ (NEW boundary)', ...
        '$\lambda_T/L = 1/2$ (OLD boundary, v3)', ...
        '$(\lambda_b+\lambda_T)/L = 0.1,\ 1,\ 5$', ...
        'DEM: gas', 'DEM: clustered'}, ...
    'Interpreter','latex','FontSize',11,'Location','southwest', ...
    'Color',[1 1 1 0.85]);

%% ---------- Figure 5: σ² order parameter (UNCHANGED) =================
fig5 = figure('Name','Fig5_Order_Param','Position',[260 60 720 540],'Color','w');

semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.6, 'MarkerSize', 12, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;

yline(1, '--', 'Color', col_grid, 'LineWidth', 2);
text(min(phi_data)*1.5, 1.3, 'Poisson (ideal gas)', ...
    'FontSize', 14, 'Color', col_grid);

xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 2);
text(phi_c_modA*1.10, max(R.sigma2_n)*0.7, ...
    sprintf('$\\phi_c \\approx %.3f$ (new)', phi_c_modA), ...
    'Interpreter','latex','FontSize',15,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter','latex','FontSize',22);

set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
box on;

%% ---------- Figure 6: T(z) profiles + fits (UNCHANGED) ===============
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
            leg_str = sprintf('$\\phi=%.3f,\\ \\lambda_T=%.2f$', ...
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

%% ---------- Figure 7: Four-witness panel (minor updates) =============
fig7 = figure('Name','Fig7_Four_Witnesses','Position',[60 60 1100 800],'Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

x_lo = min(phi_data)*0.7;
x_hi = max(phi_data)*1.4;

% --- (a) λ(φ): now includes λ_b and λ_total ---
nexttile;
loglog(phi_th, lambda_T_Noir, '-', 'Color', col_th, 'LineWidth', 2.0); hold on;
loglog(phi_th, lambda_T_modA, '-', 'Color', col_modA, 'LineWidth', 2.0);
loglog(phi_th, lambda_b, '-.', 'Color', col_ball, 'LineWidth', 1.8);
loglog(phi_th, lambda_tot_modA, '-', 'Color', col_tot, 'LineWidth', 2.6);
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
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6);
if isfinite(phi_Kn1)
    xline(phi_Kn1, ':', 'Color', col_untrust, 'LineWidth', 1.5);
end
ylabel('Lengths [cm]', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(a)', 'Units','normalized','FontSize',18,'FontWeight','bold');
text(0.02, 0.78, '$\lambda_b+\lambda_T=L/2$', 'Units','normalized', ...
    'Interpreter','latex','FontSize',12,'Color',col_grid);
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
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6);
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
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6);
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
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6);
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\tau_E$ [s]', 'Interpreter','latex','FontSize',18);
text(0.02, 0.92, '(d)', 'Units','normalized','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',15,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
xlim([x_lo x_hi]);
box on;

title(tlo, 'Four independent witnesses of the gas-to-clustering transition', ...
      'FontSize', 16, 'FontName', fnt);

%% ---------- Figure 8 (REPLACEMENT): linear-scale criterion plot ======
% Why this replaces the old Fig 8:
%   The previous design plotted λ_b / (λ_b + λ_T) vs φ. But algebraically
%       λ_b / (λ_b + λ_T) = √ξ / (√ξ + 1)
%   which has NO φ dependence in the dilute-limit form λ ∝ 1/[φχ(φ)].
%   So the old Fig 8 was three flat horizontal lines — visually useless.
%
%   The physically interesting fact hiding in that flat ratio:
%   the ballistic contribution is a φ-INDEPENDENT fraction of the total
%   energy penetration (25-30% in our setup), i.e. the v3 single-length
%   approximation undercounts λ_tot by a uniform factor at all densities,
%   not just in the dilute regime.
%
%   The replacement plot below puts λ_b/L, λ_T/L, and (λ_b+λ_T)/L on a
%   linear y-axis. The phase criterion is the horizontal line y = 1/2.
%   The shift from φ_c_old (λ_T/L crossing 1/2) to φ_c_new (λ_tot/L
%   crossing 1/2) is visually obvious.
fig8 = figure('Name','Fig8_Linear_Criterion','Position',[360 60 820 580],'Color','w');

semilogx(phi_th, lambda_T_Noir/L , ':' , 'Color', col_th  , 'LineWidth', 1.8); hold on;
semilogx(phi_th, lambda_T_modA/L , '-' , 'Color', col_modA, 'LineWidth', 2.4);
semilogx(phi_th, lambda_b/L      , '-.', 'Color', col_ball, 'LineWidth', 2.4);
semilogx(phi_th, lambda_tot_modA/L,'-' , 'Color', col_tot , 'LineWidth', 3.2);

% Criterion line
yline(0.5, ':', 'Color', col_grid, 'LineWidth', 2.5);
text(0.0015, 0.55, 'criterion: $\lambda_{\rm tot}/L = 1/2$', ...
    'Interpreter','latex','FontSize',13,'Color',col_grid);

% Old (v3) and new (v4) critical phi as vertical guides
xline(phi_c_Noir_old, ':' , 'Color', col_th , 'LineWidth', 1.5);
xline(phi_c_modA    , '-.', 'Color', col_tot, 'LineWidth', 2.0);
text(phi_c_Noir_old*0.55, 1.05, ...
    sprintf('$\\phi_c^{\\rm v3}=%.3f$', phi_c_Noir_old), ...
    'Interpreter','latex','FontSize',12,'Color',col_th);
text(phi_c_modA*1.05, 1.05, ...
    sprintf('$\\phi_c^{\\rm v4}=%.3f$', phi_c_modA), ...
    'Interpreter','latex','FontSize',12,'Color',col_tot);

% Shade the clustered region (where λ_tot < L/2)
patch([phi_c_modA 0.1 0.1 phi_c_modA], [0 0 0.5 0.5], ...
      [1 0.92 0.92], 'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
text(0.045, 0.25, '\textbf{Clustering}', 'Interpreter','latex', ...
    'FontSize',16,'Color',[0.6 0.1 0.1]);
text(0.0025, 0.85, '\textbf{Homogeneous Gas}', 'Interpreter','latex', ...
    'FontSize',16,'Color',[0.1 0.5 0.1]);

% Show DEM σ²-onset for context (if data is available)
if any(R.sigma2_n > 0)
    [~, idx_onset] = min(abs(R.sigma2_n(R.sigma2_n > 0) - 2.0));
    phi_with_sigma = phi_data(R.sigma2_n > 0);
    phi_sigma_onset = phi_with_sigma(idx_onset);
    xline(phi_sigma_onset, '--', 'Color', col_dem, 'LineWidth', 1.6);
    text(phi_sigma_onset*1.05, 0.35, ...
        sprintf('$\\phi(\\sigma^2{=}2) \\approx %.3f$', phi_sigma_onset), ...
        'Interpreter','latex','FontSize',11,'Color',col_dem);
end

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Normalized length $\lambda / L$', 'Interpreter','latex','FontSize',22);

legend({'$\lambda_T/L$ (smooth, v3 form)', ...
        sprintf('$\\lambda_T/L$ (Model A, $\\xi_{\\rm eff}=%.2f$)', xi_eff), ...
        '$\lambda_b/L$ (ballistic m.f.p.)', ...
        '$(\lambda_b+\lambda_T)/L$ (v4, criterion-relevant)'}, ...
        'Interpreter','latex','FontSize',12,'Location','northeast');
legend boxoff;
set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on');
xlim([min(phi_th) 0.1]); ylim([0 1.2]);
box on;

%% ========== 6. Summary table & save ==================================
fprintf('\n=========== Summary ===========\n');
fprintf('  φ      T_global   λ_DEM   λ_T_th  λ_b_th  λ_tot   Kn   R²   fit    T_c/T_w  σ²    trust\n');
fprintf('-----------------------------------------------------------------------------------------\n');
mode_str = {'full ','outer','poly '};
for i = 1:n_data
    sat = ' '; if R.lambda_is_lb(i), sat = '*'; end
    fm = max(R.fit_mode(i),1);
    trust_str = '  Y'; if ~R.trust_fit(i), trust_str = '  N'; end
    lam_tot_i = R.lambda_b_theory(i) + R.lambda_T_theory(i);
    fprintf('%.4f  %.3e  %6.2f%s  %5.2f   %5.2f   %5.2f   %.2f  %.2f  %s  %6.3f   %6.2f%s\n', ...
        phi_data(i), R.T_global(i), R.lambda_DEM(i), sat, ...
        R.lambda_T_theory(i), R.lambda_b_theory(i), lam_tot_i, ...
        R.Kn_th(i), R.R2(i), mode_str{fm}, ...
        R.Tc_over_Tw(i), R.sigma2_n(i), trust_str);
end
fprintf('-----------------------------------------------------------------------------------------\n');
fprintf('* = λ saturated (lower bound only)\n');
fprintf('Y/N: trusted (R² > %.2f AND Kn < %g) / untrusted\n\n', ...
        R2_thresh, Kn_thresh);

fprintf('=========== Critical φ values (final) ===========\n');
fprintf('  Old criterion  λ_T = L/2:                NEW criterion  λ_b+λ_T = L/2:\n');
fprintf('    smooth:    %.4f                          smooth:    %.4f\n', ...
        phi_c_Noir_old, phi_c_Noir);
fprintf('    full:      %.4f                          full:      %.4f\n', ...
        phi_c_full_old, phi_c_full);
fprintf('    Model A:   %.4f                          Model A:   %.4f\n', ...
        phi_c_modA_old, phi_c_modA);
fprintf('  DEM σ²-onset:  ≈ 0.020\n');
fprintf('=================================================\n\n');

save_path = fullfile(pos_vel_dir, 'Phase_Diagram_Analysis_v4.mat');
save(save_path, 'R', 'phi_data', 'phi_th', ...
     'lambda_b', 'lambda_T_Noir', 'lambda_T_full', 'lambda_T_modA', ...
     'lambda_tot_Noir', 'lambda_tot_full', 'lambda_tot_modA', ...
     'lambda_T_Noir_at_DEM', 'lambda_b_at_DEM', 'lambda_T_modA_at_DEM', ...
     'tau_lambda_th', 'tau_E_th', ...
     'phi_c_Noir', 'phi_c_full', 'phi_c_modA', ...
     'phi_c_Noir_old', 'phi_c_full_old', 'phi_c_modA_old', ...
     'xi_smooth', 'xi_full', 'xi_eff', ...
     'gamma_eff', 'g_ci', 'rms_log', ...
     'lambda_max_phys', 'sigma_thresh', 'Tratio_thresh', ...
     'saturation_handling', 'lam_fit_upper', ...
     'R2_thresh', 'Kn_thresh', 'phi_Kn1', '-v7.3');
fprintf('Results saved to: %s\n', save_path);

fprintf('\nDone! 8 figures produced.\n');


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

%% 2plot
%% =====================================================================
%  PUBLICATION-GRADE REPLACEMENT FOR ALL 8 FIGURES OF Phase_v4.m
%  ---------------------------------------------------------------------
%  USAGE
%    Run Phase_v4.m first (so all variables are in base workspace).
%    Then run this file. The original figures will be closed and
%    publication-grade versions will be drawn in their place.
%
%  PER-FIGURE CHANGES
%
%  Fig 1 (λ_b, λ_T, λ_tot vs φ):
%    REMOVED  - ξ_smooth, ξ_full reference curves
%             - old (v3) φ_c reference line and label
%             - Knudsen-regime annotation
%             - lower-bound triangle markers
%             - untrusted gray DEM points
%    KEPT     - λ_b (purple, thin), λ_T Model A (orange), λ_tot (black bold)
%             - DEM trusted (green ○) and outer-refit (red-orange ▪)
%             - L/2 horizontal + φ_c vertical + clustered region shading
%
%  Fig 2 (characteristic times τ_λ, τ_E vs φ):
%    FIXED    - τ_λ recolored to col_modA (orange) to match the Model A
%               framework used everywhere else in v4
%             - τ_E recolored to a darker amber for marker contrast
%             - φ_c vertical line now gets a label
%             - explicit handle→legend mapping
%
%  Fig 3 (λ_DEM/λ_theory + σ² dual-axis):
%    FIXED    - denominator switched from λ_T_smooth to λ_T_Model_A
%               (consistent with the framework actually validated)
%             - cleaner dual-axis label coloring
%             - trusted-only DEM points
%
%  Fig 4 (predictive phase map in (φ, e_n)):
%    FIXED    - explicit handle→label legend mapping
%               (contourf no longer steals a legend slot)
%             - removed the 0.1/1/5 white-dotted contours
%             - removed the "Homogeneous Gas"/"Clustering" overlay text
%
%  Fig 5 (σ²/⟨n⟩ order parameter):
%    FIXED    - line and marker style match Fig 3 right axis
%             - φ_c label position guaranteed inside axes
%             - lighter clustered region shading
%
%  Fig 6 (T(z)/T_wall profiles):
%    FIXED    - jet (perceptually non-uniform) replaced with parula
%             - crowded legend replaced with a single φ colorbar
%             - profiles thinned to readable density (every k-th φ shown)
%
%  Fig 7 (four-witness 2×2 panel):
%    FIXED    - panel (a) reduced to v4-only (λ_b, λ_T_modA, λ_tot)
%               matching the cleaned-up Fig 1
%             - Knudsen vertical removed
%             - x-limits unified across all four panels
%             - panel labels (a)(b)(c)(d) repositioned to a fixed
%               normalized corner so they don't drift
%
%  Fig 8 (linear-scale criterion plot):
%    FIXED    - dropped the redundant λ_T smooth line
%             - Homogeneous/Clustering labels repositioned to avoid
%               overlap with the curves
%             - σ²=2 onset annotation only drawn if non-trivial
% =====================================================================

% ============== Sanity-check required variables ======================
required_vars = {'phi_th','lambda_b','lambda_T_Noir','lambda_T_modA', ...
                 'lambda_tot_modA','phi_data','R','phi_c_modA', ...
                 'phi_c_Noir_old','half_L','chi_f','d','L','e_n', ...
                 'xi_eff','xi_smooth','fnt','tau_lambda_th','tau_E_th', ...
                 'v_T_ref','lambda_T_modA_at_DEM','n_data','dist_from_wall'};
missing = required_vars(~cellfun(@(v) ...
    evalin('base', sprintf('exist(''%s'',''var'')',v))==1, required_vars));
if ~isempty(missing)
    error(['Run Phase_v4.m first. Missing variables: ', strjoin(missing,', ')]);
end

% ============== Color palette (matches Phase_v4.m) ====================
col_th    = [0.85 0.10 0.10];
col_th2   = [0.20 0.20 0.85];
col_modA  = [0.85 0.55 0.10];
col_ball  = [0.40 0.20 0.55];
col_tot   = [0.00 0.00 0.00];
col_dem   = [0.10 0.55 0.20];
col_dem2  = [0.85 0.35 0.10];
col_grid  = [0.45 0.45 0.45];
col_amber = [0.85 0.65 0.13];

% Convenient masks
valid_m_full_tr  = R.fit_mode == 1 & ~R.lambda_is_lb &  R.trust_fit;
valid_m_outer_tr = R.fit_mode == 2 &  R.trust_fit;
valid_t          = R.tau_lambda > 0 & isfinite(R.tau_lambda) & ...
                   ~R.lambda_is_lb & R.trust_fit;
valid_s_arr      = R.sigma2_n > 0;
ok_ratio         = R.lambda_DEM > 0 & ~R.lambda_is_lb;
ok_tc            = R.Tc_over_Tw > 0;

x_lo = min(phi_data) * 0.7;
x_hi = max(phi_data) * 1.4;

%% =====================================================================
%  Figure 1: clean v4-only length-scale plot
%% =====================================================================
if exist('fig1','var') && isgraphics(fig1), close(fig1); end
fig1 = figure('Name','Fig1_v4_PubGrade','Position',[60 60 760 560],'Color','w');

h_lb = loglog(phi_th, lambda_b      , '-.', ...
              'Color', col_ball,  'LineWidth', 1.8); hold on;
h_lT = loglog(phi_th, lambda_T_modA , '-' , ...
              'Color', col_modA,  'LineWidth', 2.2);
h_lt = loglog(phi_th, lambda_tot_modA,'-' , ...
              'Color', col_tot,   'LineWidth', 3.2);

h_dem1 = draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, ...
                           valid_m_full_tr , 'o', col_dem);
h_dem2 = draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, ...
                           valid_m_outer_tr, 's', col_dem2);

yline(half_L   , ':' , 'Color', col_grid, 'LineWidth', 2.0, ...
      'HandleVisibility','off');
xline(phi_c_modA,'-.', 'Color', col_tot , 'LineWidth', 1.5, ...
      'HandleVisibility','off');

yl1 = [3e-2 5e1];
patch([phi_c_modA 0.5 0.5 phi_c_modA], ...
      [yl1(1) yl1(1) half_L half_L], ...
      [1 0.93 0.93], 'EdgeColor','none','FaceAlpha',0.55, ...
      'HandleVisibility','off');

text(min(phi_th)*1.3, half_L*1.35, '$L/2$', ...
     'Interpreter','latex','FontSize',16,'Color',col_grid);
text(phi_c_modA*1.10, yl1(2)*0.4, ...
     sprintf('$\\phi_c = %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',15,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Length scales [cm]'    , 'Interpreter','latex','FontSize',22);

leg_h = [h_lb, h_lT, h_lt];
leg_l = { ...
    '$\lambda_b$ (ballistic m.f.p.)', ...
    sprintf('$\\lambda_T$ (Model A, $\\xi_{\\rm eff}=%.2f$)', xi_eff), ...
    '$\lambda_{\rm tot}=\lambda_b+\lambda_T$'};
if isgraphics(h_dem1), leg_h(end+1)=h_dem1; leg_l{end+1}='DEM $\lambda_T$ (trusted)'; end
if isgraphics(h_dem2), leg_h(end+1)=h_dem2; leg_l{end+1}='DEM $\lambda_T$ (outer refit)'; end
legend(leg_h, leg_l, 'Interpreter','latex','FontSize',12,'Location','southwest');
legend boxoff;

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([min(phi_th) 0.3]); ylim(yl1); box on;

%% =====================================================================
%  Figure 2: characteristic times τ_λ, τ_E
%% =====================================================================
if exist('fig2','var') && isgraphics(fig2), close(fig2); end
fig2 = figure('Name','Fig2_Tau_PubGrade','Position',[110 110 740 560],'Color','w');

h_tl_th = loglog(phi_th, tau_lambda_th, '-' , ...
                 'Color', col_modA, 'LineWidth', 2.8); hold on;
h_tE_th = loglog(phi_th, tau_E_th     , '--', ...
                 'Color', col_amber, 'LineWidth', 2.0);

h_tl_dem = []; h_tE_dem = [];
if any(valid_t)
    h_tl_dem = loglog(phi_data(valid_t), R.tau_lambda(valid_t), 'o', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_modA);
    h_tE_dem = loglog(phi_data(valid_t), R.tau_E(valid_t), 's', ...
        'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_amber);
end

tau_c = half_L / v_T_ref;
yline(tau_c, ':', 'Color', col_grid, 'LineWidth', 2.0, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.5, 'HandleVisibility','off');

% annotations placed where curves don't go
text(min(phi_th)*1.5, tau_c*1.5, '$\tau_c = (L/2)/v_T$', ...
     'Interpreter','latex','FontSize',14,'Color',col_grid);
text(phi_c_modA*1.08, tau_c*0.10, sprintf('$\\phi_c=%.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Characteristic time [s]', 'Interpreter','latex','FontSize',22);

leg2_h = [h_tl_th, h_tE_th];
leg2_l = {'$\tau_\lambda = \lambda_T/v_T$ (theory)', ...
          '$\tau_E$ (Enskog collision time)'};
if ~isempty(h_tl_dem)
    leg2_h(end+1) = h_tl_dem; leg2_l{end+1} = '$\tau_\lambda$ (DEM)';
end
if ~isempty(h_tE_dem)
    leg2_h(end+1) = h_tE_dem; leg2_l{end+1} = '$\tau_E$ (DEM)';
end
legend(leg2_h, leg2_l, 'Interpreter','latex','FontSize',12, ...
       'Location','southwest');
legend boxoff;

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([min(phi_th) max(phi_th)]); box on;

%% =====================================================================
%  Figure 3: λ_DEM/λ_T,Model_A and σ² (dual-axis)
%% =====================================================================
if exist('fig3','var') && isgraphics(fig3), close(fig3); end
fig3 = figure('Name','Fig3_DualAxis_PubGrade','Position',[160 60 800 560],'Color','w');

% Use Model A in denominator (was smooth in original — inconsistent with v4)
lam_ratio_modA = R.lambda_DEM ./ lambda_T_modA_at_DEM;

yyaxis left
plot(phi_data(ok_ratio), lam_ratio_modA(ok_ratio), 'o-', ...
    'LineWidth', 2.2, 'MarkerSize', 10, ...
    'Color', col_modA, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_modA);
hold on;
yline(1, '--', 'Color', col_modA, 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel('$\lambda_{T,\rm DEM}\,/\,\lambda_{T,\rm Model\ A}$', ...
       'Interpreter','latex','FontSize',20,'Color',col_modA);
ax = gca; ax.YColor = col_modA;
ylim([0 1.8]);
text(min(phi_data)*1.5, 1.10, '\textit{theory holds}', ...
     'Interpreter','latex','FontSize',12,'Color',col_modA);

yyaxis right
semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 's-', ...
    'LineWidth', 2.2, 'MarkerSize', 10, ...
    'Color', col_th, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_th);
yline(1, '--', 'Color', col_th, 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel('$\sigma_n^2/\langle n\rangle$', ...
       'Interpreter','latex','FontSize',20,'Color',col_th);
ax.YColor = col_th;
set(gca, 'YScale', 'log');

xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6, ...
      'HandleVisibility','off');
text(phi_c_modA*1.08, 25, ...
     sprintf('$\\phi_c \\approx %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([x_lo, x_hi]);

legend({'$\lambda_{T,\rm DEM}/\lambda_{T,\rm Model\ A}$', ...
        '$\sigma_n^2/\langle n\rangle$'}, ...
    'Interpreter','latex','FontSize',12,'Location','northwest');
legend boxoff;
box on;

%% =====================================================================
%  Figure 4: clean predictive phase map (φ, e_n)
%% =====================================================================
if exist('fig4','var') && isgraphics(fig4), close(fig4); end
fig4 = figure('Name','Fig4_PhaseMap_PubGrade', ...
              'Position',[210 110 820 600],'Color','w');

phi_g = logspace(log10(0.001), log10(0.3), 360);
en_g  = linspace(0.55, 0.999, 240);
[PHI, EN] = meshgrid(phi_g, en_g);
CHI    = chi_f(PHI);
XI_g   = 1 - EN.^2;
LB_g   = d ./ (6 * PHI .* CHI);
LT_g   = d ./ (6 * PHI .* CHI .* sqrt(XI_g));
LTOT_g = LB_g + LT_g;
LTOT_norm = LTOT_g / L;
LT_norm   = LT_g   / L;

[~, h_bg] = contourf(PHI, EN, log10(LTOT_norm), 30, 'LineColor','none');
set(h_bg, 'HandleVisibility', 'off');
hold on;
colormap(gca, parula);
cb = colorbar;
cb.Label.String      = '$\log_{10}\,(\lambda_b+\lambda_T)/L$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 16;
cb.FontSize          = 13;
cb.FontName          = fnt;

[~, h_new] = contour(PHI, EN, LTOT_norm, [0.5 0.5], '-', ...
                     'Color', [0.85 0.10 0.10], 'LineWidth', 3.2);
[~, h_old] = contour(PHI, EN, LT_norm,   [0.5 0.5], '--', ...
                     'Color', [1.00 1.00 1.00], 'LineWidth', 2.2);

sigma_valid_arr = R.sigma2_n(R.sigma2_n > 0);
if isempty(sigma_valid_arr)
    clust_thr = inf;
else
    clust_thr = max(2.0, 2*median(sigma_valid_arr));
end
is_clust = R.sigma2_n > clust_thr;

h_gas  = plot(phi_data(~is_clust), e_n*ones(sum(~is_clust),1), 'o', ...
    'MarkerSize', 10, 'LineWidth', 2.0, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
h_clst = plot(phi_data( is_clust), e_n*ones(sum( is_clust),1), 's', ...
    'MarkerSize', 11, 'LineWidth', 2.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.85 0.10 0.10]);

set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on', 'Layer','top');
xlabel('Volume fraction $\phi$' , 'Interpreter','latex','FontSize',22);
ylabel('Normal restitution $e_n$', 'Interpreter','latex','FontSize',22);
xlim([min(phi_g) max(phi_g)]);
ylim([min(en_g) max(en_g)]);
box on;

legend([h_new, h_old, h_gas, h_clst], { ...
    '$(\lambda_b+\lambda_T)/L = 1/2$  (v4, new)', ...
    '$\lambda_T/L = 1/2$  (v3, old)', ...
    'DEM: homogeneous gas', ...
    'DEM: clustered'}, ...
    'Interpreter','latex','FontSize',12, ...
    'Location','southwest', ...
    'Color',[1 1 1 0.92], 'EdgeColor','none');

%% =====================================================================
%  Figure 5: σ²/⟨n⟩ order parameter
%% =====================================================================
if exist('fig5','var') && isgraphics(fig5), close(fig5); end
fig5 = figure('Name','Fig5_OrderParam_PubGrade','Position',[260 60 720 540],'Color','w');

semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.6, 'MarkerSize', 11, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;

yline(1, '--', 'Color', col_grid, 'LineWidth', 1.8, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.8, 'HandleVisibility','off');

% Position labels safely inside log-y axes
yl5 = ylim;
text(min(phi_data)*1.5, 1.18, 'Poisson (ideal gas)', ...
     'FontSize', 13, 'Color', col_grid);
text(phi_c_modA*1.08, yl5(2)*0.4, ...
     sprintf('$\\phi_c \\approx %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

% Light shading on the clustered side
patch_y_top = yl5(2);
patch([phi_c_modA, max(phi_data)*2, max(phi_data)*2, phi_c_modA], ...
      [yl5(1), yl5(1), patch_y_top, patch_y_top], ...
      [1 0.95 0.95], 'EdgeColor','none','FaceAlpha',0.4, ...
      'HandleVisibility','off');

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter','latex','FontSize',22);

set(gca,'XScale','log','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([x_lo, x_hi]);
box on;

%% =====================================================================
%  Figure 6: T(z)/T_wall profiles with parula colorbar
%% =====================================================================
if exist('fig6','var') && isgraphics(fig6), close(fig6); end
fig6 = figure('Name','Fig6_Tz_PubGrade','Position',[310 110 800 580],'Color','w');

% Use parula instead of jet; thin to readable count of profiles
n_show = min(n_data, 14);   % show at most ~14 profiles
show_idx = round(linspace(1, n_data, n_show));
cmap6 = parula(n_data);

for k = 1:n_show
    i = show_idx(k);
    T_prof = R.T_shell(i,:);
    vs = T_prof > 0;
    if ~any(vs), continue; end
    Tw_norm = T_prof(find(vs,1,'last'));
    col_i = cmap6(i,:);

    plot(dist_from_wall(vs), T_prof(vs)/Tw_norm, 'o', ...
        'MarkerSize', 5, 'Color', col_i, 'MarkerFaceColor', col_i, ...
        'HandleVisibility','off'); hold on;

    if R.lambda_DEM(i) > 0
        zd = linspace(0, half_L, 200);
        plot(zd, exp(-zd/R.lambda_DEM(i)), '-', ...
            'Color', col_i, 'LineWidth', 1.3, 'HandleVisibility','off');
    end
end

yline(1, ':', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');

xlabel('Distance from wall $z$ [cm]', 'Interpreter','latex','FontSize',22);
ylabel('$T(z) / T_{\rm wall}$', 'Interpreter','latex','FontSize',22);

% Colorbar coding φ (uses log10 for visual uniformity on log-spaced φ)
log_phi = log10(phi_data);
caxis([log_phi(1) log_phi(end)]);
colormap(parula);
cb6 = colorbar;
cb6.Label.String      = '$\log_{10}\,\phi$';
cb6.Label.Interpreter = 'latex';
cb6.Label.FontSize    = 16;
cb6.FontSize          = 13;
cb6.FontName          = fnt;

% Custom phi ticks
n_ct = 5;
tick_pos = linspace(log_phi(1), log_phi(end), n_ct);
cb6.Ticks = tick_pos;
cb6.TickLabels = arrayfun(@(x) sprintf('%.3f', 10^x), tick_pos, ...
                          'UniformOutput', false);
cb6.Label.String      = '$\phi$';

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
ylim([0 1.2]); xlim([0 half_L]); box on;

%% =====================================================================
%  Figure 7: four-witness 2×2 panel
%% =====================================================================
if exist('fig7','var') && isgraphics(fig7), close(fig7); end
fig7 = figure('Name','Fig7_FourWitnesses_PubGrade', ...
              'Position',[60 60 1100 800],'Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% --- (a) Length scales: v4-only, matching new Fig 1 ----------------------
ax_a = nexttile;
loglog(phi_th, lambda_b      , '-.', 'Color', col_ball , 'LineWidth', 1.6); hold on;
loglog(phi_th, lambda_T_modA , '-' , 'Color', col_modA , 'LineWidth', 2.0);
loglog(phi_th, lambda_tot_modA,'-' , 'Color', col_tot  , 'LineWidth', 2.8);
draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_full_tr , 'o', col_dem);
draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_outer_tr, 's', col_dem2);

yline(half_L , ':' , 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot , 'LineWidth', 1.4, 'HandleVisibility','off');
ylabel('Lengths [cm]', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(a)}', 'Units','normalized', ...
     'Interpreter','latex','FontSize',18);
text(0.04, 0.78, '$\lambda_b+\lambda_T=L/2$', 'Units','normalized', ...
     'Interpreter','latex','FontSize',11,'Color',col_grid);
set(gca,'FontSize',14,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
ylim_a = [1e-1 1e2];
xlim([x_lo x_hi]); ylim(ylim_a); box on;

% --- (b) T_center / T_wall vs φ ------------------------------------------
ax_b = nexttile;
semilogx(phi_data(ok_tc), R.Tc_over_Tw(ok_tc), 'o-', ...
    'LineWidth', 2.0, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
ylabel('$T_{\rm center}/T_{\rm wall}$', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(b)}', 'Units','normalized', ...
     'Interpreter','latex','FontSize',18);
text(0.04, 0.78, 'equilibrium', 'Units','normalized', ...
     'FontSize',11,'Color',col_grid);
set(gca,'FontSize',14,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([x_lo x_hi]); ylim([0 1.2]); box on;

% --- (c) σ²/⟨n⟩ vs φ -----------------------------------------------------
ax_c = nexttile;
loglog(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.0, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem);
hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\sigma_n^2/\langle n\rangle$', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(c)}', 'Units','normalized', ...
     'Interpreter','latex','FontSize',18);
text(0.04, 0.78, 'Poisson', 'Units','normalized', ...
     'FontSize',11,'Color',col_grid);
set(gca,'FontSize',14,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([x_lo x_hi]); box on;

% --- (d) τ_λ vs φ (replacing the τ_E-only original) ----------------------
ax_d = nexttile;
loglog(phi_th, tau_lambda_th, '-' , 'Color', col_modA, 'LineWidth', 2.0); hold on;
if any(valid_t)
    loglog(phi_data(valid_t), R.tau_lambda(valid_t), 'o', ...
        'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_modA);
end
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\tau_\lambda$ [s]', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(d)}', 'Units','normalized', ...
     'Interpreter','latex','FontSize',18);
set(gca,'FontSize',14,'FontName',fnt,'LineWidth',1.5, 'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([x_lo x_hi]); box on;

title(tlo, 'Four independent witnesses of the gas-to-clustering transition', ...
      'FontSize', 16, 'FontName', fnt, 'FontWeight', 'normal');

%% =====================================================================
%  Figure 8: linear-scale criterion plot
%% =====================================================================
if exist('fig8','var') && isgraphics(fig8), close(fig8); end
fig8 = figure('Name','Fig8_Linear_PubGrade','Position',[360 60 820 580],'Color','w');

h8_lT  = semilogx(phi_th, lambda_T_modA/L , '-' , ...
                  'Color', col_modA, 'LineWidth', 2.4); hold on;
h8_lb  = semilogx(phi_th, lambda_b/L      , '-.', ...
                  'Color', col_ball, 'LineWidth', 2.4);
h8_lt  = semilogx(phi_th, lambda_tot_modA/L,'-' , ...
                  'Color', col_tot , 'LineWidth', 3.2);

yline(0.5, ':', 'Color', col_grid, 'LineWidth', 2.5, 'HandleVisibility','off');

% v3 vs v4 critical phi
xline(phi_c_Noir_old, ':' , 'Color', col_th , 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA    , '-.', 'Color', col_tot, 'LineWidth', 2.0, 'HandleVisibility','off');

% Shade clustered region
patch([phi_c_modA 0.1 0.1 phi_c_modA], [0 0 0.5 0.5], ...
      [1 0.92 0.92], 'EdgeColor','none','FaceAlpha',0.50, ...
      'HandleVisibility','off');

% Annotations placed in clearly empty regions
text(0.0014, 0.55, 'criterion $\lambda_{\rm tot}/L = 1/2$', ...
     'Interpreter','latex','FontSize',13,'Color',col_grid);
text(phi_c_Noir_old*0.55, 1.08, ...
     sprintf('$\\phi_c^{\\rm v3}=%.3f$', phi_c_Noir_old), ...
     'Interpreter','latex','FontSize',12,'Color',col_th);
text(phi_c_modA*1.05, 1.08, ...
     sprintf('$\\phi_c^{\\rm v4}=%.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',12,'Color',col_tot);

% Phase labels — moved to corners where curves are far away
text(0.0018, 0.18, 'Homogeneous Gas', ...
     'FontSize', 13, 'Color', [0.10 0.40 0.10], ...
     'FontWeight','bold');
text(0.055, 0.18, 'Clustering', ...
     'FontSize', 13, 'Color', [0.60 0.10 0.10], ...
     'FontWeight','bold');

% σ²=2 onset (only if it's far enough from φ_c to avoid clutter)
if any(R.sigma2_n > 0)
    [~, idx_onset] = min(abs(R.sigma2_n(R.sigma2_n > 0) - 2.0));
    phi_with_sigma = phi_data(R.sigma2_n > 0);
    phi_sigma_onset = phi_with_sigma(idx_onset);
    if abs(log10(phi_sigma_onset) - log10(phi_c_modA)) > 0.05
        xline(phi_sigma_onset, '--', 'Color', col_dem, 'LineWidth', 1.5, ...
              'HandleVisibility','off');
        text(phi_sigma_onset*1.05, 0.92, ...
             sprintf('$\\phi(\\sigma^2{=}2)\\approx%.3f$', phi_sigma_onset), ...
             'Interpreter','latex','FontSize',11,'Color',col_dem);
    end
end

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Normalized length $\lambda / L$', 'Interpreter','latex','FontSize',22);

legend([h8_lT, h8_lb, h8_lt], { ...
    sprintf('$\\lambda_T/L$ (Model A, $\\xi_{\\rm eff}=%.2f$)', xi_eff), ...
    '$\lambda_b/L$ (ballistic m.f.p.)', ...
    '$(\lambda_b+\lambda_T)/L$ (v4 criterion)'}, ...
    'Interpreter','latex','FontSize',12,'Location','northeast');
legend boxoff;

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on', 'Layer','top');
xlim([min(phi_th) 0.1]); ylim([0 1.2]); box on;

fprintf('All 8 publication-grade figures produced.\n');

%% =====================================================================
%  Local helper function (must be at end of script in MATLAB ≥ R2016b)
%% =====================================================================
function h = draw_errorbar_pts(phi_v, lam_v, lci_v, mask, marker, faceColor)
    h = [];
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

    h = errorbar(phi_s, lam_s, lneg, lpos, ...
        marker, 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', faceColor, ...
        'CapSize', 5, 'LineStyle', 'none');
end