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

%% ========== 5. PLOTTING SETUP =======================================
% Color palette
col_th     = [0.85 0.10 0.10];   % normal-only reference (red)
col_th2    = [0.20 0.20 0.85];   % (reserved)
col_modA   = [0.85 0.55 0.10];   % Model A (orange)
col_dem    = [0.10 0.55 0.20];   % DEM trusted (green)
col_dem2   = [0.85 0.35 0.10];   % DEM outer-refit (red-orange)
col_grid   = [0.45 0.45 0.45];
col_untrust= [0.65 0.65 0.65];
col_ball   = [0.40 0.20 0.55];   % lambda_b (purple)
col_tot    = [0.00 0.00 0.00];   % lambda_total (black)
col_amber  = [0.85 0.65 0.13];
fnt = 'Times New Roman';

% Masks
valid_m_full_tr  = R.fit_mode == 1 & ~R.lambda_is_lb &  R.trust_fit;
valid_m_full_un  = R.fit_mode == 1 & ~R.lambda_is_lb & ~R.trust_fit;
valid_m_outer_tr = R.fit_mode == 2 &  R.trust_fit;
valid_m_outer_un = R.fit_mode == 2 & ~R.trust_fit;
valid_m_lb       = R.lambda_is_lb;
valid_t          = R.tau_lambda > 0 & isfinite(R.tau_lambda) & ...
                   ~R.lambda_is_lb & R.trust_fit;
valid_s_arr      = R.sigma2_n > 0;
ok_ratio         = R.lambda_DEM > 0 & ~R.lambda_is_lb;
ok_tc            = R.Tc_over_Tw > 0;

% Knudsen-validity boundary (kept for the summary/save block)
phi_Kn1 = phi_th(find(lambda_T_Noir <= L*Kn_thresh, 1, 'first'));
if isempty(phi_Kn1), phi_Kn1 = NaN; end

% ---- Effective rotational dissipation, extracted from the fit --------
% xi_eff (fitted to lambda_T,DEM) is e_n-independent. Decompose it into a
% normal channel (1-e_n^2) and a rotational/tangential channel:
%     xi_eff = xi_normal + xi_rot_fit
% With the corrected near-elastic e_n=0.99, xi_normal is tiny, so almost
% all of the measured dissipation lives in the rotational channel. This
% xi_rot_fit (assumed e_n-independent) is what the phase map uses to draw
% boundaries across the (phi, e_n) plane.
xi_rot_fit = max(xi_eff - (1 - e_n^2), 0);
fprintf('Rotational dissipation (effective): xi_rot_fit = %.3f\n', xi_rot_fit);
fprintf('  -> normal channel  1-e_n^2 = %.4f  (%.0f%% of xi_eff)\n', ...
        1-e_n^2, 100*(1-e_n^2)/xi_eff);
fprintf('  -> rotational channel        = %.4f  (%.0f%% of xi_eff)\n\n', ...
        xi_rot_fit, 100*xi_rot_fit/xi_eff);

% Linear x-axis viewing window (focused on the data + transition)
phi_view_lo = 0.002;
phi_view_hi = 0.07;

%% ---------- Figure 1: lambda_T, lambda_b, lambda_tot (LINEAR x) ------
fig1 = figure('Name','Fig1_lengths','Position',[60 60 760 560],'Color','w');

h_lb = semilogy(phi_th, lambda_b      , '-.', ...
                'Color', col_ball,  'LineWidth', 1.8); hold on;
h_lT = semilogy(phi_th, lambda_T_modA , '-' , ...
                'Color', col_modA,  'LineWidth', 2.2);
h_lt = semilogy(phi_th, lambda_tot_modA,'-' , ...
                'Color', col_tot,   'LineWidth', 3.2);

h_dem1 = draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, ...
                           valid_m_full_tr , 'o', col_dem);
h_dem2 = draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, ...
                           valid_m_outer_tr, 's', col_dem2);

yline(half_L   , ':' , 'Color', col_grid, 'LineWidth', 2.0, 'HandleVisibility','off');
xline(phi_c_modA,'-.', 'Color', col_tot , 'LineWidth', 1.5, 'HandleVisibility','off');

yl1 = [3e-2 5e1];
patch([phi_c_modA phi_view_hi phi_view_hi phi_c_modA], ...
      [yl1(1) yl1(1) half_L half_L], ...
      [1 0.93 0.93], 'EdgeColor','none','FaceAlpha',0.55, 'HandleVisibility','off');

text(phi_view_lo+0.0015, half_L*1.35, '$L/2$', ...
     'Interpreter','latex','FontSize',16,'Color',col_grid);
text(phi_c_modA+0.0015, yl1(2)*0.4, sprintf('$\\phi_c = %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',15,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Length scales [cm]'    , 'Interpreter','latex','FontSize',22);

leg_h = [h_lb, h_lT, h_lt];
leg_l = { '$\lambda_b$ (ballistic m.f.p.)', ...
    sprintf('$\\lambda_T$ (Model A, $\\xi_{\\rm eff}=%.2f$)', xi_eff), ...
    '$\lambda_{\rm tot}=\lambda_b+\lambda_T$'};
if isgraphics(h_dem1), leg_h(end+1)=h_dem1; leg_l{end+1}='DEM $\lambda_T$ (trusted)'; end
if isgraphics(h_dem2), leg_h(end+1)=h_dem2; leg_l{end+1}='DEM $\lambda_T$ (outer refit)'; end
legend(leg_h, leg_l, 'Interpreter','latex','FontSize',12,'Location','northeast');
legend boxoff;

set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); ylim(yl1); box on;

%% ---------- Figure 2: characteristic times (LINEAR x) ---------------
fig2 = figure('Name','Fig2_Tau','Position',[110 110 740 560],'Color','w');

h_tl_th = semilogy(phi_th, tau_lambda_th, '-' , ...
                   'Color', col_modA, 'LineWidth', 2.8); hold on;
h_tE_th = semilogy(phi_th, tau_E_th     , '--', ...
                   'Color', col_amber, 'LineWidth', 2.0);

h_tl_dem = []; h_tE_dem = [];
if any(valid_t)
    h_tl_dem = semilogy(phi_data(valid_t), R.tau_lambda(valid_t), 'o', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_modA);
    h_tE_dem = semilogy(phi_data(valid_t), R.tau_E(valid_t), 's', ...
        'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', col_amber);
end

tau_c = half_L / v_T_ref;
yline(tau_c, ':', 'Color', col_grid, 'LineWidth', 2.0, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.5, 'HandleVisibility','off');

text(phi_view_lo+0.0015, tau_c*1.6, '$\tau_c = (L/2)/v_T$', ...
     'Interpreter','latex','FontSize',14,'Color',col_grid);
text(phi_c_modA+0.0015, tau_c*0.10, sprintf('$\\phi_c=%.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('Characteristic time [s]', 'Interpreter','latex','FontSize',22);

leg2_h = [h_tl_th, h_tE_th];
leg2_l = {'$\tau_\lambda = \lambda_T/v_T$ (theory)', ...
          '$\tau_E$ (Enskog collision time)'};
if ~isempty(h_tl_dem), leg2_h(end+1)=h_tl_dem; leg2_l{end+1}='$\tau_\lambda$ (DEM)'; end
if ~isempty(h_tE_dem), leg2_h(end+1)=h_tE_dem; leg2_l{end+1}='$\tau_E$ (DEM)'; end
legend(leg2_h, leg2_l, 'Interpreter','latex','FontSize',12,'Location','northeast');
legend boxoff;

set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); box on;

%% ---------- Figure 3: lambda_DEM/lambda_modA + sigma^2 (LINEAR x) ---
fig3 = figure('Name','Fig3_DualAxis','Position',[160 60 800 560],'Color','w');

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
text(phi_view_lo+0.002, 1.10, '\textit{theory holds}', ...
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

xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.6, 'HandleVisibility','off');
text(phi_c_modA+0.0015, 25, sprintf('$\\phi_c \\approx %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]);

legend({'$\lambda_{T,\rm DEM}/\lambda_{T,\rm Model\ A}$', ...
        '$\sigma_n^2/\langle n\rangle$'}, ...
    'Interpreter','latex','FontSize',12,'Location','northwest');
legend boxoff; box on;

%% ---------- Figure 4: predictive phase map (LINEAR x) ---------------
% IMPORTANT (physics): the phase boundary uses the FULL effective
% dissipation xi = (1-e_n^2) + xi_rot_fit, NOT the normal-only 1-e_n^2.
% With near-elastic e_n=0.99 the normal channel alone (1-e_n^2 ~ 0.02)
% predicts the boundary far to the right and misclassifies the DEM points;
% including the fitted rotational channel restores agreement.
fig4 = figure('Name','Fig4_PhaseMap','Position',[210 110 820 600],'Color','w');

phi_g = linspace(phi_view_lo, 0.08, 400);
en_g  = linspace(0.55, 0.999, 260);
[PHI, EN] = meshgrid(phi_g, en_g);
CHI    = chi_f(PHI);
XI_g   = (1 - EN.^2) + xi_rot_fit;          % full effective dissipation
LB_g   = d ./ (6 * PHI .* CHI);
LT_g   = d ./ (6 * PHI .* CHI .* sqrt(XI_g));
LTOT_g = LB_g + LT_g;
LTOT_norm = LTOT_g / L;
LT_norm   = LT_g   / L;

[~, h_bg] = contourf(PHI, EN, log10(LTOT_norm), 30, 'LineColor','none');
set(h_bg, 'HandleVisibility', 'off'); hold on;
colormap(gca, parula);
cb = colorbar;
cb.Label.String      = '$\log_{10}\,(\lambda_b+\lambda_T)/L$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 16; cb.FontSize = 13; cb.FontName = fnt;

[~, h_new] = contour(PHI, EN, LTOT_norm, [0.5 0.5], '-', ...
                     'Color', [0.85 0.10 0.10], 'LineWidth', 3.2);
[~, h_old] = contour(PHI, EN, LT_norm,   [0.5 0.5], '--', ...
                     'Color', [1.00 1.00 1.00], 'LineWidth', 2.2);

sigma_valid_arr = R.sigma2_n(R.sigma2_n > 0);
if isempty(sigma_valid_arr), clust_thr = inf;
else, clust_thr = max(2.0, 2*median(sigma_valid_arr)); end
is_clust = R.sigma2_n > clust_thr;

h_gas  = plot(phi_data(~is_clust), e_n*ones(sum(~is_clust),1), 'o', ...
    'MarkerSize', 10, 'LineWidth', 2.0, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
h_clst = plot(phi_data( is_clust), e_n*ones(sum( is_clust),1), 's', ...
    'MarkerSize', 11, 'LineWidth', 2.2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.85 0.10 0.10]);

set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','XMinorTick','on','YMinorTick','on','Layer','top');
xlabel('Volume fraction $\phi$' , 'Interpreter','latex','FontSize',22);
ylabel('Normal restitution $e_n$', 'Interpreter','latex','FontSize',22);
xlim([min(phi_g) max(phi_g)]); ylim([min(en_g) max(en_g)]); box on;

legend([h_new, h_old, h_gas, h_clst], { ...
    '$(\lambda_b+\lambda_T)/L = 1/2$  (v4, new)', ...
    '$\lambda_T/L = 1/2$  (v3, old)', ...
    'DEM: homogeneous gas', 'DEM: clustered'}, ...
    'Interpreter','latex','FontSize',12,'Location','southeast', ...
    'Color',[1 1 1 0.92],'EdgeColor','none');

%% ---------- Figure 5: sigma^2 order parameter (LINEAR x) ------------
fig5 = figure('Name','Fig5_OrderParam','Position',[260 60 720 540],'Color','w');

semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.6, 'MarkerSize', 11, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;

yline(1, '--', 'Color', col_grid, 'LineWidth', 1.8, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.8, 'HandleVisibility','off');

yl5 = ylim;
text(phi_view_lo+0.002, 1.18, 'Poisson (ideal gas)', 'FontSize', 13, 'Color', col_grid);
text(phi_c_modA+0.0015, yl5(2)*0.4, sprintf('$\\phi_c \\approx %.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',14,'Color',col_tot);

patch([phi_c_modA, phi_view_hi, phi_view_hi, phi_c_modA], ...
      [yl5(1), yl5(1), yl5(2), yl5(2)], ...
      [1 0.95 0.95], 'EdgeColor','none','FaceAlpha',0.4, 'HandleVisibility','off');

xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',22);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter','latex','FontSize',22);

set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); box on;

%% ---------- Figure 6: T(z)/T_wall profiles (x already linear) -------
fig6 = figure('Name','Fig6_Tz','Position',[310 110 800 580],'Color','w');

n_show = min(n_data, 14);
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
        'MarkerSize', 5, 'Color', col_i, 'MarkerFaceColor', col_i, 'HandleVisibility','off'); hold on;
    if R.lambda_DEM(i) > 0
        zd = linspace(0, half_L, 200);
        plot(zd, exp(-zd/R.lambda_DEM(i)), '-', ...
            'Color', col_i, 'LineWidth', 1.3, 'HandleVisibility','off');
    end
end

yline(1, ':', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Distance from wall $z$ [cm]', 'Interpreter','latex','FontSize',22);
ylabel('$T(z) / T_{\rm wall}$', 'Interpreter','latex','FontSize',22);

log_phi = log10(phi_data);
caxis([log_phi(1) log_phi(end)]);
colormap(parula);
cb6 = colorbar;
cb6.Label.Interpreter = 'latex'; cb6.Label.FontSize = 16;
cb6.FontSize = 13; cb6.FontName = fnt;
n_ct = 5;
tick_pos = linspace(log_phi(1), log_phi(end), n_ct);
cb6.Ticks = tick_pos;
cb6.TickLabels = arrayfun(@(x) sprintf('%.3f', 10^x), tick_pos, 'UniformOutput', false);
cb6.Label.String = '$\phi$';

set(gca,'FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
ylim([0 1.2]); xlim([0 half_L]); box on;

%% ---------- Figure 7: four-witness 2x2 panel (LINEAR x) -------------
fig7 = figure('Name','Fig7_FourWitnesses','Position',[60 60 1100 800],'Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) length scales
nexttile;
semilogy(phi_th, lambda_b      , '-.', 'Color', col_ball , 'LineWidth', 1.6); hold on;
semilogy(phi_th, lambda_T_modA , '-' , 'Color', col_modA , 'LineWidth', 2.0);
semilogy(phi_th, lambda_tot_modA,'-' , 'Color', col_tot  , 'LineWidth', 2.8);
draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_full_tr , 'o', col_dem);
draw_errorbar_pts(phi_data, R.lambda_DEM, R.lambda_DEM_CI, valid_m_outer_tr, 's', col_dem2);
yline(half_L , ':' , 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot , 'LineWidth', 1.4, 'HandleVisibility','off');
ylabel('Lengths [cm]', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(a)}', 'Units','normalized','Interpreter','latex','FontSize',18);
text(0.04, 0.78, '$\lambda_b+\lambda_T=L/2$', 'Units','normalized', ...
     'Interpreter','latex','FontSize',11,'Color',col_grid);
set(gca,'XScale','linear','FontSize',14,'FontName',fnt,'LineWidth',1.5,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); ylim([1e-1 1e2]); box on;

% (b) T_center/T_wall
nexttile;
plot(phi_data(ok_tc), R.Tc_over_Tw(ok_tc), 'o-', ...
    'LineWidth', 2.0, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
ylabel('$T_{\rm center}/T_{\rm wall}$', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(b)}', 'Units','normalized','Interpreter','latex','FontSize',18);
text(0.04, 0.78, 'equilibrium', 'Units','normalized','FontSize',11,'Color',col_grid);
set(gca,'XScale','linear','FontSize',14,'FontName',fnt,'LineWidth',1.5,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); ylim([0 1.2]); box on;

% (c) sigma^2
nexttile;
semilogy(phi_data(valid_s_arr), R.sigma2_n(valid_s_arr), 'o-', ...
    'LineWidth', 2.0, 'MarkerSize', 9, ...
    'Color', col_dem, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_dem); hold on;
yline(1, '--', 'Color', col_grid, 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\sigma_n^2/\langle n\rangle$', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(c)}', 'Units','normalized','Interpreter','latex','FontSize',18);
text(0.04, 0.78, 'Poisson', 'Units','normalized','FontSize',11,'Color',col_grid);
set(gca,'XScale','linear','FontSize',14,'FontName',fnt,'LineWidth',1.5,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); box on;

% (d) tau_lambda
nexttile;
semilogy(phi_th, tau_lambda_th, '-' , 'Color', col_modA, 'LineWidth', 2.0); hold on;
if any(valid_t)
    semilogy(phi_data(valid_t), R.tau_lambda(valid_t), 'o', ...
        'MarkerSize', 9, 'LineWidth', 1.5, 'MarkerEdgeColor','k', 'MarkerFaceColor', col_modA);
end
xline(phi_c_modA, '-.', 'Color', col_tot, 'LineWidth', 1.4, 'HandleVisibility','off');
xlabel('Volume fraction $\phi$', 'Interpreter','latex','FontSize',18);
ylabel('$\tau_\lambda$ [s]', 'Interpreter','latex','FontSize',18);
text(0.04, 0.92, '\textbf{(d)}', 'Units','normalized','Interpreter','latex','FontSize',18);
set(gca,'XScale','linear','FontSize',14,'FontName',fnt,'LineWidth',1.5,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([phi_view_lo phi_view_hi]); box on;

title(tlo, 'Four independent witnesses of the gas-to-clustering transition', ...
      'FontSize', 16, 'FontName', fnt, 'FontWeight', 'normal');

%% ---------- Figure 8: linear-scale criterion plot (LINEAR x) -------
fig8 = figure('Name','Fig8_Criterion','Position',[360 60 820 580],'Color','w');

h8_lT  = plot(phi_th, lambda_T_modA/L , '-' , 'Color', col_modA, 'LineWidth', 2.4); hold on;
h8_lb  = plot(phi_th, lambda_b/L      , '-.', 'Color', col_ball, 'LineWidth', 2.4);
h8_lt  = plot(phi_th, lambda_tot_modA/L,'-' , 'Color', col_tot , 'LineWidth', 3.2);

yline(0.5, ':', 'Color', col_grid, 'LineWidth', 2.5, 'HandleVisibility','off');

% v3 (single-length, Model A) vs v4 critical phi — both use xi_eff
xline(phi_c_modA_old, ':' , 'Color', col_th , 'LineWidth', 1.5, 'HandleVisibility','off');
xline(phi_c_modA    , '-.', 'Color', col_tot, 'LineWidth', 2.0, 'HandleVisibility','off');

patch([phi_c_modA phi_view_hi phi_view_hi phi_c_modA], [0 0 0.5 0.5], ...
      [1 0.92 0.92], 'EdgeColor','none','FaceAlpha',0.50, 'HandleVisibility','off');

text(phi_view_lo+0.001, 0.55, 'criterion $\lambda_{\rm tot}/L = 1/2$', ...
     'Interpreter','latex','FontSize',13,'Color',col_grid);
if isfinite(phi_c_modA_old)
    text(phi_c_modA_old-0.006, 1.08, sprintf('$\\phi_c^{\\rm v3}=%.3f$', phi_c_modA_old), ...
         'Interpreter','latex','FontSize',12,'Color',col_th);
end
text(phi_c_modA+0.0015, 1.08, sprintf('$\\phi_c^{\\rm v4}=%.3f$', phi_c_modA), ...
     'Interpreter','latex','FontSize',12,'Color',col_tot);

text(phi_view_lo+0.002, 0.18, 'Homogeneous Gas', ...
     'FontSize', 13, 'Color', [0.10 0.40 0.10], 'FontWeight','bold');
text(phi_c_modA+0.018, 0.18, 'Clustering', ...
     'FontSize', 13, 'Color', [0.60 0.10 0.10], 'FontWeight','bold');

if any(R.sigma2_n > 0)
    [~, idx_onset] = min(abs(R.sigma2_n(R.sigma2_n > 0) - 2.0));
    phi_with_sigma = phi_data(R.sigma2_n > 0);
    phi_sigma_onset = phi_with_sigma(idx_onset);
    if abs(phi_sigma_onset - phi_c_modA) > 0.002
        xline(phi_sigma_onset, '--', 'Color', col_dem, 'LineWidth', 1.5, 'HandleVisibility','off');
        text(phi_sigma_onset+0.0015, 0.92, ...
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

set(gca,'XScale','linear','FontSize',18,'FontName',fnt,'LineWidth',1.8, ...
    'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','Layer','top');
xlim([0 phi_view_hi]); ylim([0 1.2]); box on;

fprintf('All 8 figures produced (linear x-axis).\n');

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

%% ========== LOCAL FUNCTIONS (must be at end of script) ===============
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