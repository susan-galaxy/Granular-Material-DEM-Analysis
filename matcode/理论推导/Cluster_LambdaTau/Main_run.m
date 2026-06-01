% =====================================================================
% Cluster_LambdaTau / Main_run.m
%
% One-shot driver that
%   (1) loads all DEM .mat files produced by Pos_Vel_Extrace_v3.m,
%   (2) extracts granular-gas observables run-by-run,
%   (3) plots three publication-quality figures :
%         Fig1 : characteristic length  lambda(phi)  - theory vs DEM
%         Fig2 : phase diagram in (phi, e_n) and (lambda*, tau*) planes
%         Fig3 : prediction figure - T_center/T_wall and T(z) profiles
%
% Edit the USER SETTINGS block, then run.
% =====================================================================

clear; close all; clc;

% Make sure helpers are on the path
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);

%% ============== USER SETTINGS =======================================
% --- where the Data_*.mat files live ---------------------------------
DATA_DIR = 'D:\Desktop\mat\matcode\Cluster_LambdaTau\demo_data';
% Examples :
%   DATA_DIR = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Pos_Vel_data/NoDissNofric';
%   DATA_DIR = 'D:\Desktop\mat\matcode\YourExtractedData';

% --- physical / numerical params (must match your DEM run) -----------
L        = 4.0;                 % cube side, cm
d        = 0.08;                % particle diameter, cm
rho      = 6.0;                 % g/cm^3 (ZrO2)
en       = 0.94;                % normal restitution coeff
TAG      = '';                  % 'Tag_2' (auto-detect if empty)

% --- where to save the figures ---------------------------------------
SAVE_DIR = fullfile(thisDir, 'figures');

% --- theory prefactor : will be auto-fitted in Fig1 if FIT_C is true
FIT_C    = true;
C_SCALE  = 1.0;                 % used only if FIT_C = false

%% ============== 1. Load + extract ===================================
fprintf('\n=== Loading DEM data from %s ===\n', DATA_DIR);
if ~exist(DATA_DIR, 'dir')
    warning('DATA_DIR does not exist : %s\nFalling back to synthetic demo data.', DATA_DIR);
    results = make_synthetic_results(L, d, rho, en);
    SAVE_DIR = fullfile(thisDir, 'figures_demo');
else
    results = dem_load_all(DATA_DIR, ...
        'L',        L, ...
        'd_force',  d, ...
        'rho',      rho, ...
        'en',       en, ...
        'TagField', TAG, ...
        'verbose',  true);
end

if isempty(results)
    error('No usable results extracted from %s.', DATA_DIR);
end

%% ============== 2. Sanity print =====================================
fprintf('\n%-30s %-10s %-12s %-12s %-12s %-10s\n', ...
    'file', 'phi', 'lam_dem(cm)', 'T_mean', 'Tc/Tw', 'cluster?');
for k = 1:numel(results)
    r = results(k);
    [~, fn] = fileparts(r.matfile);
    Tratio  = r.T_center / max(r.T_wall, eps);
    fprintf('%-30s %-10.4f %-12.3g %-12.3g %-12.3g %-10d\n', ...
        truncate_str(fn, 30), r.phi, r.lam_dem, r.T_mean, Tratio, r.cluster_flag);
end

%% ============== 3. Make figures =====================================
opt = struct('L', L, 'd', d, 'en', en, 'Cscale', C_SCALE, ...
             'savepath', SAVE_DIR, 'fit_C', FIT_C, ...
             'm', rho*(4/3)*pi*(d/2)^3, ...
             'T_ref', median([results.T_mean]));

[C_fit, phi_c_theory] = Fig1_LengthVsPhi(results, opt);
opt.Cscale = C_fit;                 % propagate the calibrated prefactor

Fig2_PhaseDiagram(results, opt);
Fig3_Prediction  (results, opt);

fprintf('\nCalibrated theory prefactor  C_scale = %.3g\n', C_fit);
if ~isnan(phi_c_theory)
    fprintf('Theory phi_c (lambda = L/2)  = %.4f\n', phi_c_theory);
end

fprintf('\nAll figures saved under %s\n', SAVE_DIR);

%% ====================================================================
function s = truncate_str(s, n)
    if length(s) > n, s = s(1:n); end
end

% ----------------------------------------------------------------------
% Synthetic demo data : used only if DATA_DIR is missing, so the user
% can quickly see what the plots look like.
function results = make_synthetic_results(L, d, rho, en)
% Synthetic demo data that mimics a realistic 3D DEM run with
% phi_c ~ 0.02 (matching the user's observed transition). We use a
% theoretical pre-factor C ~ 8 to convert default theory_lambda into the
% physical scale.
    m   = rho * (4/3) * pi * (d/2)^3;
    phis = [0.004 0.006 0.008 0.012 0.018 0.022 0.030 0.045 0.060];
    results = struct([]);
    T_wall  = 1.5e3;
    C_phys  = 8.0;            % calibrated to put phi_c at ~ 0.02
    rng(42);
    for k = 1:numel(phis)
        phi = phis(k);
        lam_th = theory_lambda(phi, en, d, 'Cscale', C_phys);
        lam_dem = lam_th * (1 + 0.10*(rand-0.5));
        Tc_over_Tw = 1 / cosh(L/(2*lam_dem));
        % Empirically (PPT)  phi_c ~ 0.02 in this experiment
        cflag = phi > 0.02;

        nBin = 25;
        re = linspace(0, L/2, nBin+1);
        rm = 0.5*(re(1:end-1)+re(2:end));
        Tprof = T_wall * cosh( (L/2 - rm)/lam_dem ) / cosh(L/(2*lam_dem)) ...
                .* (1 + 0.03*randn(1,nBin));
        n_avg = 6*phi / (pi*d^3);
        [tH, tM] = theory_tau(phi, en, d, T_wall, m);

        if cflag, s2 = 0.2 + 0.5*rand; else, s2 = 0.005 + 0.02*rand; end

        r = struct( ...
            'matfile', sprintf('synthetic_phi%.4f.mat', phi), ...
            'phi',       phi,   'L', L, 'd', d, 'm', m, 'en', en, ...
            'nFrames',   1500, ...
            'n_mean',    n_avg, 'T_mean', T_wall*0.7, ...
            'T_wall',    T_wall,'T_center', T_wall * Tc_over_Tw, ...
            'r_edges',   re,    'r_mid', rm, ...
            'T_of_rwall', Tprof,'n_of_rwall', n_avg*ones(1,nBin), ...
            'lam_dem',    lam_dem, ...
            'tau_dem_H',  tH,   'tau_dem_mfp', tM, ...
            'sigma2_n',   s2,   'cluster_flag', cflag );
        if isempty(results), results = r; else, results(end+1) = r; end %#ok<AGROW>
    end
end
