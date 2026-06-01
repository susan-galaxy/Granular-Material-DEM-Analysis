function [tau_H, tau_mfp] = theory_tau(phi, en, d, T, m, varargin)
% Theoretical characteristic time scales for the granular gas.
%
%   tau_H   : Haff cooling time (energy decay e-folding time at wall T).
%             tau_H = 1 / [ gamma0 * (1-en^2) * n * d^2 * chi(phi) * v_T ]
%             with v_T = sqrt(T/m). This is the time the system would need
%             to cool by inelastic collisions if heat injection stopped.
%
%   tau_mfp : Enskog mean free time (between collisions per particle):
%             tau_mfp = 1 / [ sqrt(2) * pi * n * d^2 * chi(phi) * v_T ]
%             (3D hard sphere Enskog).
%
% Default constants match theory_lambda.m :
%   gamma0 = 12/sqrt(pi) , so factor (1-en^2) gives ~ same prefactor as in
%   the dissipation term used in the ODE derivation.
%
% Inputs:
%   phi : volume fraction (scalar/array)
%   en  : restitution coefficient
%   d   : particle diameter
%   T   : granular temperature (energy units = m * v^2 in cgs : g*cm^2/s^2)
%   m   : particle mass (g)
%   Name-Value :
%     'gamma0' : default 12/sqrt(pi)
%
% Outputs:
%   tau_H, tau_mfp (s)

    p = inputParser;
    p.addParameter('gamma0', 12/sqrt(pi));
    p.parse(varargin{:});

    chi = chi_CS(phi);
    n   = 6 * phi ./ (pi * d.^3);          % number density
    vT  = sqrt(T ./ m);                    % thermal velocity ~ sqrt(T/m)

    tau_H   = 1 ./ ( p.Results.gamma0 .* (1 - en.^2) .* n .* d.^2 .* chi .* vT );
    tau_mfp = 1 ./ ( sqrt(2) .* pi      .* n .* d.^2 .* chi .* vT );
end
