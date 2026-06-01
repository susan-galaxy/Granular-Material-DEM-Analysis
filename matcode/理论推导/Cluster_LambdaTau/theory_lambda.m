function lam = theory_lambda(phi, en, d, varargin)
% Theoretical characteristic length (energy penetration depth) from the
% ballistic-thermalized two-fluid model.
%
% Derivation (see PPT / Word notes):
%   kappa * d^2 T / dz^2 ~ Gamma                  (bulk balance, no source)
%   kappa = kappa0 * n * d * sqrt(T)
%   Gamma = gamma0 * n^2 * d^2 * chi(phi) * (1-en^2) * T^(3/2)
%   =>  d^2 T/dz^2 = T / lambda^2
%   with  lambda^2 = kappa0 / [ gamma0 * (1-en^2) * n * d * chi(phi) ]
%   Using n*d^3 = 6 phi / pi:
%     lambda = C0 * d / sqrt( (1-en^2) * phi * chi(phi) )
%     C0 = sqrt( kappa0 * pi / (6 * gamma0) )
%
% Default hard-sphere constants (Chapman-Enskog) :
%   kappa0 = sqrt(pi)/6 ,  gamma0 = 12/sqrt(pi)
%   => C0 = pi / (12 * sqrt(3)) ~ 0.151
%
% A multiplicative tunable prefactor 'Cscale' may be passed to absorb
% Sonine corrections / 3D geometric factors. Default 1.
%
% Inputs:
%   phi    : volume fraction (scalar/array)
%   en     : restitution coefficient
%   d      : particle diameter (cm in our convention)
%   Name-Value:
%     'Cscale' : multiplicative O(1) prefactor (default 1)
%     'kappa0' : default sqrt(pi)/6
%     'gamma0' : default 12/sqrt(pi)
%
% Output:
%   lam : characteristic length (same units as d)

    p = inputParser;
    p.addParameter('Cscale', 1.0);
    p.addParameter('kappa0', sqrt(pi)/6);
    p.addParameter('gamma0', 12/sqrt(pi));
    p.parse(varargin{:});

    C0 = sqrt( p.Results.kappa0 * pi / (6 * p.Results.gamma0) );
    C  = C0 * p.Results.Cscale;

    chi = chi_CS(phi);
    lam = C .* d ./ sqrt( max((1 - en.^2) .* phi .* chi, eps) );
end
