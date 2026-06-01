function chi = chi_CS(phi)
% Carnahan-Starling pair correlation function at contact for hard spheres.
% chi(phi) = (1 - phi/2) / (1 - phi)^3
%
% Inputs:
%   phi : volume fraction (scalar or array)
% Outputs:
%   chi : Enskog correction factor at contact

    chi = (1 - phi/2) ./ (1 - phi).^3;
end
