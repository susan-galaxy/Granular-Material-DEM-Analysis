function phi = phi_from_name(name)
% Parse the volume fraction phi from a folder / file name that contains
% something like 'phi0.01', 'phi_0.012', 'phi0p02', '_phi-0.04_' etc.
%
% Returns NaN if no match is found.

    name = char(name);
    tok = regexpi(name, 'phi[_\-=]?([0-9]+(?:[\.p][0-9]+)?)', 'tokens', 'once');
    if isempty(tok)
        phi = NaN;
        return;
    end
    s = strrep(tok{1}, 'p', '.');
    phi = str2double(s);
end
