
function fk_case = classify_fk_case(Tmix, Tin, dT_fully_mixed)
% ---------------------------------------------------------------
% SUMMARY:
%   Classifies the numerical free-convection regime.
%
% INPUT:
%   Tmix            - temperature vector in mixing region
%   Tin             - inlet temperature
%   dT_fully_mixed  - threshold for fully mixed assumption [K]
%
% OUTPUT:
%   fk_case         - FK regime identifier:
%                     'internal_zmix' | 'extrapolated_zmix' | 'fully_mixed'
% ---------------------------------------------------------------

    Tmin = min(Tmix);
    Tmax = max(Tmix);

    % Fully mixed: almost uniform temperature profile
    if (Tmax - Tmin) <= dT_fully_mixed
        fk_case = 'fully_mixed';

    % Inlet temperature lies within temperature range
    elseif  Tin >= min(Tmix) && Tin <= max(Tmix)
        fk_case = 'internal_zmix';

    % Inlet temperature outside temperature range
    else
        fk_case = 'extrapolated_zmix';
    end
end