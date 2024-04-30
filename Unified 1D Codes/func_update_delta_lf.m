function [delta_lf] = func_update_delta_lf(delta_lf,k)
% This function updates delta_lf

if k<5
    delta_lf=min(1e-7,(10^(log10(delta_lf)+0.2)));
elseif k>=5 && k<=12
    delta_lf=delta_lf;
elseif k>12
    delta_lf=max(1e-24,(10^(log10(delta_lf)-0.2)));
end

end
