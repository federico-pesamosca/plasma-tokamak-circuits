function [I_c,topo] = map_coils_to_circuits(Ia,topo)
%map_coils_to_circuits Coil currents -> circuit currents.
%   I_c = map_coils_to_circuits(IA) maps the 14 coil currents
%     Ia = {CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 PF6 VS3U VS3L}
%   to the 12 circuit currents
%     Ic = {CS3U CS2U CS1 CS2L CS3L PF1 PF2 PF3 PF4 PF5 PF6 VS3}:
%     Ic_CS1 = Ia_CS1U = Ia_CS1L   (series),
%     Ic_VS3 = Ia_VS3L = -Ia_VS3U  (antiseries), others one to one.
%   IA is a 14 element vector (I_c is a column) or an N x 14 matrix with one
%   sample per row (I_c is N x 12). Coil currents that violate the series /
%   antiseries constraints (e.g. Ia_CS1U ~= Ia_CS1L) are not representable
%   by circuit currents: an error is raised if the mismatch exceeds 1e-6 of
%   the largest coil current (a warm start must be consistent).
%   [I_c,TOPO] = map_coils_to_circuits(IA,TOPO) uses another topology from
%   plasmaless_circuit_topology (default: ITER).
%
%   See also plasmaless_circuit_topology, map_powersupplies_to_circuits

if nargin < 2, topo = plasmaless_circuit_topology(); end
n = topo.n_coils;
if isvector(Ia) && numel(Ia) == n
    Ia = Ia(:)';
    single = true;
elseif size(Ia,2) == n
    single = false;
else
    error('plasmaless:coilCurrents','Ia must have %d elements or %d columns.',n,n);
end
I_c = Ia*topo.Tpinv.';
residual = Ia - I_c*topo.T.';
scale = max(1,max(abs(Ia(:))));
if max(abs(residual(:))) > 1e-6*scale
    error('plasmaless:seriesConstraint', ...
        ['Coil currents are not consistent with the circuit topology ', ...
        '(series CS1U=CS1L, antiseries VS3U=-VS3L): max mismatch %.3g A.'], ...
        max(abs(residual(:))));
end
if single, I_c = I_c(:); end
end
