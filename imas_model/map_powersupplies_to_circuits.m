function [V_c,P] = map_powersupplies_to_circuits(Vps,topo)
%map_powersupplies_to_circuits Power-supply voltages -> circuit voltages.
%   V_c = map_powersupplies_to_circuits(VPS) maps the 14 power-supply voltages
%     Vps = {CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3}
%   to the 12 circuit voltages
%     Vc  = {CS3U CS2U CS1 CS2L CS3L PF1 PF2 PF3 PF4 PF5 PF6 VS3}:
%     Vc_CS1 = Vps_CS1U + Vps_CS1L           (supplies in series)
%     Vc_PF2 = Vps_PF2 - Vps_VS1, Vc_PF3 = Vps_PF3 - Vps_VS1
%     Vc_PF4 = Vps_PF4 + Vps_VS1, Vc_PF5 = Vps_PF5 + Vps_VS1
%     Vc_VS3 = Vps_VS3, all others one to one.
%   VPS is a vector with 14 elements (V_c has the orientation of VPS) or an
%   N x 14 matrix with one sample per row, as in timeseries data (V_c is
%   N x 12). [V_c,P] = ... also returns the matrix P with Vc = P*Vps.
%   map_powersupplies_to_circuits(VPS,TOPO) uses another topology from
%   plasmaless_circuit_topology.
%
%   See also plasmaless_circuit_topology, map_coils_to_circuits

if nargin < 2, topo = plasmaless_circuit_topology(); end
P = topo.P;
n = topo.n_ps;
if iscolumn(Vps) && numel(Vps) == n
    V_c = P*Vps;
elseif size(Vps,2) == n
    V_c = Vps*P.';
else
    error('plasmaless:powerSupplies', ...
        'Vps must have %d elements or %d columns (got %s).',n,n,mat2str(size(Vps)));
end
end
