function topo = plasmaless_circuit_topology(coil_names)
%plasmaless_circuit_topology Coil <-> circuit <-> power-supply connections.
%   TOPO = plasmaless_circuit_topology() returns the ITER topology of the
%   14 active coils (order of pf_active in data/md_dd4):
%     coils    {CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 PF6 VS3U VS3L}
%     circuits {CS3U CS2U CS1 CS2L CS3L PF1 PF2 PF3 PF4 PF5 PF6 VS3}   (12)
%     supplies {CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3} (14)
%   - CS1U and CS1L are in SERIES:      Ia_CS1U = Ia_CS1L = Ic_CS1
%   - VS3U and VS3L are in ANTISERIES:  Ia_VS3U = -Ic_VS3, Ia_VS3L = +Ic_VS3
%   so a positive Ic_VS3 (and positive Vc_VS3) gives a positive VS3L current.
%
%   TOPO.T  (n_coils x n_circuits):  Ia = T*Ic
%   TOPO.P  (n_circuits x n_ps):     Vc = P*Vps (see map_powersupplies_to_circuits)
%   TOPO.Tpinv = (T'*T)\T':          Ic = Tpinv*Ia (least-squares projection)
%   With these, the coil quantities M, R of the n_coils coils become
%   Mcc = T'*Maa*T, Rcc = T'*Raa*T and the circuit voltages Vc = T'*Va_coils.
%
%   TOPO = plasmaless_circuit_topology(COIL_NAMES) uses the ITER topology if
%   COIL_NAMES equals the 14 names above (same order); any other coil set
%   (for example the one-coil fixtures of the tests) is treated with one
%   circuit and one power supply per coil (T = P = identity).
%
%   See also map_powersupplies_to_circuits, map_coils_to_circuits,
%   configure_plasmaless_model

iter_coils = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3', ...
    'PF4','PF5','PF6','VS3U','VS3L'};
if nargin < 1 || isempty(coil_names), coil_names = iter_coils; end
coil_names = cellstr(coil_names(:))';

if isequal(coil_names,iter_coils)
    circuit_names = {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF2','PF3', ...
        'PF4','PF5','PF6','VS3'};
    ps_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3', ...
        'PF4','PF5','VS1','PF6','VS3'};
    % Ia = T*Ic
    T = zeros(14,12);
    T(1,1) = 1;  T(2,2) = 1;
    T(3,3) = 1;  T(4,3) = 1;                  % CS1U, CS1L: series
    T(5,4) = 1;  T(6,5) = 1;
    T(7:12,6:11) = eye(6);                    % PF1..PF6
    T(13,12) = -1;  T(14,12) = 1;             % VS3U, VS3L: antiseries
    % Vc = P*Vps
    P = zeros(12,14);
    P(1,1) = 1;  P(2,2) = 1;
    P(3,[3 4]) = 1;                           % Vc_CS1 = Vps_CS1U + Vps_CS1L (series supplies)
    P(4,5) = 1;  P(5,6) = 1;
    P(6,7) = 1;                               % PF1
    P(7,[8 12]) = 1;                          % PF2: + Vps_VS1
    P(8,[9 12]) = 1;                          % PF3: + Vps_VS1
    P(9,10) = 1;   P(9,12) = -1;              % PF4: - Vps_VS1
    P(10,11) = 1;  P(10,12) = -1;             % PF5: - Vps_VS1
    P(11,13) = 1;                             % PF6
    P(12,14) = 1;                             % Vc_VS3 = Vps_VS3
    is_iter = true;
else
    n = numel(coil_names);
    circuit_names = coil_names;
    ps_names = coil_names;
    T = eye(n);
    P = eye(n);
    is_iter = false;
end

topo.coil_names = coil_names;
topo.circuit_names = circuit_names;
topo.ps_names = ps_names;
topo.T = T;
topo.P = P;
topo.Tpinv = (T'*T)\T';
topo.n_coils = size(T,1);
topo.n_circuits = size(T,2);
topo.n_ps = size(P,2);
topo.is_iter = is_iter;
end
