function [model,x0] = configure_vertical_instability_model(dt,equilibrium,em_coupling,pf_active,pf_passive)
%configure_vertical_instability_model Build a reduced vertical RZIp model.
%   [MODEL,X0] = configure_vertical_instability_model(DT) loads the machine
%   description and equilibrium using load_plasmaless_md and
%   extract_equilibrium_data. DT defaults to 0.001 s.
%
%   [MODEL,X0] = configure_vertical_instability_model(DT,EQ,EM,PA,PE) uses
%   supplied equilibrium and machine-description structures. EQ contains
%   RR, ZZ, Jtor [A/m^2], Ip0 [A], boundary, pf_active and pf_passive.
%   An optional Ic0 replaces the two current IDSs (active then passive).
%
%   States are perturbations [Ic;Ie] [A] of the 12 active circuits (CS1U/CS1L
%   in series, VS3U/VS3L in antiseries, see plasmaless_circuit_topology) and
%   of the passive loops. Inputs are the 14 power-supply voltages Vps [V]
%   (map_powersupplies_to_circuits) and outputs are [zIp;Ia;Ie] with the 14
%   coil currents Ia = T*Ic, zIp = Ip0*delta_z [A m]. Radial displacement and
%   plasma-current perturbations are excluded. Instantaneous vertical
%   force balance eliminates zIp using a Schur complement. X0 is zero.
%   MODEL.A/B/C/D and sysd are discrete, as in configure_plasmaless_model;
%   MODEL.Ac/Bc/Cc/Dc and sys are continuous. MODEL.M is the effective
%   inductance; M11/M12/M21/M22 expose the unreduced descriptor blocks.
%   Mxc = Mxi*blkdiag(T,I) is the circuit-to-grid flux map. An optional
%   EQ.Ic0 holds equilibrium circuit currents [Ic0;Ie0] (n_circuits+n_loops).
%
%   See also configure_plasmaless_model, extract_equilibrium_data

if nargin < 1 || isempty(dt), dt = 0.001; end
validateattributes(dt,{'numeric'},{'scalar','real','finite','positive'});
if nargin < 2 || isempty(equilibrium)
    equilibrium = extract_equilibrium_data();
end
if nargin < 3
    [em_coupling,pf_active,pf_passive] = load_plasmaless_md();
elseif nargin ~= 5
    error('vertical:inputs','Supply DT, EQ, EM, PA and PE together.');
end
base = configure_plasmaless_model(em_coupling,pf_active,pf_passive,dt);
[Mxi,RR,ZZ] = plasmaless_grid_from_ids(em_coupling);
% The state uses circuit currents Ic (CS1U/CS1L in series, VS3U/VS3L in
% antiseries), so the conductor-to-grid flux map is Mxc = Mxi*blkdiag(T,I)
% (Ia = T*Ic); radial-field derivatives are then taken for the circuits.
n_circuits = base.n_circuits;
Tf = blkdiag(base.T,eye(base.n_loops));
Mxc = Mxi*Tf;
n = n_circuits+base.n_loops;
dBrxdIc = zeros(size(Mxc));
for k = 1:n
    flux = reshape(Mxc(:,k),size(RR));
    [~,dfluxdz] = gradient(flux,RR(1,:),ZZ(:,1));
    dBrxdIc(:,k) = reshape(-dfluxdz./(2*pi*RR),[],1);
end
% Reuse RR/ZZ directly when the machine and equilibrium meshes agree.
% For a different equilibrium mesh, interpolate before reassigning RR/ZZ.
if ~isequal(RR,equilibrium.RR) || ~isequal(ZZ,equilibrium.ZZ)
    Mxc = interpolate_to_equilibrium(Mxc,RR,ZZ,equilibrium);
    dBrxdIc = interpolate_to_equilibrium(dBrxdIc,RR,ZZ,equilibrium);
end
RR = double(equilibrium.RR);
ZZ = double(equilibrium.ZZ);
J = double(equilibrium.Jtor);
Ip0 = double(equilibrium.Ip0);
validateattributes(Ip0,{'numeric'},{'scalar','finite','real','nonzero'});
if ~isequal(size(J),size(RR),size(ZZ)) || min(size(J)) < 3 ...
        || any(~isfinite([RR(:);ZZ(:);J(:)])) || any(RR(:)<=0)
    error('vertical:grid','RR, ZZ and Jtor must be finite matching grids.');
end
r = RR(1,:);
z = ZZ(:,1);
if any(diff(r)<=0) || any(diff(z)<=0) ...
        || max(abs(RR-r),[],'all') > 1e-10 ...
        || max(abs(ZZ-z),[],'all') > 1e-10
    error('vertical:grid','Use a tensor grid with Z rows and increasing R/Z.');
end
% Exclude vacuum current-density extrapolations outside the LCFS.
if isfield(equilibrium,'boundary') && ~isempty(equilibrium.boundary.r)
    inside = inpolygon(RR,ZZ,equilibrium.boundary.r,equilibrium.boundary.z);
    J(~inside) = 0;
end
area = integration_weights(z)*integration_weights(r)';
Igrid = J.*area;
Ip_grid = sum(Igrid,'all');
if abs(Ip_grid-Ip0) > 0.05*abs(Ip0)
    error('vertical:plasmaCurrent', ...
        'Grid current %.6g A differs from equilibrium Ip %.6g A by >5%%.', ...
        Ip_grid,Ip0);
end
% Correct the small quadrature discrepancy, preserving the IDS total Ip.
J = J*(Ip0/Ip_grid);
Igrid = J.*area;
% Negative spatial derivative = response to rigid vertical translation.
[~,dJdz] = gradient(J,r,z);
dIdz = -dJdz.*area;
if any(abs(J([1 end],:))>1e-8*max(abs(J(:))),'all') ...
        || any(abs(J(:,[1 end]))>1e-8*max(abs(J(:))),'all')
    error('vertical:plasmaBoundary','Plasma must fit inside the grid.');
end
% Only the plasma support and its gradient enter the integrals.
support = Igrid(:)~=0 | dIdz(:)~=0;
if any(~isfinite(Mxc(support,:)),'all') ...
        || any(~isfinite(dBrxdIc(support,:)),'all')
    error('vertical:coverage','Machine flux maps must cover the plasma.');
end
Mxc(~support,:) = 0;
dBrxdIc(~support,:) = 0;
if isfield(equilibrium,'Ic0')
    Ic0 = equilibrium.Ic0(:);   % circuit currents [Ic; Ie] (n_circuits+n_loops)
else
    active_equilibrium = equilibrium.pf_active;
    % DINA uses VSU/VSL for the same ITER circuits named VS3U/VS3L in
    % the DD4 machine description. Match this explicit naming difference.
    aliases = {'VSU','VS3U';'VSL','VS3L'};
    for k = 1:numel(active_equilibrium.coil)
        key = strtrim(active_equilibrium.coil{k}.name);
        token = regexp(key,'\(([^()]*)\)','tokens','once');
        if ~isempty(token), key = strtrim(token{1}); end
        for a = 1:size(aliases,1)
            if strcmp(key,aliases{a,1}) ...
                    && any(strcmp(base.coil_names,aliases{a,2}))
                active_equilibrium.coil{k}.name = aliases{a,2};
            end
        end
    end
    [~,Ic0] = configure_plasmaless_model(em_coupling,pf_active,pf_passive, ...
        dt,active_equilibrium,equilibrium.pf_passive);
end
validateattributes(Ic0,{'numeric'},{'column','numel',n,'finite','real'});
% Vacuum conductor field only, excluding the plasma self-field.
Br_external = reshape(dBrxdIc*Ic0,size(J));
M11 = base.M;
M12 = (dIdz(:)'*Mxc)'/Ip0;
% Direct Lorentz-force form used in the updated make_rzip_model.m.
M21 = -2*pi*(Igrid(:).*RR(:))'*dBrxdIc/Ip0;
M22 = -2*pi*sum(dIdz.*RR.*Br_external,'all')/Ip0^2;
if ~isfinite(M22) || abs(M22) <= eps*norm(M11,'fro')
    error('vertical:forceBalance','Vertical force stiffness is singular.');
end
Czr = -M21/M22;
Meff = M11-M12*(M21/M22);
% Inputs: power-supply voltages Vps, mapped to circuit voltages Vc = P*Vps
S = eye(n,n_circuits)*base.P;
Ac = -(Meff\base.R);
Bc = Meff\S;
% Outputs: zIp, the 14 coil currents Ia = T*Ic and the passive currents Ie
Cc = [Czr;Tf];
n_out = size(Tf,1);
Dc = zeros(n_out+1,base.n_ps);
sys = ss(Ac,Bc,Cc,Dc);
sys.InputName = strcat('Vps_',base.ps_names);
sys.OutputName = [{'zIp'},strcat('I_',base.coil_names), ...
    strcat('I_',base.loop_names)];
sys.StateName = [strcat('dIc_',base.circuit_names),strcat('dI_',base.loop_names)];
sysd = c2d(sys,dt,'zoh');
model = base;
model.sys = sys;
model.sysd = sysd;
model.A = sysd.A;
model.B = sysd.B;
model.C = sysd.C;
model.D = sysd.D;
model.Ac = Ac;
model.Bc = Bc;
model.Cc = Cc;
model.Dc = Dc;
model.M = Meff;
model.M11 = M11;
model.M12 = M12;
model.M21 = M21;
model.M22 = M22;
model.Czr = Czr;
model.Ip0 = Ip0;
model.Ic0 = Ic0;                 % equilibrium circuit currents [Ic0; Ie0]
model.Ia0 = Tf*Ic0;              % equilibrium coil currents [Ia0; Ie0] (single coils)
model.Mxc = Mxc;
model.dBrxdIc = dBrxdIc;
model.Jtor = J;
model.dJdzC = -dJdz;
model.RR = RR;
model.ZZ = ZZ;
model.Br_external = Br_external;
model.current_normalization = Ip0/Ip_grid;
model.ind_zIp = 1;
model.ind_Ia = 1+(1:base.n_coils);
model.ind_Ie = 1+base.n_coils+(1:base.n_loops);
model.singular_perturbation = true;
x0 = zeros(n,1);
end

function weights = integration_weights(nodes)
% Trapezoidal area weights for a possibly nonuniform tensor-product mesh.
nodes = nodes(:);
weights = [diff(nodes(1:2));nodes(3:end)-nodes(1:end-2); ...
    diff(nodes(end-1:end))]/2;
end

function interpolated = interpolate_to_equilibrium(values,RR,ZZ,equilibrium)
% Preserve flux and field interpolation for different equilibrium meshes.
interpolated = zeros(numel(equilibrium.RR),size(values,2));
for k = 1:size(values,2)
    interpolated(:,k) = reshape(interp2(RR,ZZ, ...
        reshape(values(:,k),size(RR)),equilibrium.RR,equilibrium.ZZ, ...
        'linear',NaN),[],1);
end
end
