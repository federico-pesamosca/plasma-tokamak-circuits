function equilibrium = extract_equilibrium_data(uri,time)
%extract_equilibrium_data Load an equilibrium and its conductor currents.
%   EQ = extract_equilibrium_data() reads the ITER equilibrium at 100 s
%   from discharge 105108, run 2. EQ contains RR, ZZ, Ip0, Jtor, boundary,
%   pf_active and pf_passive, as used by the vertical model configurator.
%   psi0 and psix provide the signed full-flux map and LCFS flux for plots.
%   Arrays have Z along rows and R along columns.
%
%   EQ = extract_equilibrium_data(URI,TIME) selects another IMAS entry
%   and time. Interpolation uses the existing extractor's mode 3.
%
%   See also load_plasmaless_md, configure_vertical_instability_model

if nargin < 1 || isempty(uri)
    uri = 'imas:mdsplus?path=/work/projects/dina/imasdb/iter/3/105108/2';
end
if nargin < 2 || isempty(time), time = 100; end
idx = imas_open(uri);
cleanup = onCleanup(@() imas_close(idx));
slice = ids_get_slice(idx,'equilibrium',time,3);
p = slice.time_slice{1}.profiles_2d{1};
g = slice.time_slice{1}.global_quantities;
RR = double(p.r');
ZZ = double(p.z');
Ip0 = double(g.ip);
equilibrium.RR = RR;
equilibrium.ZZ = ZZ;
equilibrium.Ip0 = Ip0;
equilibrium.boundary = slice.time_slice{1}.boundary.outline;
% Match full poloidal flux in the machine-description conductor maps.
equilibrium.psi0 = -double(p.psi');
equilibrium.psix = -double(slice.time_slice{1}.boundary.psi);
if ~isempty(p.j_phi)
    equilibrium.Jtor = double(p.j_phi');
else
    % This DINA entry stores total flux [Wb], p'(psi) and F F'(psi),
    % but not 2D current density. Grad-Shafranov gives its distribution.
    p1 = slice.time_slice{1}.profiles_1d;
    psi = double(p.psi');
    dp = interp1(double(p1.psi),double(p1.dpressure_dpsi),psi,'linear',0);
    ff = interp1(double(p1.psi),double(p1.f_df_dpsi),psi,'linear',0);
    J = 2*pi*(RR.*dp+ff./(4e-7*pi*RR));
    inside = inpolygon(RR,ZZ,equilibrium.boundary.r,equilibrium.boundary.z);
    J(~inside) = 0;
    integrated_current = trapz(RR(1,:),trapz(ZZ(:,1),J,1));
    if ~isfinite(integrated_current) || integrated_current == 0 ...
            || abs(abs(integrated_current)-abs(Ip0)) > 0.05*abs(Ip0)
        error('vertical:equilibriumCurrent', ...
            'Grad-Shafranov current disagrees with Ip by more than 5%%.');
    end
    % Infer the source flux orientation from the signed IDS plasma current.
    % Do not infer a force-stiffness sign or impose an unstable eigenvalue.
    orientation = sign(Ip0/integrated_current);
    equilibrium.Jtor = orientation*J;
end
equilibrium.pf_active = ids_get_slice(idx,'pf_active',time,3);
equilibrium.pf_passive = ids_get_slice(idx,'pf_passive',time,3);
end
