function equilibrium = extract_equilibrium_data(uri,time,opts)
%extract_equilibrium_data Load an equilibrium and its conductor currents.
%   EQ = extract_equilibrium_data() reads the ITER equilibrium at 100 s
%   from discharge 105108, run 2. EQ contains RR, ZZ, Ip0, Jtor, boundary,
%   pf_active and pf_passive, as used by the vertical model configurator.
%   psi0 and psix provide the signed full-flux map and LCFS flux for plots.
%   Arrays have Z along rows and R along columns.
%
%   EQ = extract_equilibrium_data(URI,TIME) selects another IMAS entry
%   and time. URI is any IMAS URI readable by the loaded IMAS-MATLAB
%   module (DD3 or DD4), or a plain path to an IMAS HDF5 entry directory
%   (or to its equilibrium.h5), converted to imas:hdf5?path=<dir>.
%   Interpolation uses the existing extractor's mode 3.
%
%   EQ = extract_equilibrium_data(URI,TIME,OPTS) with OPTS fields:
%     passive_currents  'equilibrium' (default): the entry's pf_passive
%                       must hold one current per machine-description loop
%                       (same names, same order), else vertical:passiveCurrents.
%                       'zero': passive currents set to 0 (explicit
%                       modelling assumption, recorded in passive_source).
%     interpolation     ids_get_slice mode (default 3, linear).
%     pf_passive_md     machine-description pf_passive for the loop names
%                       (default: load_plasmaless_md()).
%
%   psi sign: the machine-description flux maps need psi0 = -psi for a DD3
%   entry (COCOS 11) and psi0 = +psi for a DD4 entry (COCOS 17); the sign is
%   taken from the entry's DD version and checked against the signed Ip
%   with Ampere's law (vertical:fluxSign). EQ.psi_sign and EQ.ampere_ip
%   record the sign and the Ampere-law current.
%   LCFS: boundary.outline of the slice when it is a closed curve inside the
%   grid. An empty outline (linear interpolation drops it when neighbouring
%   samples have different sizes) or an open one (touching the grid edge,
%   as in the PDS scenario entries) is replaced by the closed contour of the
%   slice's psi at global_quantities.psi_boundary around the magnetic axis;
%   EQ.boundary_source says which was used.
%   Also returned: uri, time, dd_version, passive_source, boundary_source
%   and current_centre_entry ([r z] stored in the entry, NaN when absent).
%
%   See also load_plasmaless_md, configure_vertical_instability_model

if nargin < 1 || isempty(uri)
    uri = 'imas:mdsplus?path=/work/projects/dina/imasdb/iter/3/105108/2';
end
if nargin < 2 || isempty(time), time = 100; end
if nargin < 3, opts = struct(); end
opts = check_options(opts);
uri = entry_uri(uri);
idx = imas_open(uri);
cleanup = onCleanup(@() imas_close(idx));
slice = ids_get_slice(idx,'equilibrium',time,opts.interpolation);
if isempty(slice.time_slice) || isempty(slice.time_slice{1}.profiles_2d)
    error('vertical:equilibrium','%s: no equilibrium profiles_2d at t = %g.',uri,time);
end
p = slice.time_slice{1}.profiles_2d{1};
g = slice.time_slice{1}.global_quantities;
RR = double(p.r');
ZZ = double(p.z');
Ip0 = double(g.ip);
equilibrium.uri = uri;
equilibrium.time = time;
equilibrium.RR = RR;
equilibrium.ZZ = ZZ;
equilibrium.Ip0 = Ip0;
[equilibrium.boundary,equilibrium.boundary_source] = lcfs_outline(slice,uri);
% Match full poloidal flux in the machine-description conductor maps.
[psi_sign,equilibrium.dd_version] = flux_sign(slice,uri);
equilibrium.psi_sign = psi_sign;
equilibrium.psi0 = psi_sign*double(p.psi');
equilibrium.psix = psi_sign*double(slice.time_slice{1}.boundary.psi);
equilibrium.ampere_ip = ampere_current(RR,ZZ,equilibrium.psi0, ...
    equilibrium.boundary);
if sign(equilibrium.ampere_ip) ~= sign(Ip0)
    error('vertical:fluxSign', ['%s: psi0 = %+d*psi gives an Ampere-law ', ...
        'current %.6g A, the entry Ip is %.6g A.'],uri,psi_sign, ...
        equilibrium.ampere_ip,Ip0);
end
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
equilibrium.current_centre_entry = [entry_value(g.current_centre.r), ...
    entry_value(g.current_centre.z)];
equilibrium.pf_active = ids_get_slice(idx,'pf_active',time,opts.interpolation);
pf_passive = ids_get_slice(idx,'pf_passive',time,opts.interpolation);
md_names = cellfun(@(l) strtrim(l.name),opts.pf_passive_md.loop(:)', ...
    'UniformOutput',false);
reason = passive_currents_missing(pf_passive,md_names);
switch opts.passive_currents
    case 'equilibrium'
        if ~isempty(reason)
            error('vertical:passiveCurrents', ['%s: pf_passive %s. Pass ', ...
                'opts.passive_currents = ''zero'' to assume zero passive ', ...
                'currents.'],uri,reason);
        end
        equilibrium.pf_passive = pf_passive;
        equilibrium.passive_source = 'equilibrium entry';
    case 'zero'
        % Explicit modelling assumption: no vessel currents at the equilibrium.
        equilibrium.pf_passive = zero_passive(opts.pf_passive_md,time);
        if isempty(reason)
            equilibrium.passive_source = 'zero (requested, entry has currents)';
        else
            equilibrium.passive_source = 'zero (not in entry)';
        end
end
end

function opts = check_options(opts)
% Options with their defaults; unknown fields are errors.
if ~isstruct(opts) || ~isscalar(opts)
    error('vertical:options','OPTS must be a scalar struct.');
end
unknown = setdiff(fieldnames(opts)',{'passive_currents','interpolation','pf_passive_md'});
if ~isempty(unknown)
    error('vertical:options','Unknown option(s) {%s}.',strjoin(unknown,','));
end
if ~isfield(opts,'passive_currents'), opts.passive_currents = 'equilibrium'; end
opts.passive_currents = char(opts.passive_currents);
if ~any(strcmp(opts.passive_currents,{'equilibrium','zero'}))
    error('vertical:options', ...
        'passive_currents must be ''equilibrium'' or ''zero'' (got ''%s'').', ...
        opts.passive_currents);
end
if ~isfield(opts,'interpolation'), opts.interpolation = 3; end
if ~isfield(opts,'pf_passive_md')
    [~,~,opts.pf_passive_md] = load_plasmaless_md();
end
end

function uri = entry_uri(uri)
% IMAS URI as given, or imas:hdf5 for a plain entry directory path.
uri = char(uri);
if startsWith(uri,'imas:'), return; end
if isfile(uri)
    [folder,~,ext] = fileparts(uri);
    if ~strcmp(ext,'.h5')
        error('vertical:uri','%s is not an IMAS URI, entry directory or .h5 file.',uri);
    end
    uri = folder;
end
if ~isfolder(uri) || ~isfile(fullfile(uri,'equilibrium.h5'))
    error('vertical:uri','%s is not an IMAS URI or an entry directory with equilibrium.h5.',uri);
end
uri = ['imas:hdf5?path=' uri];
end

function [boundary,source] = lcfs_outline(slice,uri)
% Stored LCFS outline if closed inside the grid, else the closed psi contour
% at psi_boundary that encloses the magnetic axis.
ts = slice.time_slice{1};
p = ts.profiles_2d{1};
r = double(p.r(:,1)); z = double(p.z(1,:))';
boundary = ts.boundary.outline;
br = double(boundary.r(:)); bz = double(boundary.z(:));
if numel(br) > 2 && all(br > min(r) & br < max(r) & bz > min(z) & bz < max(z))
    source = 'boundary.outline';
    return
end
if isempty(br)
    reason = 'empty';
else
    reason = 'open (touches the grid edge)';
end
g = ts.global_quantities;
level = double(g.psi_boundary);
axis_rz = double([g.magnetic_axis.r g.magnetic_axis.z]);
if numel(level) ~= 1 || ~isfinite(level) || level == -9e40 ...
        || numel(axis_rz) ~= 2 || any(~isfinite(axis_rz)) || any(axis_rz == -9e40)
    error('vertical:boundary', ['%s: boundary.outline is %s and psi_boundary ', ...
        'or the magnetic axis is missing.'],uri,reason);
end
C = contourc(r,z,double(p.psi'),[level level]);
boundary = [];
k = 1;
while k < size(C,2)
    n = C(2,k);
    seg = C(:,k+1:k+n);
    if isequal(seg(:,1),seg(:,end)) && inpolygon(axis_rz(1),axis_rz(2),seg(1,:),seg(2,:))
        boundary = struct('r',seg(1,:)','z',seg(2,:)');
        break
    end
    k = k+n+1;
end
if isempty(boundary)
    error('vertical:boundary', ['%s: boundary.outline is %s and no closed psi ', ...
        'contour at psi_boundary = %g encloses the magnetic axis.'],uri,reason,level);
end
source = sprintf('psi contour at psi_boundary = %.6g Wb (boundary.outline %s)',level,reason);
end

function [psi_sign,dd] = flux_sign(slice,uri)
% DD3 entries are COCOS 11, DD4 entries COCOS 17: psi changes sign between
% them; the machine-description maps need -psi of a DD3 entry.
dd = strtrim(char(slice.ids_properties.version_put.data_dictionary));
major = sscanf(dd,'%d',1);
if isequal(major,3)
    psi_sign = -1;
elseif isequal(major,4)
    psi_sign = 1;
else
    error('vertical:ddVersion','%s: unknown data-dictionary version ''%s''.',uri,dd);
end
end

function Ip_flux = ampere_current(RR,ZZ,psi0,boundary)
% Ampere's law inside the LCFS, as the check in plot_vertical_instability.
r = RR(1,:); z = ZZ(:,1);
[dpsiR,dpsiZ] = gradient(psi0,r,z);
[dpsiRR,~] = gradient(dpsiR,r,z);
[~,dpsiZZ] = gradient(dpsiZ,r,z);
J = -(dpsiRR-dpsiR./RR+dpsiZZ)./(2*pi*4e-7*pi*RR);
J(~inpolygon(RR,ZZ,boundary.r,boundary.z)) = 0;
Ip_flux = trapz(r,trapz(z,J,1));
end

function v = entry_value(v)
% Stored scalar, NaN when missing.
if ~isnumeric(v) || numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
    v = NaN;
end
v = double(v);
end

function reason = passive_currents_missing(pf_passive,md_names)
% '' when pf_passive has one current per machine-description loop.
reason = '';
if numel(pf_passive.loop) ~= numel(md_names)
    reason = sprintf('has %d loops, the machine description %d', ...
        numel(pf_passive.loop),numel(md_names));
    return
end
for k = 1:numel(md_names)
    name = strtrim(pf_passive.loop{k}.name);
    c = pf_passive.loop{k}.current;
    if ~strcmp(name,md_names{k})
        reason = sprintf('loop %d is ''%s'', the machine description ''%s''', ...
            k,name,md_names{k});
        return
    end
    if numel(c) ~= 1 || c == -9e40 || ~isfinite(c)
        reason = sprintf('loop %s has no current',name);
        return
    end
end
end

function pf_passive = zero_passive(pf_passive_md,time)
% One-slice pf_passive with the machine-description loops and zero currents.
pf_passive = ids_init('pf_passive');
pf_passive.ids_properties.homogeneous_time = 1;
pf_passive.time = time;
pf_passive.loop = ids_allocate('pf_passive','loop',numel(pf_passive_md.loop));
for k = 1:numel(pf_passive_md.loop)
    pf_passive.loop{k}.name = pf_passive_md.loop{k}.name;
    pf_passive.loop{k}.current = 0;
end
end
