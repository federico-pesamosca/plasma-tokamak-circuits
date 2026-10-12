function [em_coupling, pf_active, pf_passive, wall] = load_plasmaless_md(md_uri)
% read the DD4 machine description (em_coupling, pf_active, pf_passive, wall)
% default: the .h5 entry in data/md_dd4 of this repository
% wall (optional 4th output) holds the limiter outlines used for plotting

if nargin < 1 || isempty(md_uri)
    repo_root = fileparts(fileparts(mfilename('fullpath')));   % this file is in imas_model/
    md_uri = ['imas:hdf5?path=' fullfile(repo_root,'data','md_dd4')];
end
idx = imas_open(md_uri,40);
em_coupling = ids_get(idx,'em_coupling');
pf_active = ids_get(idx,'pf_active');
pf_passive = ids_get(idx,'pf_passive');
if nargout > 3
    wall = ids_get(idx,'wall');
end
imas_close(idx);
end
