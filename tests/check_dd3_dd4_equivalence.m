% Data-equivalence check between the DD3 sources of the plasmaless machine
% description and the DD4 entry data/md_dd4 (made from them with
% tools/convert_md_dd3_to_dd4.py and read by imas_model/load_plasmaless_md.m):
% the conductor mutual inductances M = [Maa Mea'; Mea Mee], the coil-to-grid
% flux maps Mxi = [Mxa Mxe] with their (R,Z) mesh, the coil and loop names and
% resistances, and the model built by configure_plasmaless_model. Exact equality
% is expected (the converter copies the data). Two MATLAB sessions (one module
% stack each), run from the repository root:
%
% part 1 (DD3 sources: DINA em_coupling/pf_active and the equilibrium grid,
%         ITER_MD pf_passive; saves ref_file). The flux maps mutual_grid_active/
%         passive were removed from the DD in 3.40.0, so a DD <= 3.39 module is needed:
%   module purge; module load IMAS/3.39.0-foss-2023b
%   module unload MATLAB/2023b-r5-GCCcore-13.2.0; module load MATLAB/2024a-r6-GCCcore-13.2.0
%   matlab -singleCompThread -nodisplay -batch "addpath('tests'); part=1; check_dd3_dd4_equivalence"
% part 2 (DD4 entry data/md_dd4, compared with ref_file):
%   module purge; module load IMAS-MATLAB/5.5.0-intel-2023b-DD-4.1.0
%   module unload MATLAB/2023b-r5-GCCcore-13.2.0; module load MATLAB/2024a-r6-GCCcore-13.2.0
%   matlab -singleCompThread -nodisplay -batch "addpath('tests'); part=2; check_dd3_dd4_equivalence"
% (called as a script on the path, not with run(): see the note below)

if ~exist('part','var')
    error('Set part = 1 (DD <= 3.39 module) or part = 2 (DD4 module) before running this script.');
end
% repository root from the current folder (tests/ if the script was started
% from there). No mfilename/which/run() here: with MATLAB R2024a (IMAS-MATLAB
% 5.5.0 modules) these file lookups, like load(), bring in MATLAB's own HDF5
% 1.8 library, and the IMAS access layer (HDF5 1.14) then fails with a fatal
% version mismatch if its entry has not been read yet; part 2 therefore reads
% the DD4 entry before load(ref_file).
rootDir = pwd;
if ~isfolder(fullfile(rootDir,'imas_model'))
    rootDir = fileparts(rootDir);
end
if ~isfolder(fullfile(rootDir,'imas_model'))
    error('Run this script from the repository root.');
end
addpath(rootDir, fullfile(rootDir,'imas_model'));

dt = 0.001;
tol = 1e-15;   % data copied, not recomputed: exact equality expected
ref_file = fullfile('/scratch/users',getenv('USER'),'plasmaless','dd3_reference.mat');
dina_uri = 'imas:mdsplus?path=/work/projects/dina/SRO_JINTRAC/15MA_10perc/output';
md_passive_uri = 'imas:mdsplus?path=/work/imas/shared/imasdb/ITER_MD/3/115005/3';
reldiff = @(X,Y) max(abs(X-Y),[],'all')/max(abs(Y),[],'all');

if part == 1
    % DD3 em_coupling, pf_active and the equilibrium grid of the DINA entry
    idx = imas_open(dina_uri);
    em = ids_get(idx,'em_coupling');
    pf_active = ids_get(idx,'pf_active');
    if ~isfield(em,'mutual_grid_active')
        error(['em_coupling has no mutual_grid_active: part 1 needs a DD <= 3.39 module ', ...
            '(e.g. IMAS/3.39.0-foss-2023b).']);
    end
    Maa = em.mutual_active_active;
    Mea = em.mutual_passive_active;
    Mee = em.mutual_passive_passive;
    M3 = [Maa Mea'; Mea Mee];
    Mxi3 = [em.mutual_grid_active, em.mutual_grid_passive];
    % mesh of the flux maps: equilibrium profiles_2d{1} grid (65 R x 129 Z, Z
    % varying fastest in the map rows), transposed to nZ x nR as in plot_plasmaless_model.m
    s = ids_get_slice(idx,'equilibrium',1,1);
    RR3 = s.time_slice{1}.profiles_2d{1}.r';
    ZZ3 = s.time_slice{1}.profiles_2d{1}.z';
    imas_close(idx);
    idx = imas_open(md_passive_uri);
    pf_passive = ids_get(idx,'pf_passive');
    imas_close(idx);

    n_coils = numel(pf_active.coil);
    n_loops = numel(pf_passive.loop);
    coil_identifiers3 = cell(1,n_coils); coil_names3_raw = cell(1,n_coils);
    coil_names3 = cell(1,n_coils); Rcoil3 = zeros(1,n_coils);
    for ii = 1:n_coils
        coil_identifiers3{ii} = strtrim(pf_active.coil{ii}.identifier);
        coil_names3_raw{ii} = strtrim(pf_active.coil{ii}.name);
        % short key as in the model: identifier when present, else name
        if ~isempty(coil_identifiers3{ii})
            coil_names3{ii} = coil_identifiers3{ii};
        else
            coil_names3{ii} = coil_names3_raw{ii};
        end
        Rcoil3(ii) = pf_active.coil{ii}.resistance;
    end
    loop_names3 = cell(1,n_loops); Rloop3 = zeros(1,n_loops);
    for ii = 1:n_loops
        loop_names3{ii} = strtrim(pf_passive.loop{ii}.name);
        Rloop3(ii) = pf_passive.loop{ii}.resistance;
    end
    R3 = diag([Rcoil3 Rloop3]);

    fprintf('DD3 sources: %d coils, %d loops; Maa %s, Mea %s, Mee %s, M %s\n', n_coils, n_loops, ...
        mat2str(size(Maa)), mat2str(size(Mea)), mat2str(size(Mee)), mat2str(size(M3)));
    fprintf('DD3 grid maps: Mxi %s, mesh %s (R %.7g..%.7g, Z %.7g..%.7g)\n', ...
        mat2str(size(Mxi3)), mat2str(size(RR3)), min(RR3,[],'all'), max(RR3,[],'all'), ...
        min(ZZ3,[],'all'), max(ZZ3,[],'all'));
    fprintf('DD3 coil keys: %s\n', strjoin(coil_names3,' '));
    assert(isequal(size(M3),[n_coils+n_loops n_coils+n_loops]), 'DD3 mutual matrices and coil/loop counts differ.');
    assert(size(Mxi3,1) == numel(RR3), 'DD3 grid maps and equilibrium mesh sizes differ.');
    assert(size(Mxi3,2) == n_coils+n_loops, 'DD3 grid maps and coil/loop counts differ.');

    if ~exist(fileparts(ref_file),'dir')
        mkdir(fileparts(ref_file));
    end
    save(ref_file,'M3','Maa','Mea','Mee','Mxi3','RR3','ZZ3','R3', ...
        'coil_names3','coil_identifiers3','coil_names3_raw','loop_names3','dina_uri','md_passive_uri','-v7');
    fprintf('Saved %s\n', ref_file);

elseif part == 2
    [em_coupling, pf_active, pf_passive] = load_plasmaless_md();   % data/md_dd4, before load()
    ref = load(ref_file);
    Maa = get_coupling_matrix(em_coupling,'mutual_active_active','pf_active/coil','pf_active/coil');
    Mea = get_coupling_matrix(em_coupling,'mutual_passive_active','pf_passive/loop','pf_active/coil');
    Mee = get_coupling_matrix(em_coupling,'mutual_passive_passive','pf_passive/loop','pf_passive/loop');
    M = [Maa Mea'; Mea Mee];
    tic;
    [Mxi, RR, ZZ] = plasmaless_grid_from_ids(em_coupling);
    fprintf('plasmaless_grid_from_ids: %.2f s, Mxi %s, mesh %s\n', toc, mat2str(size(Mxi)), mat2str(size(RR)));
    model = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt);

    assert(isequal(size(Maa),size(ref.Maa)) && isequal(size(Mea),size(ref.Mea)) ...
        && isequal(size(Mee),size(ref.Mee)), 'DD4 mutual matrix sizes differ from DD3.');
    assert(isequal(size(Mxi),size(ref.Mxi3)) && isequal(size(RR),size(ref.RR3)) ...
        && isequal(size(ZZ),size(ref.ZZ3)), 'DD4 grid maps or mesh sizes differ from DD3.');
    fprintf('DD4 vs DD3 mutuals: Maa %.3g, Mea %.3g, Mee %.3g, M %.3g (max relative difference), M exact: %d\n', ...
        reldiff(Maa,ref.Maa), reldiff(Mea,ref.Mea), reldiff(Mee,ref.Mee), reldiff(M,ref.M3), isequal(M,ref.M3));
    fprintf('DD4 vs DD3 grid: Mxi %.3g, RR %.3g, ZZ %.3g (max relative difference), Mxi exact: %d, mesh exact: %d\n', ...
        reldiff(Mxi,ref.Mxi3), reldiff(RR,ref.RR3), reldiff(ZZ,ref.ZZ3), ...
        isequal(Mxi,ref.Mxi3), isequal(RR,ref.RR3) && isequal(ZZ,ref.ZZ3));
    fprintf('DD4 model vs DD3: model.M %.3g, model.R %.3g (max relative difference)\n', ...
        reldiff(model.M,ref.M3), reldiff(model.R,ref.R3));
    fprintf('DD4 coil keys: %s\n', strjoin(model.coil_names,' '));
    fprintf('Coil names equal: %d, loop names equal: %d (%d coils, %d loops)\n', ...
        isequal(model.coil_names,ref.coil_names3), isequal(model.loop_names,ref.loop_names3), ...
        model.n_coils, model.n_loops);

    assert(isequal(model.coil_names,ref.coil_names3), 'Coil names/order differ between DD3 and DD4.');
    assert(isequal(model.loop_names,ref.loop_names3), 'Loop names/order differ between DD3 and DD4.');
    assert(reldiff(M,ref.M3) <= tol, 'DD4 mutual inductances differ from DD3.');
    assert(reldiff(Mxi,ref.Mxi3) <= tol && reldiff(RR,ref.RR3) <= tol && reldiff(ZZ,ref.ZZ3) <= tol, ...
        'DD4 grid maps or mesh differ from DD3.');
    assert(reldiff(model.M,ref.M3) <= tol, 'configure_plasmaless_model M differs from the DD3 matrices.');
    assert(reldiff(model.R,ref.R3) <= tol, 'configure_plasmaless_model R differs from the DD3 resistances.');
    fprintf('DD3/DD4 equivalence check passed.\n');
else
    error('part must be 1 or 2.');
end
