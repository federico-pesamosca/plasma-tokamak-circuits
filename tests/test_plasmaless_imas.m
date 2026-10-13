classdef test_plasmaless_imas < matlab.unittest.TestCase
    % Unit tests for the IMAS interface of the plasmaless model (imas_model/).
    %
    % The machine description is the DD4 entry data/md_dd4 of this repository
    % (load_plasmaless_md default), or the IMAS URI in the environment variable
    % PLASMALESS_MD_URI when set. Run with an IMAS-MATLAB DD-4.x module from the
    % repository root: runtests('tests/test_plasmaless_imas.m'). With MATLAB
    % R2024a (IMAS-MATLAB 5.5.0 modules) read an IMAS entry first, e.g.
    % add_plasmaless_paths; load_plasmaless_md(); runtests(...): the test
    % framework otherwise loads MATLAB's own HDF5 1.8 before the access
    % layer's HDF5 1.14 and MATLAB exits with a fatal HDF5 version error.
    % The Simulink comparison needs a Simulink licence; sim() resolves A, B, C, D, x0, TS,
    % dt, tstart and tend from the base workspace (as run_plasmaless_model.m),
    % so they are pushed there with assignin(). The equivalence of the DD3
    % sources with the DD4 entry is checked by tests/check_dd3_dd4_equivalence.m.
    % Three tests
    % 1) test_simulink_equivalence (IMAS step loop vs simulator_plasmaless_model_simulink.slx)
    % 2) test_contract_errors (invalid inputs raise errors)
    % 3) test_resistive_coil (expected step response for VS3U/VS3L, as in test_plasmaless_model)
    properties
        modelName = 'simulator_plasmaless_model_simulink';
        rootDir
        em_coupling
        pf_active
        pf_passive
    end

    methods (TestClassSetup)
        function setup(testCase)
            % tests/ is one level below the repository root; the IMAS functions
            % are in imas_model/ and plasmaless_step_ids calls
            % plasmaless_timestep_forward.m in the root
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(testCase.rootDir, fullfile(testCase.rootDir,'imas_model'));
            md_uri = getenv('PLASMALESS_MD_URI');   % empty: data/md_dd4 default
            [testCase.em_coupling, testCase.pf_active, testCase.pf_passive] = ...
                load_plasmaless_md(md_uri);
        end
    end

    methods (TestMethodSetup)
        function freshBaseWorkspace(testCase)
            evalin('base','clear');
            testCase.addTeardown(@() evalin('base','clear'));
        end
    end

    methods (Test)

        function test_simulink_equivalence(testCase)
            % Same A,B,C,D and input in the Simulink model and in the IMAS step
            % loop. simout at sample k is the state before the update with
            % u_k, i.e. the IMAS output of step k-1.
            tstart = 0; tend = 28; dt = 0.01;
            model = testCase.buildModel(dt);
            TS = timeseries_plasmaless_model(tstart,tend,dt);
            x0 = zeros(model.n_coils+model.n_loops,1);

            % variables read by simulator_plasmaless_model_simulink.slx (and
            % the referenced plasmaless_timestep.slx) from the base workspace
            assignin('base','A',model.A);
            assignin('base','B',model.B);
            assignin('base','C',model.C);
            assignin('base','D',model.D);
            assignin('base','x0',x0);
            assignin('base','TS',TS);
            assignin('base','dt',dt);
            assignin('base','tstart',tstart);
            assignin('base','tend',tend);
            testCase.verifyWarningFree( ...
                @() evalin('base',sprintf("out = sim('%s.slx');",testCase.modelName)), ...
                'Simulink simulation raised a warning.');
            out = evalin('base','out');

            X = test_plasmaless_imas.runImasLoop(model, TS.Time, TS.Data, x0);
            testCase.assertEqual(out.simout.Time, TS.Time, 'AbsTol', dt*1e-9);
            testCase.assertSize(out.simout.Data, size(X));
            relErr = max(abs(X - out.simout.Data),[],'all') / max(abs(out.simout.Data),[],'all');
            fprintf('IMAS loop vs Simulink: max relative difference %.3g\n', relErr);
            testCase.verifyLessThanOrEqual(relErr, 1e-9, ...
                sprintf('IMAS loop differs from Simulink (relative %.3g).', relErr));
        end

        function test_contract_errors(testCase)
            dt = 0.01;
            model = testCase.buildModel(dt);
            x = zeros(model.n_coils+model.n_loops,1);
            pf_in = test_plasmaless_imas.inputSlice(model, 0, zeros(1,model.n_coils));
            [~, ~, x1] = plasmaless_step_ids(model, pf_in, x);
            testCase.verifyEqual(x1, x);

            bad = pf_in; bad.coil{3}.voltage.data = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:voltage');
            bad = pf_in; bad.coil{3}.voltage.data = [1 2];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:voltage');
            bad = pf_in; bad.coil(5) = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:coils');
            bad = pf_in; bad.coil{5}.name = 'XX1';
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:coils');
            bad = pf_in; bad.time = [0 dt];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:slice');
            bad = pf_in; bad.ids_properties.homogeneous_time = 0;
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:slice');
            testCase.verifyError(@() plasmaless_step_ids(model, pf_in, x(1:end-1)), 'plasmaless:state');

            bad = testCase.pf_active; bad.coil{13}.resistance = -9e40;
            testCase.verifyError(@() configure_plasmaless_model(testCase.em_coupling, ...
                bad, testCase.pf_passive, dt), 'plasmaless:resistance');
            bad = testCase.pf_passive; bad.loop(end) = [];
            testCase.verifyError(@() configure_plasmaless_model(testCase.em_coupling, ...
                testCase.pf_active, bad, dt), 'plasmaless:size');
            bad = testCase.em_coupling;
            names = cellfun(@(c) strtrim(c.name), bad.coupling_matrix, 'UniformOutput', false);
            bad.coupling_matrix(strcmp(names,'mutual_active_active')) = [];
            testCase.verifyError(@() configure_plasmaless_model(bad, ...
                testCase.pf_active, testCase.pf_passive, dt), 'plasmaless:emCoupling');

            % warm start from a one-slice output pair reproduces the state
            [pfa_out, pfp_out, x2] = plasmaless_step_ids(model, ...
                test_plasmaless_imas.inputSlice(model, 0, ones(1,model.n_coils)), x);
            [~, x0] = configure_plasmaless_model(testCase.em_coupling, ...
                testCase.pf_active, testCase.pf_passive, dt, pfa_out, pfp_out);
            testCase.verifyEqual(x0, x2);
            pfa_out.coil{2}.current.data = [];
            testCase.verifyError(@() configure_plasmaless_model(testCase.em_coupling, ...
                testCase.pf_active, testCase.pf_passive, dt, pfa_out, pfp_out), ...
                'plasmaless:data');
            testCase.verifyError(@() configure_plasmaless_model(testCase.em_coupling, ...
                testCase.pf_active, testCase.pf_passive, dt, pfa_out), 'plasmaless:nargin');
        end

        function test_resistive_coil(testCase)
            % Drive VS3U/VS3L with 5 V from zero currents: the tail
            % currents settle to V/R (as test_plasmaless_model).
            V = 5; tstart = 0; tend = 3; dt = 0.01;
            tol_rel = 0.01; tol_abs = 1e-6;
            model = testCase.buildModel(dt);
            drivenCoils = find(ismember(model.coil_names, {'VS3U','VS3L'}));
            testCase.assertNumElements(drivenCoils, 2);
            resistance = diag(model.R(drivenCoils,drivenCoils))';
            testCase.assertTrue(all(isfinite(resistance) & resistance > 0), ...
                'Driven coils must have finite positive resistance.');
            expectedCurrent = V ./ resistance;

            t = (tstart:dt:tend)';
            u = zeros(numel(t),model.n_coils);
            u(:,drivenCoils) = V;
            X = test_plasmaless_imas.runImasLoop(model, t, u, ...
                zeros(model.n_coils+model.n_loops,1));
            tailMask = t >= tstart + 0.8*(tend-tstart);
            for k = 1:numel(drivenCoils)
                maxError = max(abs(X(tailMask,drivenCoils(k)) - expectedCurrent(k)));
                testCase.verifyLessThanOrEqual(maxError, ...
                    tol_abs + tol_rel*abs(expectedCurrent(k)), ...
                    sprintf('Coil %s: tail current not within 1%% of V/R = %.6g A (max error %.6g A).', ...
                    model.coil_names{drivenCoils(k)}, expectedCurrent(k), maxError));
            end
        end

    end

    methods (Access = private)
        function model = buildModel(testCase, dt)
            model = configure_plasmaless_model(testCase.em_coupling, ...
                testCase.pf_active, testCase.pf_passive, dt);
        end
    end

    methods (Static, Access = private)
        function pf_in = inputSlice(model, t, u)
            % one-slice pf_active with coil voltages u at time t
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.coil = ids_allocate('pf_active','coil',model.n_coils);
            for ii = 1:model.n_coils
                pf_in.coil{ii}.name = model.coil_names{ii};
                pf_in.coil{ii}.voltage.data = u(ii);
            end
        end

        function X = runImasLoop(model, t, u, x0)
            % state history (one row per time) from plasmaless_step_ids,
            % reading the currents from the output IDSs
            X = zeros(numel(t), model.n_coils+model.n_loops);
            X(1,:) = x0';
            x = x0;
            for k = 1:numel(t)-1
                pf_in = test_plasmaless_imas.inputSlice(model, t(k), u(k,:));
                [pfa_out, pfp_out, x] = plasmaless_step_ids(model, pf_in, x);
                X(k+1,:) = [cellfun(@(c) c.current.data, pfa_out.coil)' ...
                    cellfun(@(l) l.current, pfp_out.loop)'];
            end
        end
    end
end
