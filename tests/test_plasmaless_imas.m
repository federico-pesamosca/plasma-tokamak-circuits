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
    % The model input is the 14 power-supply voltages (pf_active.supply, names
    % model.supply_names); a slice with coil voltages only is mapped with
    % model.coil_to_input (VS1 = 0). Three tests
    % 1) test_simulink_equivalence (IMAS step loop vs simulator_plasmaless_model_simulink.slx)
    % 2) test_contract_errors (invalid supply/coil inputs raise errors, coil-voltage
    %    fallback = supply slice with coil_to_input*v, supply output of a step, warm start)
    % 3) test_resistive_coil (expected step response of the VS3 supply on VS3U/VS3L,
    %    as in test_plasmaless_model, from a supply slice and from a coil slice)
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
            x0 = zeros(model.n_states,1);

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

            % the Simulink output is C*x = [Ia; Ie] (14 coil + 102 loop currents),
            % runImasLoop reads the same [Ia Ie] from the output IDSs
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
            x = zeros(model.n_states,1);
            pf_in = test_plasmaless_imas.inputSlice(model, 0, zeros(1,model.n_inputs));
            [~, ~, x1] = plasmaless_step_ids(model, pf_in, x);
            testCase.verifyEqual(x1, x);

            % supply slice: voltages and supply set
            bad = pf_in; bad.supply{3}.voltage.data = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:voltage');
            bad = pf_in; bad.supply{3}.voltage.data = [1 2];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:voltage');
            bad = pf_in; bad.supply(5) = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:supplies');
            bad = pf_in; bad.supply(end+1) = bad.supply(5);
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:supplies');
            bad = pf_in; bad.supply{5}.name = 'XX1-1-2';
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:supplies');
            bad = pf_in; bad.time = [0 dt];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:slice');
            bad = pf_in; bad.ids_properties.homogeneous_time = 0;
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:slice');
            testCase.verifyError(@() plasmaless_step_ids(model, pf_in, x(1:end-1)), 'plasmaless:state');
            % coil slice (no supplies): coil set and voltages
            pf_coil = test_plasmaless_imas.inputSliceCoils(model, 0, zeros(1,model.n_coils));
            bad = pf_coil; bad.coil(5) = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:coils');
            bad = pf_coil; bad.coil{5}.name = 'XX1';
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:coils');
            bad = pf_coil; bad.coil{3}.voltage.data = [];
            testCase.verifyError(@() plasmaless_step_ids(model, bad, x), 'plasmaless:voltage');

            % coil-voltage fallback: a coil slice with v gives the same step as
            % a supply slice with coil_to_input*v; the VS1 input is then 0
            v = (1:model.n_coils)*0.5;
            u = model.coil_to_input*v';
            testCase.verifyEqual(u(strcmp(model.input_names,'VS1')), 0);
            testCase.verifyEqual(u(strcmp(model.input_names,'VS3')), ...
                v(strcmp(model.coil_names,'VS3L')) - v(strcmp(model.coil_names,'VS3U')));
            [pfa_c, pfp_c, x_c] = plasmaless_step_ids(model, ...
                test_plasmaless_imas.inputSliceCoils(model, 0, v), x);
            [pfa_s, pfp_s, x_s] = plasmaless_step_ids(model, ...
                test_plasmaless_imas.inputSlice(model, 0, u'), x);
            testCase.verifyEqual(x_c, x_s);
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_c.coil), ...
                cellfun(@(c) c.current.data, pfa_s.coil));
            testCase.verifyEqual(cellfun(@(l) l.current, pfp_c.loop), ...
                cellfun(@(l) l.current, pfp_s.loop));
            % the supply output carries the applied voltages with model.supply_names;
            % the coil voltages are echoed only for the coil slice
            for pfa = {pfa_c, pfa_s}
                testCase.verifyEqual(numel(pfa{1}.supply), model.n_inputs);
                testCase.verifyEqual(cellfun(@(s) strtrim(s.name), pfa{1}.supply(:)', ...
                    'UniformOutput', false), model.supply_names);
                testCase.verifyEqual(cellfun(@(s) s.voltage.data, pfa{1}.supply(:)), u);
                testCase.verifyEqual(pfa{1}.time, dt);
            end
            testCase.verifyEqual(cellfun(@(c) c.voltage.data, pfa_c.coil(:)), v');
            testCase.verifyTrue(all(cellfun(@(c) isempty(c.voltage.data), pfa_s.coil)));
            % the machine-description supply names ('CS3U-1-2', ...) are used
            testCase.verifyEqual(model.supply_names, cellfun(@(s) strtrim(s.name), ...
                testCase.pf_active.supply(:)', 'UniformOutput', false));

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

            % warm start from a one-slice output pair reproduces the state (the
            % step output coil currents satisfy the series/antiseries constraints
            % by construction, Ia = T*Ic)
            [pfa_out, pfp_out, x2] = plasmaless_step_ids(model, ...
                test_plasmaless_imas.inputSlice(model, 0, ones(1,model.n_inputs)), x);
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
            % Drive the VS3 supply with 5 V from zero currents: the tail current
            % of the VS3 circuit settles to V/Rcc (Rcc = R_VS3U + R_VS3L, as
            % test_plasmaless_model), VS3L -> +V/Rcc and VS3U -> -V/Rcc
            % (antiseries). Same result from a coil slice with VS3L = 5 V,
            % VS3U = 0 V (Vps_VS3 = V_VS3L - V_VS3U).
            V = 5; tstart = 0; tend = 3; dt = 0.01;
            tol_rel = 0.01; tol_abs = 1e-6;
            model = testCase.buildModel(dt);
            circ = find(strcmp(model.circuit_names, 'VS3'));
            inp = find(strcmp(model.input_names, 'VS3'));
            iU = find(strcmp(model.coil_names, 'VS3U'));
            iL = find(strcmp(model.coil_names, 'VS3L'));
            testCase.assertNumElements([circ inp iU iL], 4);
            Rcc = model.R(circ,circ);
            testCase.assertEqual(Rcc, model.R_coils(iU,iU) + model.R_coils(iL,iL), ...
                'RelTol', 1e-12, 'Rcc must be the sum of the VS3U and VS3L resistances.');
            testCase.assertTrue(isfinite(Rcc) && Rcc > 0, ...
                'The driven circuit must have finite positive resistance.');
            drivenCoils = [iL iU];
            expectedCurrent = [V/Rcc, -V/Rcc];

            t = (tstart:dt:tend)';
            u = zeros(numel(t),model.n_inputs);
            u(:,inp) = V;
            X = test_plasmaless_imas.runImasLoop(model, t, u, zeros(model.n_states,1));
            tailMask = t >= tstart + 0.8*(tend-tstart);
            for k = 1:numel(drivenCoils)
                maxError = max(abs(X(tailMask,drivenCoils(k)) - expectedCurrent(k)));
                testCase.verifyLessThanOrEqual(maxError, ...
                    tol_abs + tol_rel*abs(expectedCurrent(k)), ...
                    sprintf('Coil %s: tail current not within 1%% of V/Rcc = %.6g A (max error %.6g A).', ...
                    model.coil_names{drivenCoils(k)}, expectedCurrent(k), maxError));
            end
            % coil slice: VS3L = 5 V, VS3U = 0 V gives the same currents
            v = zeros(numel(t),model.n_coils);
            v(:,iL) = V;
            Xc = test_plasmaless_imas.runImasLoop(model, t, v, zeros(model.n_states,1), true);
            testCase.verifyEqual(Xc, X);
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
            % one-slice pf_active with the supply voltages u (model.n_inputs,
            % names model.supply_names) at time t, no coils
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.supply = ids_allocate('pf_active','supply',model.n_inputs);
            for jj = 1:model.n_inputs
                pf_in.supply{jj}.name = model.supply_names{jj};
                pf_in.supply{jj}.voltage.data = u(jj);
            end
        end

        function pf_in = inputSliceCoils(model, t, v)
            % one-slice pf_active with coil voltages v at time t (coils by
            % name, no supplies): mapped with model.coil_to_input
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.coil = ids_allocate('pf_active','coil',model.n_coils);
            for ii = 1:model.n_coils
                pf_in.coil{ii}.name = model.coil_names{ii};
                pf_in.coil{ii}.voltage.data = v(ii);
            end
        end

        function X = runImasLoop(model, t, u, x0, coils)
            % current history [Ia Ie] (one row per time) from plasmaless_step_ids,
            % reading the currents from the output IDSs; u is nt x n_inputs
            % (supply slices) or, with coils = true, nt x n_coils (coil slices)
            if nargin < 5, coils = false; end
            X = zeros(numel(t), model.n_coils+model.n_loops);
            X(1,:) = [(model.op.coil_current + model.x_to_coil*x0)' ...
                (model.op.loop_current + model.x_to_loop*x0)'];
            x = x0;
            for k = 1:numel(t)-1
                if coils
                    pf_in = test_plasmaless_imas.inputSliceCoils(model, t(k), u(k,:));
                else
                    pf_in = test_plasmaless_imas.inputSlice(model, t(k), u(k,:));
                end
                [pfa_out, pfp_out, x] = plasmaless_step_ids(model, pf_in, x);
                X(k+1,:) = [cellfun(@(c) c.current.data, pfa_out.coil)' ...
                    cellfun(@(l) l.current, pfp_out.loop)'];
            end
        end
    end
end
