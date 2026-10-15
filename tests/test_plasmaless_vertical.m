classdef test_plasmaless_vertical < matlab.unittest.TestCase
    % IMAS interface of the vertical-instability model
    % (vertical_instability/configure_vertical_instability_model.m with the
    % generalised extract_equilibrium_data.m).
    %
    % Run from the repository root with an IMAS-MATLAB DD-4.x module, reading an
    % IMAS entry first with MATLAB R2024a (see tests/test_plasmaless_imas.m):
    % add_plasmaless_paths; load_plasmaless_md(); runtests('tests/test_plasmaless_vertical.m')
    % Data: DD4 PDS scenario entry 105073 at 84.4176 s (j_phi stored, no passive
    % currents, open stored LCFS) and the DD3 DINA entry 105108/2 at 100 s (the
    % model author's default).
    % Seven tests
    % 1) test_build_dd4 (zero passive currents: one unstable pole, common
    %    model fields, discrete C row, zc0 vs the entry's current_centre.z)
    % 2) test_build_dd3_default (author's default entry, passive currents of the entry)
    % 3) test_step_loop (plasmaless_step_ids vs the direct recursion; Ip, rc
    %    constant, zc = zc0 + C row * x, loop currents = Ic0 + deviation)
    % 4) test_observables_round_trip (Ip, rc, zc through the equilibrium IDS)
    % 5) test_passive_currents_option (zero passive currents only on request)
    % 6) test_entry_path (plain directory / equilibrium.h5 path -> imas:hdf5 URI)
    % 7) test_circuit_linear_unchanged (op.loop_current zero, same currents)
    properties
        rootDir
        em_coupling
        pf_active
        pf_passive
        eq4
        model4
    end

    properties (Constant)
        dd4Dir = '/work/projects/pds/pds-scenarios/105073/data/in'
        dd4Time = 84.4176
        dt = 0.001
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(testCase.rootDir, fullfile(testCase.rootDir,'imas_model'), ...
                fullfile(testCase.rootDir,'vertical_instability'));
            [testCase.em_coupling, testCase.pf_active, testCase.pf_passive] = load_plasmaless_md();
            testCase.eq4 = extract_equilibrium_data(testCase.dd4Dir, testCase.dd4Time, ...
                struct('passive_currents','zero','pf_passive_md',testCase.pf_passive));
            testCase.model4 = configure_vertical_instability_model(testCase.dt, testCase.eq4, ...
                testCase.em_coupling, testCase.pf_active, testCase.pf_passive);
        end
    end

    methods (Test)

        function test_build_dd4(testCase)
            eq = testCase.eq4; model = testCase.model4;
            testCase.verifyEqual(eq.uri, ['imas:hdf5?path=' testCase.dd4Dir]);
            testCase.verifyEqual(eq.passive_source, 'zero (not in entry)');
            % DD4 (COCOS 17): psi0 = +psi, checked by Ampere's law
            testCase.verifyEqual(eq.psi_sign, 1);
            testCase.verifyEqual(sign(eq.ampere_ip), sign(eq.Ip0));
            testCase.verifyTrue(startsWith(eq.boundary_source, 'psi contour'));
            n = model.n_coils + model.n_loops;
            testCase.verifyEqual(model.kind, 'vertical');
            testCase.verifyEqual([model.n_coils model.n_loops model.n_states], [14 102 116]);
            testCase.verifyEqual(model.x_to_coil, [eye(14) zeros(14,102)]);
            testCase.verifyEqual(model.x_to_loop, [zeros(102,14) eye(102)]);
            testCase.verifyEqual(model.y_names, {'Ip','rc','zc'});
            testCase.verifyEqual(model.op.y_known, true(3,1));
            testCase.verifyEqual(model.op.y, [eq.Ip0; model.rc0; model.zc0]);
            testCase.verifyEqual(model.op.coil_current, model.Ic0(1:14));
            testCase.verifyEqual(model.op.loop_current, zeros(102,1));
            testCase.verifyEqual(model.op.voltage, zeros(14,1));
            testCase.verifyFalse(any(model.op.voltage_given));
            % discrete C: zIp row = sysd.C(1,:); C_y = [0; 0; zIp row / Ip0]
            testCase.verifyEqual(model.C(model.ind_zIp,:), model.sysd.C(1,:));
            testCase.verifyEqual(model.C_y, [zeros(2,n); model.sysd.C(1,:)/model.Ip0]);
            poles = eig(model.Ac);
            tol = 1e-8*max(1,max(abs(poles)));
            testCase.verifyEqual(nnz(real(poles) > tol), 1);
            testCase.verifyEqual(nnz(abs(eig(model.A)) > 1+1e-8), 1);
            growth = max(real(poles));
            dz = model.zc0 - eq.current_centre_entry(2);
            testCase.log(1, sprintf(['DD4 105073 at %.10g s: growth rate %.6g 1/s, zc0 %.6f m, ', ...
                'entry current_centre.z %.6f m (difference %.3g m), rc0 %.6f m, entry r %.6f m'], ...
                testCase.dd4Time, growth, model.zc0, eq.current_centre_entry(2), dz, ...
                model.rc0, eq.current_centre_entry(1)));
            fprintf('test_build_dd4: growth rate %.6g 1/s, zc0 - current_centre.z = %.3g m\n', growth, dz);
            testCase.verifyLessThan(abs(dz), 0.01);
        end

        function test_build_dd3_default(testCase)
            eq = extract_equilibrium_data([], [], struct('pf_passive_md',testCase.pf_passive));
            testCase.verifyEqual(eq.passive_source, 'equilibrium entry');
            testCase.verifyEqual(eq.psi_sign, -1);
            testCase.verifyEqual(eq.boundary_source, 'boundary.outline');
            model = configure_vertical_instability_model(testCase.dt, eq, ...
                testCase.em_coupling, testCase.pf_active, testCase.pf_passive);
            poles = eig(model.Ac);
            tol = 1e-8*max(1,max(abs(poles)));
            testCase.verifyEqual(nnz(real(poles) > tol), 1);
            testCase.verifyNotEqual(model.op.loop_current, zeros(102,1));
            dz = model.zc0 - eq.current_centre_entry(2);
            fprintf('test_build_dd3_default: growth rate %.6g 1/s, zc0 - current_centre.z = %.3g m\n', ...
                max(real(poles)), dz);
            testCase.verifyLessThan(abs(dz), 0.01);
        end

        function test_step_loop(testCase)
            model = testCase.model4;
            n_steps = 300;
            % 1 V on VS3U and -1 V on VS3L from the first step (deviation input)
            u = zeros(model.n_coils,1);
            u(strcmp(model.coil_names,'VS3U')) = 1;
            u(strcmp(model.coil_names,'VS3L')) = -1;
            x = zeros(model.n_states,1); xd = x;
            zc = zeros(n_steps,1);
            for k = 1:n_steps
                pf_in = test_plasmaless_vertical.inputSlice(model, (k-1)*model.dt, u');
                [pfa, pfp, x, y_out] = plasmaless_step_ids(model, pf_in, x);
                xd = model.A*xd + model.B*u;
                testCase.assertEqual(x, xd);
                testCase.assertEqual(cellfun(@(c) c.current.data, pfa.coil(:)), ...
                    model.op.coil_current + xd(1:14), 'AbsTol', 1e-9);
                testCase.assertEqual(cellfun(@(l) l.current, pfp.loop(:)), ...
                    model.op.loop_current + xd(15:end), 'AbsTol', 1e-9);
                testCase.assertEqual(y_out.names, {'Ip','rc','zc'});
                testCase.assertEqual(y_out.values(1:2), [model.Ip0; model.rc0]);
                zIp = model.C(model.ind_zIp,:)*xd;
                testCase.assertEqual(y_out.values(3), model.zc0 + zIp/model.Ip0, 'AbsTol', 1e-12);
                zc(k) = y_out.values(3);
            end
            testCase.verifyEqual(numel(pfp.loop), 102);
            testCase.verifyTrue(startsWith(pfa.ids_properties.comment, 'plasmaless vertical-instability model'));
            % open loop: the displacement grows
            dz = abs(zc - model.zc0);
            testCase.verifyGreaterThan(dz(end), dz(n_steps/2));
            testCase.verifyGreaterThan(dz(end), 0);
        end

        function test_observables_round_trip(testCase)
            names = {'Ip','rc','zc'};
            y.names = names; y.values = [-1.5e7; 6.12; 0.59]; y.known = true(3,1);
            eq = ids_init('equilibrium');
            eq.ids_properties.homogeneous_time = 1;
            eq.time = 1;
            eq.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq.time_slice{1}.time = 1;
            eq = plasmaless_observables_to_equilibrium(eq, y);
            s = eq.time_slice{1};
            testCase.verifyEqual(s.global_quantities.ip, -1.5e7);
            testCase.verifyEqual(s.global_quantities.current_centre.r, 6.12);
            testCase.verifyEqual(s.global_quantities.current_centre.z, 0.59);
            [values, known] = plasmaless_observables_from_equilibrium(s, names);
            testCase.verifyEqual(known, true(3,1));
            testCase.verifyEqual(values, y.values);
            % serialised and read back, as sent by the actor
            eq2 = imas_deserialize(imas_serialize(eq,'equilibrium'),'equilibrium');
            [v2, k2] = plasmaless_observables_from_equilibrium(eq2.time_slice{1}, names);
            testCase.verifyEqual(k2, true(3,1));
            testCase.verifyEqual(v2, y.values);
        end

        function test_passive_currents_option(testCase)
            % the entry has no passive currents: error unless 'zero' is explicit
            testCase.verifyError(@() extract_equilibrium_data(testCase.dd4Dir, testCase.dd4Time, ...
                struct('pf_passive_md',testCase.pf_passive)), 'vertical:passiveCurrents');
            testCase.verifyError(@() extract_equilibrium_data(testCase.dd4Dir, testCase.dd4Time, ...
                struct('passive_currents','none','pf_passive_md',testCase.pf_passive)), 'vertical:options');
            testCase.verifyError(@() extract_equilibrium_data(testCase.dd4Dir, testCase.dd4Time, ...
                struct('passive','zero')), 'vertical:options');
        end

        function test_entry_path(testCase)
            eq = extract_equilibrium_data(fullfile(testCase.dd4Dir,'equilibrium.h5'), testCase.dd4Time, ...
                struct('passive_currents','zero','pf_passive_md',testCase.pf_passive));
            testCase.verifyEqual(eq.uri, ['imas:hdf5?path=' testCase.dd4Dir]);
            testCase.verifyEqual(eq.Jtor, testCase.eq4.Jtor);
            testCase.verifyError(@() extract_equilibrium_data(testCase.rootDir, testCase.dd4Time, ...
                struct('passive_currents','zero','pf_passive_md',testCase.pf_passive)), 'vertical:uri');
        end

        function test_circuit_linear_unchanged(testCase)
            step = 0.01;
            circuit = configure_plasmaless_model(testCase.em_coupling, testCase.pf_active, ...
                testCase.pf_passive, step);
            testCase.verifyEqual(circuit.op.loop_current, zeros(circuit.n_loops,1));
            linear = configure_linear_model(fullfile(testCase.rootDir, ...
                'linearModelSimple_105073_t84.417598'), step, testCase.pf_active);
            testCase.verifyEqual(linear.op.loop_current, zeros(linear.n_loops,1));
            TS = timeseries_plasmaless_model(0, 20*step, step);
            for model = {circuit, linear}
                m = model{1};
                x = zeros(m.n_states,1);
                for k = 1:20
                    pf_in = test_plasmaless_vertical.inputSlice(m, TS.Time(k), TS.Data(k,:));
                    [pfa, pfp, x_new] = plasmaless_step_ids(m, pf_in, x);
                    testCase.assertEqual(x_new, m.A*x + m.B*(TS.Data(k,:)' - m.op.voltage));
                    testCase.assertEqual(cellfun(@(c) c.current.data, pfa.coil(:)), ...
                        m.op.coil_current + m.x_to_coil*x_new);
                    testCase.assertEqual(cellfun(@(l) l.current, pfp.loop(:)), m.x_to_loop*x_new);
                    x = x_new;
                end
            end
        end
    end

    methods (Static)
        function pf_in = inputSlice(model, t, u)
            % one-slice pf_active with coil voltages u at time t (coils by name)
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.coil = ids_allocate('pf_active','coil',model.n_coils);
            for ii = 1:model.n_coils
                pf_in.coil{ii}.name = model.coil_names{ii};
                pf_in.coil{ii}.voltage.data = u(ii);
            end
        end
    end
end
