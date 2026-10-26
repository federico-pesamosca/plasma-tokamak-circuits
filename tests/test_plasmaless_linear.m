classdef test_plasmaless_linear < matlab.unittest.TestCase
    % Unit tests for the external linear model (imas_model/configure_linear_model.m)
    % with the received NICE linearisation linearModelSimple_105073_t84.417598
    % (pulse 105073, confirmed by the user on 2026-10-14, delivered under the label 105078;
    % descriptor with a sign-convention block).
    %
    % Run from the repository root with an IMAS-MATLAB DD-4.x module, reading an
    % IMAS entry first with MATLAB R2024a (see tests/test_plasmaless_imas.m):
    % add_plasmaless_paths; load_plasmaless_md(); runtests('tests/test_plasmaless_linear.m')
    % The model inputs are the descriptor's 14 generator (power-supply) voltages,
    % B in input space; input_map (= coil_to_input) is the coil-voltage fallback.
    % Eight tests
    % 1) test_build (sizes, input-space B, generator-to-coil input map, supply names,
    %    state mapping, sign convention of the descriptor (B unchanged, Ip row of C
    %    unchanged, geometric rows negated), operating point of the descriptor (IMAS
    %    convention, voltage in input space), D = 0, stability, read_npy vs the .txt files)
    % 2) test_step_loop (plasmaless_step_ids vs the direct recursion x_{k+1} = A x_k + B u_k)
    % 3) test_operating_point (du = u - op.voltage, output offsets, supply output)
    % 4) test_descriptor_errors (broken descriptors raise the documented errors,
    %    incl. the convention block and the operating_point.voltage key sets; without
    %    convention the file matrices are used as is)
    % 5) test_input_map_compat (descriptor without input_map: coil inputs, permutation)
    % 6) test_observables_round_trip (equilibrium IDS write/read)
    % 7) test_circuit_unchanged (circuit model selections and 4-output step)
    % 8) test_one_pulse_per_generator (input identification: one bipolar pulse per
    %    generator on the continuous file matrices, which currents react)
    properties
        rootDir
        modelDir
        em_coupling
        pf_active
        pf_passive
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(testCase.rootDir, fullfile(testCase.rootDir,'imas_model'));
            testCase.modelDir = fullfile(testCase.rootDir,'linearModelSimple_105073_t84.417598');
            md_uri = getenv('PLASMALESS_MD_URI');   % empty: data/md_dd4 default
            [testCase.em_coupling, testCase.pf_active, testCase.pf_passive] = ...
                load_plasmaless_md(md_uri);
        end
    end

    methods (Test)

        function test_build(testCase)
            dt = 0.01;
            [model, x0] = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            testCase.verifyEqual(model.kind, 'linear');
            testCase.verifyEqual(model.n_states, 17);
            testCase.verifyEqual(model.n_coils, 14);
            testCase.verifyEqual(model.n_loops, 4);
            testCase.verifyEqual(numel(model.y_names), 9);
            testCase.verifySize(model.A, [17 17]);
            testCase.verifySize(model.B, [17 14]);
            testCase.verifySize(model.C_y, [9 17]);
            testCase.verifySize(model.x_to_coil, [14 17]);
            testCase.verifySize(model.x_to_loop, [4 17]);
            testCase.verifySize(x0, [17 1]);
            testCase.verifyEqual(x0, zeros(17,1));
            testCase.verifyEqual(model.coil_names, {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','PF6','VS3U','VS3L'});
            % inputs: the 14 generator voltages (model provider, 2026-10-14), the
            % model inputs; supply names from the machine description by key
            testCase.verifyEqual(model.n_inputs, 14);
            testCase.verifyEqual(model.input_names, {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'});
            testCase.verifyEqual(model.supply_names, cellfun(@(s) strtrim(s.name), ...
                testCase.pf_active.supply(:)', 'UniformOutput', false));
            testCase.verifyEqual(model.supply_names{12}, 'VS1-1-2-3-4-5-6');
            testCase.verifyEqual(model.supply_names{14}, 'VS3-1');
            testCase.verifyEqual(model.state_names, {'CS3U','CS2U','CS1','CS2L','CS3L', ...
                'PF1','PF2','PF3','PF4','PF5','PF6','VS3','-VS3','PS_1','PS_2','VV_1','VV_2'});
            testCase.verifyEqual(model.loop_names, {'PS_1','PS_2','VV_1','VV_2'});
            iU = find(strcmp(model.coil_names,'CS1U')); iL = find(strcmp(model.coil_names,'CS1L'));
            testCase.verifyEqual(model.x_to_coil(iU,:), model.x_to_coil(iL,:));
            testCase.verifyEqual(model.x_to_coil(strcmp(model.coil_names,'VS3U'),:), ...
                double(strcmp(model.state_names,'VS3')));
            testCase.verifyEqual(model.x_to_coil(strcmp(model.coil_names,'VS3L'),:), ...
                double(strcmp(model.state_names,'-VS3')));
            testCase.verifyEqual(sum(model.x_to_coil,2), ones(14,1));
            testCase.verifyEqual(sum(model.x_to_loop,2), ones(4,1));
            testCase.verifyEqual(model.D, zeros(9,14));
            testCase.verifyEqual(model.dt, dt);
            testCase.verifyLessThanOrEqual(max(real(eig(model.A_c))), 1e-9);
            % the same matrix as readmatrix of the 6-digit .txt files
            A_txt = readmatrix(fullfile(testCase.modelDir,'linModSimpleA.txt'));
            A_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleA.npy'));
            testCase.verifyEqual(A_npy, model.A_c);
            testCase.verifyLessThanOrEqual(max(abs(A_npy - A_txt),[],'all'), ...
                1e-4*max(abs(A_txt),[],'all'));
            C_txt = readmatrix(fullfile(testCase.modelDir,'linModSimpleC.txt'));
            C_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleC.npy'));
            testCase.verifyEqual(model.C_file, C_npy);
            testCase.verifyLessThanOrEqual(max(abs(C_npy - C_txt),[],'all'), ...
                1e-4*max(abs(C_txt),[],'all'));
            B_txt = readmatrix(fullfile(testCase.modelDir,'linModSimpleB.txt'));
            B_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleB.npy'));
            testCase.verifyLessThanOrEqual(max(abs(B_npy - B_txt),[],'all'), ...
                1e-4*max(abs(B_txt),[],'all'));
            % B in input space (generator columns): B_c = B_c_input = B_npy
            testCase.verifyEqual(model.B_c_input, B_npy);
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifySize(model.B_c, [17 14]);
            % input map (coil-voltage fallback, u = input_map*v_coil): CS/PF
            % generator rows on their coil, PF6 = row 13 on coil 12, VS3 (row 14)
            % = -VS3U +VS3L (generator voltage V(VS3L) - V(VS3U), provider's
            % positive-V-positive-I rule on B, 2026-10-14), VS1 (row 12) on no coil
            testCase.verifySize(model.input_map, [14 14]);
            testCase.verifyEqual(model.coil_to_input, model.input_map);
            row = @(name) model.input_map(strcmp(model.input_names,name),:);
            for name = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6'}
                testCase.verifyEqual(row(name{1}), double(strcmp(model.coil_names,name{1})), name{1});
            end
            testCase.verifyEqual(row('VS3'), double(strcmp(model.coil_names,'VS3L')) ...
                - double(strcmp(model.coil_names,'VS3U')));
            testCase.verifyEqual(model.input_map(12,:), zeros(1,14));   % VS1 undriven
            testCase.verifyEqual(sum(model.input_map ~= 0, 2), [ones(11,1); 0; 1; 2]);
            % sign convention of the descriptor (2026-10-14): state_sign -1,
            % input_sign -1, output_sign Ip -1 and +1 for the geometric outputs.
            % B_conv = (input_sign/state_sign) B = B: the input-space B_c is the file
            % B with no sign change; C_conv(k,:) = (state_sign/output_sign(k)) C(k,:):
            % the Ip row is unchanged, the geometric rows (rmag, zmag, rX, zX,
            % gap1..3, Dist2ndSep) are negated.
            testCase.verifyEqual(model.convention.state_sign, -1);
            testCase.verifyEqual(model.convention.input_sign, -1);
            testCase.verifyEqual(model.convention.output_sign, [1;1;1;1;-1;1;1;1;1]);
            testCase.verifyEqual(model.A_file, A_npy);
            testCase.verifyEqual(model.B_file, B_npy);
            testCase.verifyEqual(model.B_c_input, model.B_file);
            testCase.verifyEqual(model.B_c, model.B_file);
            iIp = find(strcmp(model.y_names,'Ip'));
            testCase.verifyEqual(iIp, 5);
            testCase.verifyEqual(model.C_c(iIp,:), C_npy(5,:));
            testCase.verifyEqual(model.C_c(strcmp(model.y_names,'rmag'),:), -C_npy(1,:));
            testCase.verifyEqual(model.C_c, [-C_npy(1:4,:); C_npy(5,:); -C_npy(6:9,:)]);
            testCase.verifyEqual(model.C_y, model.C_c);   % zero-order hold keeps C
            % building emits the undriven-input warning once (VS1)
            testCase.verifyWarning(@() configure_linear_model(testCase.modelDir, dt, testCase.pf_active), ...
                'plasmaless:inputUndriven');
            % operating point of the descriptor in the IMAS convention: the negated
            % provider values (model provider, 2026-10-14; frame derived from the DINA
            % 105073 data, 2026-10-14)
            op_cur = -[23248.2 -9514.24 -34331.9 -34331.9 -6293.26 37480.8 21873.6 -29058.6 ...
                -19884.8 -40412.3 -10061.6 19263.3 203.564 -203.564]';
            testCase.verifyEqual(model.op.coil_current, op_cur);
            testCase.verifyEqual(model.op.coil_current_given, true(14,1));
            % op.voltage in input space: the descriptor's voltages are keyed by
            % coil (all 0), mapped with input_map; VS1 is not given (zero row)
            testCase.verifyEqual(model.op.voltage, zeros(14,1));
            testCase.verifyEqual(model.op.voltage_given, [true(11,1); false; true; true]);
            y_known = ismember(model.y_names, {'rmag','zmag','Ip'})';
            testCase.verifyEqual(model.op.y_known, y_known);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'rmag')), 6.29704);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'zmag')), 0.596501);
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), -15063200); % provider's 1.50632e7 A, IMAS sign (DINA 105073)
            testCase.verifyEqual(model.source.shot, 105073);
            testCase.verifyEqual(model.op.y(~y_known), zeros(6,1));
        end

        function test_step_loop(testCase)
            dt = 0.01; n_steps = 200;
            model = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            % plain deviation model: zero operating point for this test
            model.op.coil_current(:) = 0; model.op.y(:) = 0; model.op.y_known(:) = false;
            TS = timeseries_plasmaless_model(0, n_steps*dt, dt);
            x = zeros(model.n_states,1);
            X = zeros(n_steps+1, model.n_states);
            Ia = zeros(n_steps+1, model.n_coils);
            Ie = zeros(n_steps+1, model.n_loops);
            Yo = zeros(n_steps+1, numel(model.y_names));
            for k = 1:n_steps
                pf_in = test_plasmaless_linear.inputSlice(model, TS.Time(k), TS.Data(k,:));
                [pfa_out, pfp_out, x, y_out] = plasmaless_step_ids(model, pf_in, x);
                X(k+1,:) = x';
                Ia(k+1,:) = cellfun(@(c) c.current.data, pfa_out.coil);
                Ie(k+1,:) = cellfun(@(l) l.current, pfp_out.loop);
                Yo(k+1,:) = y_out.values';
                testCase.assertEqual(numel(pfp_out.loop), 4);
                testCase.assertEqual(y_out.names, model.y_names);
                testCase.assertEqual(y_out.known, false(9,1));
            end
            testCase.verifyEqual(cellfun(@(l) l.name, pfp_out.loop(:)', 'UniformOutput', false), ...
                model.loop_names);
            % direct recursion with the discrete matrices
            Xd = zeros(n_steps+1, model.n_states);
            xd = zeros(model.n_states,1);
            for k = 1:n_steps
                xd = model.A*xd + model.B*TS.Data(k,:)';
                Xd(k+1,:) = xd';
            end
            scale = max(abs(Xd),[],'all');
            testCase.verifyGreaterThan(scale, 0);
            testCase.verifyLessThanOrEqual(max(abs(X - Xd),[],'all'), 1e-12*scale);
            testCase.verifyLessThanOrEqual(max(abs(Ia - Xd*model.x_to_coil'),[],'all'), 1e-12*scale);
            testCase.verifyLessThanOrEqual(max(abs(Ie - Xd*model.x_to_loop'),[],'all'), 1e-12*scale);
            Yd = Xd*model.C_y';
            testCase.verifyLessThanOrEqual(max(abs(Yo - Yd),[],'all'), 1e-12*max(abs(Yd),[],'all'));
            iU = find(strcmp(model.coil_names,'CS1U')); iL = find(strcmp(model.coil_names,'CS1L'));
            testCase.verifyEqual(Ia(:,iU), Ia(:,iL));
            testCase.verifyTrue(all(isfinite(X(:))));
        end

        function test_operating_point(testCase)
            dt = 0.01;
            model = configure_linear_model(testCase.modelDir, dt, testCase.pf_active);
            model.op.coil_current = (1:14)'*1e3;
            model.op.voltage = (1:14)'*0.5;
            model.op.y = (1:9)'*0.1;
            model.op.y_known = true(9,1);
            x = zeros(model.n_states,1);
            % u = op.voltage: du = 0, the state stays zero and the outputs are op
            pf_in = test_plasmaless_linear.inputSlice(model, 0, model.op.voltage');
            [pfa_out, pfp_out, x1, y_out] = plasmaless_step_ids(model, pf_in, x);
            testCase.verifyEqual(x1, zeros(model.n_states,1));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), model.op.coil_current);
            testCase.verifyEqual(cellfun(@(l) l.current, pfp_out.loop(:)), zeros(4,1));
            testCase.verifyEqual(y_out.values, model.op.y);
            testCase.verifyEqual(y_out.known, true(9,1));
            % the applied supply voltages are written to the output (no coil
            % voltages for a supply slice)
            testCase.verifyEqual(cellfun(@(s) s.voltage.data, pfa_out.supply(:)), model.op.voltage);
            testCase.verifyEqual(cellfun(@(s) strtrim(s.name), pfa_out.supply(:)', 'UniformOutput', false), ...
                model.supply_names);
            testCase.verifyTrue(all(cellfun(@(c) isempty(c.voltage.data), pfa_out.coil)));
            % u = op.voltage + v: same state as the plain model driven by v
            v = 10*ones(1,14);
            pf_in = test_plasmaless_linear.inputSlice(model, 0, model.op.voltage' + v);
            [pfa_out, ~, x2, y_out] = plasmaless_step_ids(model, pf_in, x);
            x_ref = model.B*v';
            testCase.verifyEqual(x2, x_ref, 'AbsTol', 1e-12*max(abs(x_ref)));
            testCase.verifyEqual(cellfun(@(c) c.current.data, pfa_out.coil(:)), ...
                model.op.coil_current + model.x_to_coil*x_ref, 'AbsTol', 1e-9);
            testCase.verifyEqual(y_out.values, model.op.y + model.C_y*x_ref, 'AbsTol', 1e-12);
            % coil-voltage fallback: a coil slice with v_coil is the input
            % input_map*v_coil (VS1 = 0), then du = u - op.voltage
            v_coil = (1:14)*2;
            pf_in = test_plasmaless_linear.inputSliceCoils(model, 0, v_coil);
            [pfa_out, ~, x3] = plasmaless_step_ids(model, pf_in, x);
            u = model.input_map*v_coil';
            testCase.verifyEqual(u(12), 0);
            testCase.verifyEqual(x3, model.B*(u - model.op.voltage), 'AbsTol', 1e-12*max(abs(x_ref)));
            testCase.verifyEqual(cellfun(@(s) s.voltage.data, pfa_out.supply(:)), u);
            testCase.verifyEqual(cellfun(@(c) c.voltage.data, pfa_out.coil(:)), v_coil');
        end

        function test_descriptor_errors(testCase)
            dt = 0.01;
            base = fullfile('/scratch/users', getenv('USER'), 'plasmaless');
            if ~isfolder(base), mkdir(base); end
            d = jsondecode(fileread(fullfile(testCase.modelDir,'linear_model.json')));

            % duplicate input name
            bad = d; bad.inputs{3} = 'CS1L';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % input count differs from the columns of B
            bad = d; bad.inputs{end+1} = 'VS1b';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:size');
            % input_map: unknown coil
            bad = d; bad.input_map.CS3U = struct('XX1', 1);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % input_map: key that is not an input
            bad = d; bad.input_map.XX1 = struct('CS3U', 1);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % input_map: non-finite coefficient (Infinity in the JSON text), then a string
            bad = d; bad.input_map.CS3U = struct('CS3U', Inf);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.input_map.CS3U = struct('CS3U', 'one');
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % input_map not an object
            bad = d; bad.input_map = 1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % missing coil in coil_states
            bad = d; bad.coil_states = rmfield(bad.coil_states,'PF6');
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % coil mapped to an unknown state
            bad = d; bad.coil_states.PF6 = 'PF7';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % unclassified state (VV_2 removed from passive_states)
            bad = d; bad.passive_states = bad.passive_states(1:3);
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % state both coil and passive
            bad = d; bad.passive_states{end+1} = 'PF6';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % missing required field
            bad = rmfield(d,'outputs');
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % unknown output name
            bad = d; bad.outputs{1} = 'rgeo';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:observable');
            % wrong-size matrix
            bad = d; bad.files.A = 'wrongA.txt';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt, 'wrongA.txt', zeros(16,17)), ...
                'plasmaless:size');
            % missing matrix file
            bad = d; bad.files.B = 'missing.npy';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:file');
            % operating point with a key that is neither an input nor a coil
            bad = d; bad.operating_point.voltage.XX1 = 1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % operating point voltage mixing coil and input names
            bad = d; bad.operating_point.voltage.VS1 = 1;   % VS3U/VS3L (coils) + VS1 (input)
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % convention: invalid sign value (not +1/-1), unknown output name in
            % output_sign, invalid output sign, non-object block, unknown field
            bad = d; bad.convention.state_sign = 2;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.input_sign = 'minus';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.output_sign.rgeo = -1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.output_sign.Ip = 0;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention = 'flip';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            bad = d; bad.convention.current_sign = -1;
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
            % no convention block: all signs +1, the file matrices are used unchanged
            none_conv = rmfield(d,'convention');
            model = testCase.buildFrom(base, none_conv, dt);
            B_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleB.npy'));
            C_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleC.npy'));
            testCase.verifyEqual(model.convention.state_sign, 1);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.convention.output_sign, ones(9,1));
            testCase.verifyEqual(model.B_c_input, B_npy);
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifyEqual(model.C_c, C_npy);
            testCase.verifyEqual(model.C_y, C_npy);
            % partial block (only state_sign): B and all C rows negated
            part = none_conv; part.convention = struct('state_sign', -1);
            model = testCase.buildFrom(base, part, dt);
            testCase.verifyEqual(model.convention.input_sign, 1);
            testCase.verifyEqual(model.B_c_input, -B_npy);
            testCase.verifyEqual(model.C_c, -C_npy);
            % missing descriptor
            testCase.verifyError(@() configure_linear_model(tempname(base), dt, testCase.pf_active), ...
                'plasmaless:descriptor');
            % valid operating point entries are applied and flagged (partial object;
            % voltage keyed by a coil name, mapped to input space)
            good = d;
            good.operating_point = struct('voltage', struct('CS3U', 2.5), ...
                'coil_current', struct('PF1', 1e4), 'outputs', struct('Ip', 1.2e7));
            model = testCase.buildFrom(base, good, dt);
            testCase.verifyEqual(model.op.voltage(1), 2.5);
            testCase.verifyEqual(model.op.voltage_given, [true; false(13,1)]);
            % voltage keyed by coil names through the input map: VS3L = 2 V gives
            % the VS3 input (row -VS3U +VS3L) 2 V, given
            good.operating_point.voltage = struct('VS3L', 2);
            model = testCase.buildFrom(base, good, dt);
            testCase.verifyEqual(model.op.voltage, [zeros(13,1); 2]);
            testCase.verifyEqual(model.op.voltage_given, [false(13,1); true]);
            % voltage keyed by input names (direct): VS1 = 3 V is input 12, given
            good.operating_point.voltage = struct('VS1', 3);
            model = testCase.buildFrom(base, good, dt);
            testCase.verifyEqual(model.op.voltage(12), 3);
            testCase.verifyEqual(model.op.voltage_given, [false(11,1); true; false(2,1)]);
            testCase.verifyEqual(model.op.coil_current(strcmp(model.coil_names,'PF1')), 1e4);
            testCase.verifyEqual(model.op.coil_current_given, strcmp(model.coil_names,'PF1')');
            testCase.verifyEqual(model.op.y(strcmp(model.y_names,'Ip')), 1.2e7);
            testCase.verifyEqual(model.op.y_known, strcmp(model.y_names,'Ip')');
            % no operating point: zeros, nothing given
            none = rmfield(d,'operating_point');
            model = testCase.buildFrom(base, none, dt);
            testCase.verifyEqual(model.op.coil_current, zeros(14,1));
            testCase.verifyFalse(any(model.op.coil_current_given));
            testCase.verifyFalse(any(model.op.voltage_given));
            testCase.verifyEqual(model.op.y_known, false(9,1));
        end

        function test_input_map_compat(testCase)
            % descriptor without input_map: inputs must be the coil set and the
            % map is the permutation to machine-description order
            dt = 0.01;
            base = fullfile('/scratch/users', getenv('USER'), 'plasmaless');
            if ~isfolder(base), mkdir(base); end
            d = jsondecode(fileread(fullfile(testCase.modelDir,'linear_model.json')));
            d = rmfield(d,'input_map');
            coils = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5', ...
                'PF6','VS3U','VS3L'};
            d.inputs = fliplr(coils);          % reversed order: a real permutation
            model = testCase.buildFrom(base, d, dt);
            B_npy = read_npy(fullfile(testCase.modelDir,'linModSimpleB.npy'));
            testCase.verifyEqual(model.input_names, fliplr(coils));
            testCase.verifyEqual(model.n_inputs, 14);
            testCase.verifyEqual(model.coil_names, coils);
            % B stays in input space (descriptor order); the permutation is the
            % coil-voltage fallback input_map
            testCase.verifyEqual(model.B_c, B_npy);
            testCase.verifyEqual(model.input_map, fliplr(eye(14)));
            testCase.verifyEqual(model.coil_to_input, fliplr(eye(14)));
            % the CS/PF coil keys match the machine-description supplies
            % (keys of names like 'CS3U-1-2') and take their names; VS3U/VS3L
            % match none and keep their key
            testCase.verifyEqual(model.supply_names, ...
                plasmaless_supply_names(fliplr(coils), testCase.pf_active));
            testCase.verifySize(model.B, [17 14]);
            % inputs that are not the coil set are an error without input_map
            bad = d; bad.inputs{1} = 'VS3';
            testCase.verifyError(@() testCase.buildFrom(base, bad, dt), 'plasmaless:descriptor');
        end

        function test_observables_round_trip(testCase)
            names = {'rmag','zmag','rX','zX','Ip','gap1','gap2','gap3','Dist2ndSep'};
            y.names = names;
            y.values = [6.2; 0.5; 5.1; -3.3; 1.5e7; 0.1; 0.2; 0.3; 0.05];
            y.known = true(9,1);
            eq = ids_init('equilibrium');
            eq.ids_properties.homogeneous_time = 1;
            eq.time = 1;
            eq.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq.time_slice{1}.time = 1;
            eq = plasmaless_observables_to_equilibrium(eq, y);
            s = eq.time_slice{1};
            testCase.verifyEqual(s.global_quantities.magnetic_axis.r, 6.2);
            testCase.verifyEqual(s.global_quantities.ip, 1.5e7);
            testCase.verifyEqual(s.contour_tree.node{1}.critical_type, int32(1));
            testCase.verifyEqual(s.contour_tree.node{1}.r, 5.1);
            testCase.verifyEqual(numel(s.boundary.gap), 3);
            testCase.verifyEqual(s.boundary.gap{3}.name, 'gap3');
            testCase.verifyEqual(s.boundary.gap{3}.value, 0.3);
            [values, known] = plasmaless_observables_from_equilibrium(s, names);
            % Dist2ndSep has no DD4 field: written nowhere, read back unknown
            expected_known = [true(8,1); false];
            testCase.verifyEqual(known, expected_known);
            testCase.verifyEqual(values(expected_known), y.values(expected_known));
            testCase.verifyTrue(isnan(values(end)));
            % partial mask: only the known entries are written
            y2 = y; y2.known = false(9,1); y2.known([2 5 7]) = true;
            eq2 = ids_init('equilibrium');
            eq2.time_slice = ids_allocate('equilibrium','time_slice',1);
            eq2 = plasmaless_observables_to_equilibrium(eq2, y2);
            [v2, k2] = plasmaless_observables_from_equilibrium(eq2.time_slice{1}, names);
            testCase.verifyEqual(k2, y2.known);
            testCase.verifyEqual(v2(k2), y.values(k2));
            testCase.verifyEqual(numel(eq2.time_slice{1}.boundary.gap), 2);
            % the fresh slice has nothing known
            [~, k0] = plasmaless_observables_from_equilibrium( ...
                ids_allocate('equilibrium','time_slice',1), names);
            testCase.verifyEqual(k0, false(9,1));
            testCase.verifyError(@() plasmaless_observables_from_equilibrium(s, {'rgeo'}), ...
                'plasmaless:observable');
        end

        function test_circuit_unchanged(testCase)
            dt = 0.01;
            model = configure_plasmaless_model(testCase.em_coupling, testCase.pf_active, ...
                testCase.pf_passive, dt);
            n = model.n_circuits + model.n_loops;
            testCase.verifyEqual(model.kind, 'circuit');
            testCase.verifyEqual(model.n_states, n);
            testCase.verifyEqual(model.state_names, [model.circuit_names model.loop_names]);
            testCase.verifyEqual(model.x_to_coil, [model.T zeros(model.n_coils,model.n_loops)]);
            testCase.verifyEqual(model.x_to_loop, [zeros(model.n_loops,model.n_circuits) eye(model.n_loops)]);
            testCase.verifyEmpty(model.y_names);
            testCase.verifySize(model.C_y, [0 n]);
            testCase.verifyEqual(model.n_inputs, model.n_ps);
            testCase.verifyEqual(model.input_names, model.ps_names);
            testCase.verifyEqual(model.coil_to_input, model.topology.coil_to_ps);
            testCase.verifyEqual(model.op.voltage, zeros(model.n_inputs,1));
            x = zeros(n,1);
            TS = timeseries_plasmaless_model(0, 20*dt, dt);
            for k = 1:20
                pf_in = test_plasmaless_linear.inputSlice(model, TS.Time(k), TS.Data(k,:));
                [pfa3, pfp3, x3] = plasmaless_step_ids(model, pf_in, x);
                [pfa4, pfp4, x4, y_out] = plasmaless_step_ids(model, pf_in, x);
                testCase.assertEqual(x4, x3);
                testCase.assertEqual(cellfun(@(c) c.current.data, pfa4.coil), ...
                    cellfun(@(c) c.current.data, pfa3.coil));
                testCase.assertEqual(cellfun(@(l) l.current, pfp4.loop), ...
                    cellfun(@(l) l.current, pfp3.loop));
                % exactly the selections of the state (coil currents T*Ic, loop currents Ie)
                testCase.assertEqual(cellfun(@(c) c.current.data, pfa4.coil(:)), model.x_to_coil*x4);
                testCase.assertEqual(cellfun(@(l) l.current, pfp4.loop(:)), model.x_to_loop*x4);
                testCase.assertEqual(x4, model.A*x + model.B*TS.Data(k,:)');
                testCase.assertEmpty(y_out.values);
                x = x4;
            end
            testCase.verifyEqual(pfa4.ids_properties.comment, 'plasmaless coil-vessel model: coil currents');
        end

        function test_one_pulse_per_generator(testCase)
            % input identification in the provider's generator space: the
            % continuous file matrices (A_file, B_file, provider frame; the
            % checks are on relative patterns, so no sign convention is needed),
            % c2d at dt = 1e-3, x0 = 0; per generator +10 V for 0.5 s, -10 V
            % for 0.5 s, then 0.5 s of zero
            dt = 1e-3;
            model = configure_linear_model(testCase.modelDir, 0.01, testCase.pf_active);
            [M, F, X] = test_plasmaless_linear.pulsePerGenerator(model.A_file, model.B_file, dt);
            gen = model.input_names; st = model.state_names;
            % 1) each CS/PF generator drives its own state hardest among the 11
            %    CS/PF coil states, with a positive initial response
            own = {'CS3U','CS3U'; 'CS2U','CS2U'; 'CS1U','CS1'; 'CS1L','CS1'; 'CS2L','CS2L'; ...
                'CS3L','CS3L'; 'PF1','PF1'; 'PF2','PF2'; 'PF3','PF3'; 'PF4','PF4'; ...
                'PF5','PF5'; 'PF6','PF6'};
            for r = 1:size(own,1)
                j = find(strcmp(gen, own{r,1})); i = find(strcmp(st, own{r,2}));
                [~, imax] = max(M(1:11,j));
                testCase.verifyEqual(st{imax}, own{r,2}, sprintf( ...
                    'generator %s: largest CS/PF response on state %s, expected %s', ...
                    own{r,1}, st{imax}, own{r,2}));
                testCase.verifyGreaterThan(F(i,j), 0, sprintf( ...
                    'generator %s: initial response of state %s is not positive (%g A)', ...
                    own{r,1}, own{r,2}, F(i,j)));
                testCase.verifyGreaterThan(model.B_file(i,j), 0, sprintf( ...
                    'generator %s: B_file(%s,%s) is not positive', own{r,1}, own{r,2}, own{r,1}));
            end
            % 2) CS1U and CS1L feed the same state: identical trajectories
            jU = find(strcmp(gen,'CS1U')); jL = find(strcmp(gen,'CS1L'));
            scale = max(abs(X{jU}),[],'all');
            testCase.verifyLessThanOrEqual(max(abs(X{jU} - X{jL}),[],'all'), 1e-12*scale, ...
                'generators CS1U and CS1L: state trajectories differ');
            % 3) VS1 (star point in series with PF2-PF5): measured max |deviation|
            %    (2026-10-14) PF2 13.13, PF3 2.905, PF4 4.504, PF5 3.608 A, initial
            %    signs - - + +; CS coils at most 0.869 A (CS2U), PF1 0.592, PF6 0.473 A
            %    (VS3 83.6, -VS3 82.0 A, not part of this check)
            jV = find(strcmp(gen,'VS1'));
            iStar = find(ismember(st, {'PF2','PF3','PF4','PF5'}));
            iOther = find(ismember(st, {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF6'}));
            for i = iStar
                for o = iOther
                    testCase.verifyGreaterThan(M(i,jV), M(o,jV), sprintf( ...
                        'generator VS1: state %s (%g A) does not exceed state %s (%g A)', ...
                        st{i}, M(i,jV), st{o}, M(o,jV)));
                end
            end
            % 4) VS3: states VS3 and -VS3 respond most among the 13 coil states,
            %    VS3 negative and -VS3 positive at first (B(12,14) < 0 < B(13,14))
            j3 = find(strcmp(gen,'VS3'));
            [~, order] = sort(M(1:13,j3), 'descend');
            testCase.verifyEqual(sort(order(1:2))', [12 13], sprintf( ...
                'generator VS3: largest coil-state responses on %s and %s, expected VS3 and -VS3', ...
                st{order(1)}, st{order(2)}));
            testCase.verifyLessThan(F(12,j3), 0, 'generator VS3: initial response of state VS3 is not negative');
            testCase.verifyGreaterThan(F(13,j3), 0, 'generator VS3: initial response of state -VS3 is not positive');
        end

    end

    methods (Access = private)
        function model = buildFrom(testCase, base, desc, dt, extra_name, extra_matrix)
            % build from a temporary copy of the model directory with descriptor desc
            dir_tmp = tempname(base);
            copyfile(testCase.modelDir, dir_tmp);
            fid = fopen(fullfile(dir_tmp,'linear_model.json'),'w');
            % Inf/NaN written as Infinity/NaN (jsondecode reads them back)
            fprintf(fid,'%s', jsonencode(desc,'ConvertInfAndNaN',false));
            fclose(fid);
            if nargin > 4
                writematrix(extra_matrix, fullfile(dir_tmp,extra_name), 'Delimiter',' ', 'FileType','text');
            end
            cleanup = onCleanup(@() rmdir(dir_tmp,'s'));
            model = configure_linear_model(dir_tmp, dt, testCase.pf_active);
        end
    end

    methods (Static, Access = private)
        function [M, F, X] = pulsePerGenerator(A, B, dt)
            % one bipolar pulse per generator column of B on the continuous
            % (A, B), discretised with c2d (zero-order hold) at dt, from x0 = 0:
            % +10 V for 0.5 s, -10 V for 0.5 s, 0.5 s of zero. M(i,j) = max |x_i|
            % for generator j, F(:,j) = x after the first step, X{j} = trajectory
            n = size(A,1); m = size(B,2);
            sysd = c2d(ss(A, B, eye(n), zeros(n,m)), dt);
            n1 = round(0.5/dt);
            u = [10*ones(1,n1), -10*ones(1,n1), zeros(1,n1)];
            M = zeros(n,m); F = zeros(n,m); X = cell(1,m);
            for j = 1:m
                x = zeros(n,1); Xj = zeros(numel(u), n);
                for k = 1:numel(u)
                    x = sysd.A*x + sysd.B(:,j)*u(k);
                    Xj(k,:) = x';
                end
                M(:,j) = max(abs(Xj),[],1)';
                F(:,j) = Xj(1,:)';
                X{j} = Xj;
            end
        end

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
            % one-slice pf_active with coil voltages v at time t (coils by name,
            % no supplies): mapped with model.coil_to_input (= input_map)
            pf_in = ids_init('pf_active');
            pf_in.ids_properties.homogeneous_time = 1;
            pf_in.time = t;
            pf_in.coil = ids_allocate('pf_active','coil',model.n_coils);
            for ii = 1:model.n_coils
                pf_in.coil{ii}.name = model.coil_names{ii};
                pf_in.coil{ii}.voltage.data = v(ii);
            end
        end
    end
end
