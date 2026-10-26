classdef test_plasmaless_model < matlab.unittest.TestCase
    % Unit tests for the plasmaless coil-vessel circuit model.
    %
    % Simulink (sim()) resolves variables (bus objects, state-space
    % matrices, timeseries, ...) from the BASE workspace, not from the
    % test method's local workspace. Every variable normally created by
    % run_plasmaless_model.m is therefore pushed to the base workspace
    % with assignin(), and the simulation/plot routines are executed
    % with evalin('base', ...) so they see exactly the same context as
    % when run manually.
    %
    % Test coverage:
    % 1) test_basic_run: configure, simulate and plot without warnings.
    % 2) test_matlab_simulink_equivalence: compare MATLAB and Simulink
    %    using timeseries_plasmaless_model and the run script's timing
    % 3) test_superconductive_coil: check CS3U/PF1 step currents (input: the
    %    matching power supply) and the superconductive and full-circuit
    %    bounds on current slope, using the circuit matrices Mcc/Mce/Mee.
    % 4) test_resistive_coil: apply 5 V to the VS3 power supply (VS3U/VS3L
    %    in antiseries) for 3 s; check settling of I_VS3 to V/Rcc, the
    %    antiseries relation and constant tail currents in the other circuits.
    % 5) test_circuit_mapping: series/antiseries and power-supply mappings.

    properties
        modelName = 'simulator_plasmaless_model_simulink';
        rootDir
    end

    methods (TestClassSetup)
        function addModelToPath(testCase)
            % tests/ is one level below the project root
            testCase.rootDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(testCase.rootDir));
        end
    end

    methods (TestMethodSetup)
        function freshBaseWorkspace(testCase)
            % Start every test from a clean base workspace and close any
            % figures left over from a previous test/run.
            evalin('base','clear');
            close all force;
            testCase.addTeardown(@() close('all','force'));
            testCase.addTeardown(@() evalin('base','clear'));
        end
    end

    methods (Test)

        function test_basic_run(testCase)
            % Basic smoke test: configure the model, run the simulation
            % and the plotting routine with a short time horizon, and
            % verify that nothing throws an error/warning.

            % --- same initialization as run_plasmaless_model.m ---
            assignin('base','tstart',0);
            assignin('base','tend',28);   % short horizon -> fast test
            assignin('base','dt',0.01);

            test_plasmaless_model.configureBase(0.01);

            evalin('base', ...
                'TS = timeseries_plasmaless_model(tstart,tend,dt);');

            % --- run model ---
            testCase.verifyWarningFree( ...
                @() evalin('base', ...
                    "out = sim('simulator_plasmaless_model_simulink.slx');"), ...
                'Simulink simulation raised a warning.');

            % sanity check on simulation output
            out = evalin('base','out');
            testCase.verifyClass(out.simout.Data, 'double');
            testCase.verifyGreaterThan(numel(out.simout.Time), 1);

            % --- plot model (reduced number of frames for speed) ---
            evalin('base','Nframes = 3;'); %#ok<NASGU>
            testCase.verifyWarningFree( ...
                @() evalin('base','plot_plasmaless_model;'), ...
                'Plotting routine raised a warning.');
        end

        function test_matlab_simulink_equivalence(testCase)
            % Match the timing in run_plasmaless_model.m. Both simulators
            % receive the same sampled voltages, matrices and initial state.
            assignin('base','tstart',0);
            assignin('base','tend',28);
            assignin('base','dt',0.001);
            test_plasmaless_model.configureBase(0.001);
            evalin('base', ...
                'TS = timeseries_plasmaless_model(tstart,tend,dt);');

            testCase.verifyWarningFree( ...
                @() evalin('base', ...
                    "out_simulink = sim('simulator_plasmaless_model_simulink.slx');"), ...
                'Simulink simulation raised a warning in the comparison test.');
            testCase.verifyWarningFree( ...
                @() evalin('base', ...
                    'out_matlab = simulator_plasmaless_model_matlab(TS,model,x0);'), ...
                'MATLAB simulation raised a warning in the comparison test.');

            outSimulink = evalin('base','out_simulink');
            outMatlab = evalin('base','out_matlab');
            TS = evalin('base','TS');
            simulinkTime = outSimulink.simout.Time(:);
            matlabTime = outMatlab.simout.Time(:);
            simulinkData = outSimulink.simout.Data;
            matlabData = outMatlab.simout.Data;

            % Compare the original sample grids: no interpolation or shifts
            % that could hide an off-by-one timestep error.
            testCase.assertEqual(size(simulinkTime), size(TS.Time));
            testCase.assertEqual(size(matlabTime), size(TS.Time));
            testCase.assertEqual(simulinkTime, TS.Time, 'AbsTol', 1e-10);
            testCase.assertEqual(matlabTime, TS.Time, 'AbsTol', 1e-10);
            testCase.assertEqual(size(simulinkData), size(matlabData));
            testCase.assertEqual(size(matlabData), ...
                [numel(TS.Time), evalin('base','size(C,1)')]);
            testCase.assertTrue(all(isfinite(simulinkData(:))), ...
                'Simulink output contains nonfinite currents.');
            testCase.assertTrue(all(isfinite(matlabData(:))), ...
                'MATLAB output contains nonfinite currents.');

            % Check every sample of every active and passive current.
            % Absolute slack protects zero crossings; relative slack allows
            % small floating-point differences in larger current values.
            tol_abs = 1e-8; % [A]
            tol_rel = 1e-8;
            for k = 1:size(matlabData,2)
                error = abs(simulinkData(:,k) - matlabData(:,k));
                tolerance = tol_abs + tol_rel*abs(matlabData(:,k));
                testCase.verifyLessThanOrEqual(max(error ./ tolerance), 1, ...
                    sprintf(['Output %d: MATLAB/Simulink mismatch ', ...
                    '(max absolute error %.6g A).'], k, max(error)));
            end
        end

        function test_superconductive_coil(testCase)
            % Coils 1-12 are superconductive; coils 13/14 and the
            % passive structures are resistive. With zero initial
            % currents and a positive voltage step on one coil only,
            % the driven current slope decreases from V*(M\e_k)(k)
            % toward V/LeffSC, where LeffSC is the Schur complement
            % over the other superconductive coils. Both bounds are
            % at least V/Lkk. The tail need not have settled by 20s.
            % Drive CS3U and PF1 separately and check current and slope.

            assignin('base','tstart',0);
            assignin('base','tend',20);
            assignin('base','dt',0.05); % coarser dt -> fast test

            test_plasmaless_model.configureBase(0.05);

            V = 1;               % [V] step amplitude
            tol_rel = 1e-6;      % relative numerical slack on slope bounds
            debugme = 0;        % Set to 1 to plot; close the figure to resume.

            if debugme == 1
                debugFig = figure('Name', ...
                    'Superconductive coil currents - close to resume test', ...
                    'NumberTitle', 'off');
                debugLayout = tiledlayout(debugFig, 2, 1);
            end

            for coilName = {'CS3U','PF1'}
                % circuit index (state), power-supply index (input), coil index (output)
                [circIdx, Lkk, M, psIdx, coilIdx] = ...
                    test_plasmaless_model.getCoilInfo(coilName{1});
                testCase.assertLessThanOrEqual(circIdx, 11, ...
                    sprintf(['Circuit %s is expected to be superconductive ', ...
                    '(circuit index <= 11); check model assumptions.'], coilName{1}));

                tstart = evalin('base','tstart');
                tend   = evalin('base','tend');
                dt     = evalin('base','dt');
                t = (tstart:dt:tend)';
                u = zeros(numel(t),14);
                u(:,psIdx) = V;
                TS = timeseries(u,t); %#ok<NASGU>
                TS.Name = 'InputSignals';
                assignin('base','TS',TS);

                testCase.verifyWarningFree( ...
                    @() evalin('base', ...
                        "out = sim('simulator_plasmaless_model_simulink.slx');"), ...
                    sprintf('Simulation raised a warning for coil %s.', coilName{1}));

                out = evalin('base','out');
                Ik = out.simout.Data(:,coilIdx);
                time = out.simout.Time;

                % --- check 1: current at tend exceeds the naive,
                % uncoupled estimate Imin = (V/Lkk)*tend ---
                Imin = (V/Lkk)*tend;
                testCase.verifyGreaterThan(Ik(end), Imin, ...
                    sprintf(['Coil %s: I(end) = %.4g is not above the ', ...
                    'uncoupled estimate Imin = %.4g A. Coupling with ', ...
                    'other superconductive circuits was expected to ', ...
                    'reinforce the current.'], coilName{1}, Ik(end), Imin));

                % --- check 2: tail slope obeys the RL step bounds ---
                % (state/matrix indices are those of the circuits)
                otherSC = setdiff(1:11, circIdx);
                m = M(otherSC,circIdx);
                LeffSC = Lkk - M(circIdx,otherSC)*(M(otherSC,otherSC)\m);
                testCase.assertGreaterThan(LeffSC, 0, ...
                    'Superconductive effective inductance must be positive.');

                % Solve using the full matrix, including the resistive VS3
                % circuit and passive structures, without explicitly inverting it.
                ek = zeros(size(M,1),1);
                ek(circIdx) = 1;
                initialResponse = M\ek;
                slope_lower = V/LeffSC;
                slope_upper = V*initialResponse(circIdx);

                if debugme == 1
                    ax = nexttile(debugLayout);
                    % Integrate the slope bounds from zero initial current.
                    elapsed = time - tstart;
                    plot(ax, time, Ik, 'LineWidth', 1.5, ...
                        'DisplayName', 'Simulated current');
                    hold(ax, 'on');
                    plot(ax, time, slope_lower*elapsed, '--', ...
                        'LineWidth', 1.2, 'DisplayName', ...
                        'Lower bound: (V/L_{eff,SC})(t-t_0)');
                    plot(ax, time, slope_upper*elapsed, ':', ...
                        'LineWidth', 1.2, 'DisplayName', ...
                        'Upper bound: V(M^{-1})_{aa}(t-t_0)');
                    hold(ax, 'off');
                    grid(ax, 'on');
                    xlabel(ax, 'Time [s]');
                    ylabel(ax, 'Current [A]');
                    title(ax, sprintf('%s: %g V step', coilName{1}, V));
                    legend(ax, 'Location', 'best');
                    drawnow;
                end

                tailMask = time >= (tstart + 0.8*(tend-tstart));
                p = polyfit(time(tailMask), Ik(tailMask), 1);
                slope = p(1);
                testCase.verifyGreaterThan(slope, V/Lkk, ...
                    sprintf(['Coil %s: tail slope = %.6g A/s must exceed ', ...
                    'the uncoupled slope V/Lkk = %.6g A/s.'], ...
                    coilName{1}, slope, V/Lkk));
                testCase.verifyGreaterThanOrEqual( ...
                    slope, slope_lower - tol_rel*abs(slope_lower), ...
                    sprintf(['Coil %s: tail slope = %.6g A/s is below ', ...
                    'the superconductive bound V/LeffSC = %.6g A/s.'], ...
                    coilName{1}, slope, slope_lower));
                testCase.verifyLessThanOrEqual( ...
                    slope, slope_upper + tol_rel*abs(slope_upper), ...
                    sprintf(['Coil %s: tail slope = %.6g A/s exceeds ', ...
                    'the full-circuit initial slope = %.6g A/s.'], ...
                    coilName{1}, slope, slope_upper));
            end

            if debugme == 1
                % Keep the plot available before test teardown closes figures.
                waitfor(debugFig);
            end
        end


        function test_resistive_coil(testCase)
            % Drive the resistive VS3 circuit (VS3U/VS3L in antiseries) from
            % zero initial currents with the power supply Vps_VS3. At steady
            % state its circuit current is V/Rcc, with Rcc = R_VS3U + R_VS3L,
            % and Ia_VS3L = +Ic_VS3, Ia_VS3U = -Ic_VS3. The other circuits
            % are superconductive: constant, generally nonzero induced currents.
            debugme = 0; % Set to 1 to plot; close the figure to resume.
            V = 5;       % [V]
            tstart = 0;
            tend = 3;    % [s], voltage stays on through the final sample
            dt = 0.01;
            tol_rel = 0.01; % 1% steady-state and tail-variation tolerance
            tol_abs = 1e-6; % [A], for currents close to zero
            psVS3 = 14;     % Vps = {... PF5 VS1 PF6 VS3}
            coilVS3U = 13; coilVS3L = 14; circVS3 = 12;

            assignin('base','tstart',tstart);
            assignin('base','tend',tend);
            assignin('base','dt',dt);
            model = test_plasmaless_model.configureBase(dt);

            testCase.assertEqual(model.circuit_names{circVS3}, 'VS3');
            % circuit resistance of VS3 from the model matrix Rcc
            Rcc = model.R(circVS3,circVS3);
            testCase.assertEqual(Rcc, model.R_coils(coilVS3U,coilVS3U) ...
                + model.R_coils(coilVS3L,coilVS3L), 'RelTol', 1e-12, ...
                'Rcc(VS3) must be the sum of the two coil resistances.');
            testCase.assertTrue(isfinite(Rcc) && Rcc > 0, ...
                'The VS3 circuit must have finite positive resistance.');
            expectedCurrent = V / Rcc;
            drivenCoils = [coilVS3U coilVS3L];
            expectedCoil = [-expectedCurrent expectedCurrent];

            t = (tstart:dt:tend)';
            u = zeros(numel(t),14);
            u(:,psVS3) = V;
            TS = timeseries(u,t);
            TS.Name = 'InputSignals';
            assignin('base','TS',TS);
            testCase.verifyWarningFree( ...
                @() evalin('base', ...
                    "out = sim('simulator_plasmaless_model_simulink.slx');"), ...
                'Simulation raised a warning for the resistive-coil step.');

            out = evalin('base','out');
            time = out.simout.Time;
            currents = out.simout.Data;
            testCase.assertEqual(time(end), tend, 'AbsTol', dt*1e-6, ...
                'Simulation must cover the full 3-second voltage step.');
            tailMask = time >= tstart + 0.8*(tend-tstart);
            testCase.assertGreaterThanOrEqual(nnz(tailMask), 2, ...
                'At least two tail samples are required to check settling.');

            % Positive Vps_VS3 must give a positive VS3L current (and the
            % antiseries VS3U current of opposite sign) at every sample.
            testCase.verifyGreaterThan(currents(end,coilVS3L), 0, ...
                'A positive VS3 voltage must give a positive VS3L current.');
            testCase.verifyEqual(currents(:,coilVS3U), -currents(:,coilVS3L), ...
                'AbsTol', 1e-9, 'VS3U and VS3L must be in antiseries.');

            % Require the entire tail to be near V/Rcc, rather than merely
            % crossing the expected current at the final sample.
            for k = 1:numel(drivenCoils)
                coilIdx = drivenCoils(k);
                maxError = max(abs(currents(tailMask,coilIdx) - expectedCoil(k)));
                testCase.verifyLessThanOrEqual(maxError, ...
                    tol_abs + tol_rel*abs(expectedCoil(k)), ...
                    sprintf(['Coil %d: tail current has not settled to ', ...
                    '%.6g A (V/Rcc) within 1%% by 3 s (max error %.6g A).'], ...
                    coilIdx, expectedCoil(k), maxError));
            end

            % Measure peak-to-peak tail variation relative to each coil's
            % own peak current, with an absolute floor for unexcited coils.
            for coilIdx = 1:12
                tailCurrent = currents(tailMask,coilIdx);
                variation = max(tailCurrent) - min(tailCurrent);
                currentScale = max(abs(currents(:,coilIdx)));
                testCase.verifyLessThanOrEqual(variation, ...
                    tol_abs + tol_rel*currentScale, ...
                    sprintf(['Coil %d: current is not constant after the ', ...
                    'transient (tail variation %.6g A).'], coilIdx, variation));
            end

            if debugme == 1
                debugFig = figure('Name', ...
                    'VS3 circuit currents - close to resume test', ...
                    'NumberTitle', 'off');
                debugLayout = tiledlayout(debugFig, 2, 1);
                for k = 1:numel(drivenCoils)
                    coilIdx = drivenCoils(k);
                    ax = nexttile(debugLayout);
                    plot(ax, time, currents(:,coilIdx), 'LineWidth', 1.5, ...
                        'DisplayName', 'Simulated current');
                    hold(ax, 'on');
                    yline(ax, expectedCoil(k), '--', ...
                        'DisplayName', 'Steady state: V/Rcc');
                    hold(ax, 'off');
                    grid(ax, 'on');
                    xlabel(ax, 'Time [s]');
                    ylabel(ax, 'Current [A]');
                    title(ax, sprintf('Coil %d: Vps_VS3 = %g V step', coilIdx, V));
                    legend(ax, 'Location', 'best');
                end
                drawnow;
                waitfor(debugFig);
            end
        end

        function test_circuit_mapping(testCase)
            % Circuit topology: series/antiseries coils and power-supply mapping.
            topo = plasmaless_circuit_topology();
            testCase.assertEqual(size(topo.T), [14 12]);
            testCase.assertEqual(size(topo.P), [12 14]);

            % coil currents from circuit currents and back
            Ic = (1:12)';
            Ia = topo.T*Ic;
            testCase.verifyEqual(Ia(3), Ia(4), 'CS1U and CS1L in series.');
            testCase.verifyEqual(Ia(13), -Ic(12), 'Ia_VS3U = -Ic_VS3.');
            testCase.verifyEqual(Ia(14), Ic(12), 'Ia_VS3L = +Ic_VS3.');
            testCase.verifyEqual(map_coils_to_circuits(Ia), Ic, 'AbsTol', 1e-12);
            testCase.verifyEqual(map_coils_to_circuits([Ia Ia]'), [Ic Ic]', ...
                'AbsTol', 1e-12);
            badIa = Ia; badIa(4) = badIa(4)+1;
            testCase.verifyError(@() map_coils_to_circuits(badIa), ...
                'plasmaless:seriesConstraint');

            % power-supply voltages -> circuit voltages
            Vps = (1:14)';   % CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3
            Vc = map_powersupplies_to_circuits(Vps);
            expected = [1 2 3+4 5 6 7 8-12 9-12 10+12 11+12 13 14]';
            testCase.verifyEqual(Vc, expected);
            testCase.verifyEqual(map_powersupplies_to_circuits(Vps'), expected');
            testCase.verifyEqual(map_powersupplies_to_circuits([Vps Vps]'), [expected expected]');
        end

    end

    methods (Static, Access = private)
        function model = configureBase(dt)
            % Build the model from the DD4 machine description
            % (load_plasmaless_md + configure_plasmaless_model) and put
            % model, x0 and A,B,C,D in the base workspace (needed by the
            % Simulink model and by evalin-based calls), as in
            % run_plasmaless_model.m.
            [em_coupling,pf_active,pf_passive] = load_plasmaless_md();
            [model,x0] = configure_plasmaless_model(em_coupling,pf_active,pf_passive,dt);
            assignin('base','model',model);
            assignin('base','x0',x0);
            assignin('base','A',model.A);
            assignin('base','B',model.B);
            assignin('base','C',model.C);
            assignin('base','D',model.D);
        end

        function [circIdx, Lkk, M, psIdx, coilIdx] = getCoilInfo(name)
            % For a single-coil circuit (CS3U, PF1, ...): its index in the
            % circuit state, its self-inductance Lcc, the circuit matrix M
            % ([Mcc Mce; Mec Mee]) of configure_plasmaless_model, the index of
            % its power supply in Vps and of the coil in the 14 coil outputs.
            [em_coupling,pf_active,pf_passive] = load_plasmaless_md();
            model = configure_plasmaless_model(em_coupling,pf_active,pf_passive,0.01);
            M = model.M;
            circIdx = find(strcmp(model.circuit_names, name), 1);
            psIdx = find(strcmp(model.ps_names, name), 1);
            coilIdx = find(strcmp(model.coil_names, name), 1);
            if isempty(circIdx) || isempty(psIdx) || isempty(coilIdx)
                error('test_plasmaless_model:CoilNotFound', ...
                    'Circuit %s not found in the model.', name);
            end
            Lkk = M(circIdx,circIdx);
        end
    end
end
