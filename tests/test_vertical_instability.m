classdef test_vertical_instability < matlab.unittest.TestCase
    %test_vertical_instability Check the reduced vertical RZIp model.

    methods (TestClassSetup)
        function addProjectPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'imas_model')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'vertical_instability')));
        end
    end

    methods (Test)
        function test_circuit_states_and_power_supply_inputs(testCase)
            % ITER model: 12 circuit + 102 passive states, 14 power-supply
            % inputs, outputs zIp + 14 coil currents + 102 passive currents.
            model = configure_vertical_instability_model(0.001);
            testCase.verifyEqual(size(model.A),[114 114]);
            testCase.verifyEqual(size(model.B),[114 14]);
            testCase.verifyEqual(size(model.C),[117 114]);
            testCase.verifyEqual(size(model.D),[117 14]);
            testCase.verifyEqual(numel(model.circuit_names),12);
            % warm-start equilibrium currents respect series/antiseries
            Ia0 = model.Ia0(1:14);
            testCase.verifyEqual(Ia0(3),Ia0(4),'AbsTol',1e-6);
            testCase.verifyEqual(Ia0(13),-Ia0(14),'AbsTol',1e-6);
            testCase.verifyEqual(model.Ic0(12),Ia0(14),'AbsTol',1e-6);
            % outputs: Ia_CS1U = Ia_CS1L = Ic_CS1, Ia_VS3U = -Ia_VS3L
            Cia = model.C(model.ind_Ia,:);
            testCase.verifyEqual(Cia(3,:),Cia(4,:));
            testCase.verifyEqual(Cia(13,:),-Cia(14,:));
            % Vps_VS1 (input 12) = +PF2 +PF3 -PF4 -PF5 supply voltages (inputs 8..11)
            Bc = model.Bc;
            testCase.verifyEqual(Bc(:,12), ...
                Bc(:,8)+Bc(:,9)-Bc(:,10)-Bc(:,11), ...
                'AbsTol',1e-12*max(abs(Bc(:))), ...
                'Vps_VS1 must act as +PF2 +PF3 -PF4 -PF5.');
            testCase.verifyEqual(model.topology.P(7:10,12),[1;1;-1;-1]);
        end

        function test_single_unstable_eigenmode(testCase)
            model = configure_vertical_instability_model(0.001);
            continuous_poles = eig(model.Ac);
            discrete_poles = eig(model.A);
            tolerance = 1e-8*max(1,max(abs(continuous_poles)));
            testCase.verifyEqual(nnz(real(continuous_poles)>tolerance),1, ...
                'The continuous matrix must have exactly one positive real eigenvalue.');
            testCase.verifyEqual(nnz(abs(discrete_poles)>1+1e-8),1, ...
                'The sampled model must have exactly one growing mode.');
            testCase.verifyTrue(all(isfinite(continuous_poles)));
            testCase.verifyEqual(imag(continuous_poles( ...
                real(continuous_poles)>tolerance)),0,'AbsTol',tolerance);
        end

        function test_force_balance_and_voltage_step(testCase)
            [eq,em,pa,pe] = test_vertical_instability.synthetic_equilibrium();
            [model,x0] = configure_vertical_instability_model(0.01,eq,em,pa,pe);
            time = (0:0.01:10)';
            voltage = ones(numel(time),1);
            [y,~,x] = lsim(model.sysd,voltage,time,x0);
            force_residual = model.M21*x'+model.M22*y(:,1)';
            testCase.verifyEqual(nnz(real(eig(model.Ac))>1e-8),1, ...
                sprintf('M12=%s, M21=%s, M22=%g, poles=%s', ...
                mat2str(model.M12),mat2str(model.M21),model.M22, ...
                mat2str(eig(model.Ac))));
            testCase.verifyLessThan(max(abs(force_residual)),1e-10);
            testCase.verifyEqual(y(:,2:end),x,'AbsTol',1e-10);
            testCase.verifyEqual(model.D,zeros(3,1),'AbsTol',1e-14);
            testCase.verifyGreaterThan(abs(y(end,1)), ...
                abs(y(floor(numel(time)/2),1)));
        end
    end

    methods (Static, Access = private)
        function [eq,em,pa,pe] = synthetic_equilibrium()
            % One active and one passive circuit, with an unstable vertical
            % stiffness. This fixture needs no IMAS database or toolkit.
            [rr,zz] = meshgrid(linspace(1,2,41),linspace(-1,1,81));
            density = exp(-((rr-1.5)/0.1).^2-(zz/0.2).^2);
            density([1 end],:) = 0;
            density(:,[1 end]) = 0;
            eq.RR = rr; eq.ZZ = zz;
            eq.Jtor = density/trapz(rr(1,:),trapz(zz(:,1),density,1));
            eq.Ip0 = 1;
            eq.Ic0 = [1;0];
            pa.coil = {struct('name','active','resistance',1)};
            pe.loop = {struct('name','passive','resistance',1)};
            node = 'grid_ggd(1)/space(1)/objects_per_dimension(1)/object';
            maps = {1,0,1,2*zz(:)+zz(:).^2/2,0.5*zz(:)};
            names = {'mutual_active_active','mutual_passive_active', ...
                'mutual_passive_passive','mutual_grid_active', ...
                'mutual_grid_passive'};
            rows = {'pf_active/coil','pf_passive/loop', ...
                'pf_passive/loop',node,node};
            cols = {'pf_active/coil','pf_active/coil', ...
                'pf_passive/loop','pf_active/coil','pf_passive/loop'};
            for k = 1:numel(maps)
                em.coupling_matrix{k} = struct('name',names{k}, ...
                    'data',maps{k},'quantity',struct('index',1), ...
                    'rows_uri',{{[rows{k} '(:)']}}, ...
                    'columns_uri',{{[cols{k} '(:)']}});
            end
            space.coordinates_type = {struct('index',4),struct('index',3)};
            opd.geometry_content.index = 1;
            for k = 1:numel(rr)
                opd.object{k}.geometry = [rr(k);zz(k)];
            end
            space.objects_per_dimension = {opd};
            em.grid_ggd = {struct('space',{{space}})};
        end
    end
end
