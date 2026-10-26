%run_vertical_instability Run and plot an open-loop active-coil voltage step.
%   Run this script after add_plasmaless_paths. The reduced states and
%   outputs are perturbations about the loaded equilibrium. The default
%   step is 0.01 V on the power supply best coupled to the unstable mode.
%   do_plot defaults to 1; set it to 0 to run without creating figures.

clearvars -except do_plot
if ~exist('do_plot','var'), do_plot = 0; end
dt = 0.001;
equilibrium = extract_equilibrium_data();
[model,x0] = configure_vertical_instability_model(dt,equilibrium);
[~,poles,left_modes] = eig(model.Ac);
poles = diag(poles);
pole_tolerance = 1e-8*max(1,max(abs(poles)));
unstable = find(real(poles)>pole_tolerance);
assert(isscalar(unstable),'vertical:unstableModes', ...
    'Expected one unstable eigenmode; found %d.',numel(unstable));
growth_rate = real(poles(unstable));
[modal_input,step_coil] = max(abs(left_modes(:,unstable)'*model.Bc));
assert(modal_input>0,'vertical:unreachableMode', ...
    'The unstable mode is not reachable from the power-supply voltages.');
% Extend the open-loop run to a visible 5 cm displacement. This deliberately
% explores a larger perturbation than the earlier near-equilibrium demo.
step_time = 0.2/growth_rate;
dt = min(dt,1/(100*growth_rate));
model.sysd = c2d(model.sys,dt,'zoh');
model.A = model.sysd.A;
model.B = model.sysd.B;
model.C = model.sysd.C;
model.D = model.sysd.D;
model.dt = dt;
target_displacement = 0.08; % [m]
step_amplitude = 0.01; % [V]
tend = step_time+6/growth_rate;
for attempt = 1:4
    time = (0:dt:tend)';
    Va = zeros(numel(time),model.n_ps);   % power-supply voltages Vps
    Va(time>=step_time,step_coil) = step_amplitude;
    [y,time,x] = lsim(model.sysd,Va,time,x0);
    displacement = abs(y(:,model.ind_zIp)/model.Ip0);
    final_sample = find(displacement>=target_displacement,1);
    if ~isempty(final_sample), break; end
    tend = tend+max(1,log(target_displacement/max(displacement(end),eps)))/growth_rate;
end
assert(~isempty(final_sample),'vertical:stepDuration', ...
    'The step response did not reach the requested displacement.');
time = time(1:final_sample);
Va = Va(1:final_sample,:);
y = y(1:final_sample,:);
x = x(1:final_sample,:);
tend = time(end);
zIp = y(:,model.ind_zIp);
Ia = y(:,model.ind_Ia);
Ie = y(:,model.ind_Ie);
fprintf('Vertical growth rate: %.6g 1/s; %g V step on power supply %s.\n', ...
    growth_rate,step_amplitude,model.ps_names{step_coil});
fprintf('Final vertical displacement: %.6g m at t = %.6g s (linear model).\n', ...
    zIp(end)/model.Ip0,time(end));
if do_plot
    [figures,animation] = plot_vertical_instability(model,time,Va,y,equilibrium);
end
