% Plot the results of a run with an external linear model (run_plasmaless_imas.m
% or run_plasmaless_model.m with linear_model_dir set, or a loaded
% plasmaless_*linear_results.mat): reads time, Ia, Ie, Y, coil_names,
% loop_names, y_names from the workspace.
%   figure 1: coil-current deviations (one line per coil)
%   figure 2: passive-conductor currents (one line per passive state)
%   figure 3: observables, one subplot per y_name

%% figure 1: coil currents
figure(1); clf
plot(time, Ia, 'LineWidth', 1);
grid on
xlabel('Time (s)'); ylabel('Coil current deviation (A)');
title('Linear model: coil currents');
legend(coil_names, 'Location', 'eastoutside');

%% figure 2: passive currents
figure(2); clf
plot(time, Ie, 'LineWidth', 1);
grid on
xlabel('Time (s)'); ylabel('Passive current (A)');
title('Linear model: passive conductor currents');
legend(loop_names, 'Location', 'eastoutside', 'Interpreter', 'none');

%% figure 3: observables
n_y = numel(y_names);
if n_y > 0
    figure(3); clf
    n_col = ceil(sqrt(n_y));
    n_row = ceil(n_y/n_col);
    for k = 1:n_y
        subplot(n_row, n_col, k);
        plot(time, Y(:,k), 'LineWidth', 1);
        grid on
        xlabel('Time (s)');
        title(y_names{k}, 'Interpreter', 'none');
    end
end
