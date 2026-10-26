function TS = timeseries_plasmaless_model(tstart,tend,dt)
% timeseries of the 14 power-supply voltages Vps (one +/-A pulse on each supply):
% Vps = {CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3}
% (see map_powersupplies_to_circuits)

% Parameters
N = 14;           % Number of power supplies
A = 10;           % Pulse amplitude
pulseWidth = 1;   % Pulse duration (seconds)
breakTime = 1;    % Break between channels (s)
% dt = 0.01;        % Sampling interval
debugme = 0;

% Time vector
t = (tstart:dt:tend)';

% Initialize 14 channels
u = zeros(length(t), N);

% Generate one shifted pulse per channel
for k = 1:N
    t_start = (k-1)*(pulseWidth + breakTime);
    t_mid = t_start + pulseWidth/2;
    t_end = t_start + pulseWidth;

    % Positive phase
    u(t >= t_start & t < t_mid, k) = A;

    % Negative phase
    u(t >= t_mid & t < t_end, k) = -A;
end
% Create a single multichannel timeseries
TS = timeseries(u, t);
TS.Name = 'InputSignals';

% Plot all channels separately
if debugme
    figure;
    for k = 1:N
        subplot(N,1,k);
        stairs(t, u(:,k), 'LineWidth', 1.2);
        ylabel(['Ch ' num2str(k)]);
        ylim([-1 11]);
        xlim([0 28]);
        grid on;
    end
    xlabel('Time (s)');
end
