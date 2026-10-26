% Plot voltage/current histories and animate the magnetic flux.
%
% Optimized figure(2) animation:
%   - coil rectangles (r_a,z_a,width,height) are created ONCE as patch
%     objects and updated in-place (FaceVertexCData) every frame instead
%     of being re-created, colored blue->red, amplitude clamped to
%     [-20 20] on a LINEAR scale, using out.simout.Data(indt,1:12)
%   - passive structure scatter is created ONCE and updated in-place
%     (CData only) every frame using a clamped linear scale
%   - contour is updated in-place (ZData) instead of being re-plotted
%   - no clf() inside the animation loop -> much smoother animation

% after running simulation, plot evolution of magnetic flux
% on a reduced timebase
Nframes = 500;
t_subsample = linspace(tstart,tend,Nframes);

% passive elements location
[em_coupling,pf_active,pf_passive,wall] = load_plasmaless_md();
[Mxi,RR,ZZ] = plasmaless_grid_from_ids(em_coupling);
% wall outlines (limiter unit 1 = first wall, unit 2 = divertor)
r_lim = wall.description_2d{1}.limiter.unit{1}.outline.r;
z_lim = wall.description_2d{1}.limiter.unit{1}.outline.z;
r_div = flip(wall.description_2d{1}.limiter.unit{2}.outline.r);
z_div = flip(wall.description_2d{1}.limiter.unit{2}.outline.z);

r_e = zeros(102,1);
z_e = zeros(102,1);
for ii = 3:102
    r_e(ii) = pf_passive.loop{ii}.element{1}.geometry.oblique.r;
    z_e(ii) = pf_passive.loop{ii}.element{1}.geometry.oblique.z;
end
r_e(1:2) = []; % empty element 2
z_e(1:2) = [];

% plot labels (from the DD4 machine description)
coil_names = cell(1,numel(pf_active.coil));
for ii = 1:numel(coil_names)
    key = strtrim(pf_active.coil{ii}.name);
    tok = regexp(key,'\(([^()]*)\)','tokens','once');
    if ~isempty(tok), key = strtrim(tok{1}); end
    coil_names{ii} = key;
end
loop_names = cellfun(@(l) strtrim(l.name), pf_passive.loop(:)', 'UniformOutput', false);
% inputs are the 14 power-supply voltages Vps (VS1, VS3 included), outputs the
% 14 coil currents Ia (all coils) and the passive loop currents Ie
topo = plasmaless_circuit_topology(coil_names);
in_label = strcat('Vps',topo.ps_names);
out_label = [strcat('Ia',coil_names) strcat('Ie',loop_names)];

% active elements location
Ncoils = 12;
r_a = zeros(Ncoils,1);
z_a = zeros(Ncoils,1);
width = zeros(Ncoils,1);
height = zeros(Ncoils,1);
for ii = 1:Ncoils
    r_a(ii) = pf_active.coil{ii}.element{1}.geometry.rectangle.r;
    z_a(ii) = pf_active.coil{ii}.element{1}.geometry.rectangle.z;
    width(ii) = pf_active.coil{ii}.element{1}.geometry.rectangle.width;
    height(ii) = pf_active.coil{ii}.element{1}.geometry.rectangle.height;
end

% prepare blue-white-red colormap
n = 256;
half = n/2;
blueToWhite = [linspace(0,1,half)', ...
               linspace(0,1,half)', ...
               ones(half,1)];
whiteToRed = [ones(half,1), ...
              linspace(1,0,half)', ...
              linspace(1,0,half)'];
cmap = [blueToWhite; whiteToRed];

% recover signals to plot
t = out.simout.Time;
y = squeeze(out.simout.Data);
u = TS.Data;
if size(y,1) ~= numel(t)
    y = y.';
end

%% actual plot
%% figure 1: voltages and currents history
fig1 = figure(1); clf
% Use separate windows and screen-relative outer bounds, including borders.
set(fig1, 'WindowStyle', 'normal', 'WindowState', 'normal', ...
    'Units', 'normalized', 'OuterPosition', [0.02 0.08 0.47 0.84]);

ax1 = subplot(121);
ax1.XTickLabel = [];
ax1.YTickLabel = [];
pos = ax1.Position;
N = 14;
h = pos(4)/N;
ax = gobjects(N,1);
for k = 1:N
    ax(k) = axes('Position', ...
        [pos(1), pos(2)+(N-k)*h, pos(3), h]);
    plot(ax(k), t, u(:,k), 'LineWidth', 1);
    xlim(ax(k), [t(1), t(end)]);
    ylim([-18 18]);
    ylabel(ax(k), in_label(k));
    if k < N
        ax(k).XTickLabel = [];
    end
    if k ==1
        title('Input voltage Vps [V] (power supplies)')
    end
end
xlabel(ax(14),'Time (s)');

ax2 = subplot(122);
ax2.XTickLabel = [];
ax2.YTickLabel = [];
pos = ax2.Position;
N = 14;
h = pos(4)/N;
ax = gobjects(N,1);
for k = 1:N
    ax(k) = axes('Position', ...
        [pos(1), pos(2)+(N-k)*h, pos(3), h]);

    plot(ax(k), t, y(:,k), 'LineWidth', 1);
    xlim(ax(k), [t(1), t(end)]);
    if k<13
        ylim([-16 16]);
    else
        ylim([-1550 1550]);
    end
    ylabel(ax(k), out_label(k));
    if k < N
        ax(k).XTickLabel = [];
    end

    if k ==1
        title('Output current Ia [A]')
    end
end
xlabel(ax(14),'Time (s)');

%% figure 2: dynamical equilibrium (optimized animation)
fig2 = figure(2); clf
set(fig2, 'WindowStyle', 'normal', 'WindowState', 'normal', ...
    'Units', 'normalized', 'OuterPosition', [0.51 0.08 0.47 0.84]);
ax2d = axes(fig2);
hold(ax2d,'on');

% static geometry, plotted once
plot(ax2d, r_lim, z_lim, 'k');
plot(ax2d, r_div, z_div, 'k');
axis(ax2d,'equal');
xlabel(ax2d,'R [m]')
ylabel(ax2d,'Z [m]')
colormap(ax2d,cmap);

% contour handle, created once, ZData updated every frame
% Use an explicit level list (rather than an automatically-chosen count)
% so that MATLAB does not collapse to only a handful of visible lines
% when the psi range/spacing at a given frame looks "flat" to its
% auto-leveling heuristic.
[~,indt0] = min(abs(out.simout.Time-t_subsample(2)));
psi0 = reshape(Mxi*[out.simout.Data(indt0,:)]',size(RR,1),size(RR,2));

Nlevels = 60;
% Cap the level range using only the first 12 (PF) coils contribution to
% psi, plus a margin (delta). Coils 13/14 (e.g. CS/plasma-like elements)
% induce much larger flux and would otherwise saturate the level range,
% hiding the finer structure induced by the first 12 coils.
psi_pf12 = Mxi(:,1:12)*out.simout.Data(:,1:12)';   % [Ngrid x Ntime]
psi_min0 = min(psi_pf12(:));
psi_max0 = max(psi_pf12(:));
delta = 0.1*(psi_max0 - psi_min0);
psi_min = psi_min0 - delta;
psi_max = psi_max0 + delta;
levelList = linspace(psi_min, psi_max, Nlevels);

[~,hcontour] = contour(ax2d,RR,ZZ,psi0,levelList,'k','LineWidth',0.75);

% coil rectangles (r_a,z_a,width,height) built ONCE as a single patch
% object using NaN-separated faces for efficiency, colored by linear
% amplitude clamped to [-20 20]
coilVerts = zeros(4*Ncoils,2);
coilFaces = nan(Ncoils,4);
for j = 1:Ncoils
    idx = (4*(j-1)+1):(4*j);
    coilVerts(idx,1) = [r_a(j)-width(j)/2, r_a(j)+width(j)/2, ...
                         r_a(j)+width(j)/2, r_a(j)-width(j)/2];
    coilVerts(idx,2) = [z_a(j)-height(j)/2, z_a(j)-height(j)/2, ...
                         z_a(j)+height(j)/2, z_a(j)+height(j)/2];
    coilFaces(j,:) = idx;
end
hp = patch(ax2d, 'Faces',coilFaces, 'Vertices',coilVerts, ...
    'FaceVertexCData', zeros(Ncoils,1), ...
    'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1);

% passive structure scatter, created ONCE, CData updated every frame
% (linear scale)
hs = scatter(ax2d, r_e, z_e, 5, zeros(size(r_e)), 'filled');

cmax_coil = 20; % linear scale clamp requested: [-20 20]
% Saturate the passive-current color scale using only the first half of
% the simulation, to avoid the much higher values observed towards the
% end dominating the color range.
Nhalf = ceil(size(out.simout.Data,1)/2);
cmax_passive = max(max(abs(out.simout.Data(1:Nhalf,17:end))));

clim(ax2d, [-cmax_passive cmax_passive]);
cb = colorbar(ax2d);
cb.Label.String = 'Ie [A]';

titleHandle = title(ax2d, sprintf('t = %2.3f s',t_subsample(1)));

for tt = 1:Nframes
    [~,indt] = min(abs(out.simout.Time-t_subsample(tt)));

    % update contour (magnetic flux)
    psi = reshape(Mxi*[out.simout.Data(indt,:)]',size(RR,1),size(RR,2));
    if indt>1  && range(psi(:)) > eps(max(abs(psi(:))))
        hcontour.ZData = psi;
        hcontour.Visible = 'on';
    else
        hcontour.Visible = 'off';
    end

    % update coil rectangles amplitude (linear, clamped to [-20 20])
    coilAmp = out.simout.Data(indt,1:12);
    coilAmp = max(min(coilAmp,cmax_coil),-cmax_coil);
    hp.FaceVertexCData = coilAmp(:);

    % update passive structure scatter (linear scale, clamped to
    % +/- cmax_passive computed from the first half of the simulation)
    passiveAmp = out.simout.Data(indt,17:end);
    passiveAmp = max(min(passiveAmp,cmax_passive),-cmax_passive);
    hs.CData = passiveAmp(:);

    titleHandle.String = sprintf('t = %2.3f s',t_subsample(tt));
    box on
    drawnow limitrate
    pause(0.001)
end

