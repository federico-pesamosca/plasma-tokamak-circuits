function [figures,animation] = plot_vertical_instability(model,time,Va,y,equilibrium,Nframes)
%plot_vertical_instability Plot histories and animate ITER flux and LCFS.
%   plot_vertical_instability(MODEL,TIME,VA,Y,EQ) creates the response
%   figure and an ITER cross-section animation with conductor-current
%   perturbations, wall, flux contours, and a tracked X-point/LCFS.
%
%   [FIGURES,ANIMATION] = plot_vertical_instability(...,NFRAMES) controls
%   the number of frames (default 120) and returns figure handles and
%   X-point diagnostics. EQ.psi0 and EQ.psix must use the same full-flux
%   sign as the conductor maps. TIME, VA and Y are lsim results; VA holds
%   the n_ps power-supply voltages, outputs Y are ordered [zIp;Ia;Ie] with
%   the 14 coil currents Ia (all coils are plotted). Currents are perturbations.
%
%   Each frame includes every conductor current and the first-order
%   rigid plasma shift: delta_psi_p = -delta_z*dpsi_p0/dZ, where
%   psi_p0 = psi0-Mxi*Ia0 (Ia0: single-coil equilibrium currents, the flux
%   maps are per coil, MODEL.Ia0). This is consistent with delta_J =
%   -delta_z*dJ/dZ; it does not solve a nonlinear equilibrium.
%
%   See also run_vertical_instability, extract_equilibrium_data

if nargin < 6 || isempty(Nframes), Nframes = 120; end
validateattributes(Nframes,{'numeric'},{'scalar','integer','positive'});
time = time(:);
n = model.n_coils+model.n_loops;
if size(y,1) ~= numel(time) || size(y,2) ~= n+1 ...
        || ~isequal(size(Va),[numel(time),model.n_ps])
    error('vertical:plotDimensions','TIME, VA and Y dimensions disagree.');
end
RR = equilibrium.RR;
ZZ = equilibrium.ZZ;
psi0 = equilibrium.psi0;
psix0 = equilibrium.psix;
if ~isequal(size(psi0),size(RR),size(ZZ)) ...
        || any(~isfinite(psi0(:))) || ~isscalar(psix0) || ~isfinite(psix0)
    error('vertical:plotFlux','Supply finite psi0 and psix on the EQ mesh.');
end
check_flux_sign(RR,ZZ,psi0,psix0,equilibrium);
[em,pf_active,pf_passive,wall] = load_plasmaless_md();
[Mxi,gridR,gridZ] = plasmaless_grid_from_ids(em);
if size(Mxi,2) ~= n
    error('vertical:plotConductors','Flux maps and model conductor counts differ.');
end
% The shipped conductor maps omit the first R/Z rows of the EQ grid.
% Restrict contours to their common domain rather than extrapolate flux.
r_nodes = RR(1,:); z_nodes = ZZ(:,1);
keepR = r_nodes>=min(gridR(:)) & r_nodes<=max(gridR(:));
keepZ = z_nodes>=min(gridZ(:)) & z_nodes<=max(gridZ(:));
RR = RR(keepZ,keepR);
ZZ = ZZ(keepZ,keepR);
psi0 = psi0(keepZ,keepR);
if min(size(RR))<3
    error('vertical:plotCoverage','Flux and equilibrium meshes do not overlap.');
end
% Use the equilibrium mesh for both the baseline and conductor flux.
if ~isequal(gridR,RR) || ~isequal(gridZ,ZZ)
    interpolated = zeros(numel(RR),n);
    for k = 1:n
        interpolated(:,k) = reshape(interp2(gridR,gridZ, ...
            reshape(Mxi(:,k),size(gridR)),RR,ZZ,'linear',NaN),[],1);
    end
    Mxi = interpolated;
end
if any(~isfinite(Mxi(:)))
    error('vertical:plotCoverage','Conductor flux maps must cover the EQ mesh.');
end
Ic = y(:,[model.ind_Ia model.ind_Ie]);
zIp = y(:,model.ind_zIp);
% Subtract the equilibrium vacuum flux to isolate the plasma self-flux.
% A vertical translation of J translates its flux by the same amount.
% Single-coil equilibrium currents [Ia0;Ie0] (the flux maps are per coil).
psi_p0 = psi0-reshape(Mxi*model.Ia0,size(RR));
[~,dpsi_p0_dZ] = gradient(psi_p0,RR(1,:),ZZ(:,1));
frame_indices = unique(round(linspace(1,numel(time),min(Nframes,numel(time)))));
Nframes = numel(frame_indices);

% Locate the initial saddle near the supplied separatrix flux.
[xpoint,psix] = recompute_xpoint_flux(RR,ZZ,psi0,[],psix0);
if abs(psix-psix0) > 0.01*(max(psi0(:))-min(psi0(:)))
    error('vertical:initialXpoint','Initial saddle flux disagrees with EQ.psix.');
end

figures = gobjects(2,1);
figures(1) = figure('Name','Open-loop vertical instability', ...
    'NumberTitle','off','WindowStyle','normal','Units','normalized', ...
    'OuterPosition',[0.02 0.08 0.47 0.84]);
layout = tiledlayout(figures(1),3,1);
ax = nexttile(layout);
plot(ax,time,zIp,'LineWidth',1.5);
ylabel(ax,'z I_{p0} [A m]'); grid(ax,'on');
title(ax,sprintf('Open-loop voltage step; growth rate %.3g s^{-1}', ...
    max(real(eig(model.Ac)))));
ax = nexttile(layout);
plot(ax,time,Ic(:,1:model.n_coils));
ylabel(ax,'\Delta I_a [A]'); grid(ax,'on');
ax = nexttile(layout);
plot(ax,time,Ic(:,model.n_coils+1:end));
ylabel(ax,'\Delta I_e [A]'); xlabel(ax,'Time [s]'); grid(ax,'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 2: animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures(2) = figure('Name','ITER vertical instability: flux and LCFS', ...
    'NumberTitle','off','WindowStyle','normal','Units','normalized', ...
    'OuterPosition',[0.51 0.08 0.47 0.84]);
ax = axes(figures(2)); hold(ax,'on');
for k = 1:numel(wall.description_2d{1}.limiter.unit)
    outline = wall.description_2d{1}.limiter.unit{k}.outline;
    plot(ax,outline.r,outline.z,'k','LineWidth',1.2,'HandleVisibility','off');
end
% All active-coil elements, including annular VS conductors.
coil_patches = gobjects(0,1);
coil_indices = zeros(0,1);
for k = 1:model.n_coils
    for j = 1:numel(pf_active.coil{k}.element)
        [vertices,~] = conductor_outline(pf_active.coil{k}.element{j}.geometry);
        if isempty(vertices), continue; end
        coil_patches(end+1,1) = patch(ax,vertices(:,1),vertices(:,2),0, ...
            'FaceColor','flat','EdgeColor','k','HandleVisibility','off'); %#ok<AGROW>
        coil_indices(end+1,1) = k; %#ok<AGROW>
    end
end
passive_centres = zeros(0,2);
passive_indices = zeros(0,1);
for k = 1:model.n_loops
    for j = 1:numel(pf_passive.loop{k}.element)
        [~,centre] = conductor_outline(pf_passive.loop{k}.element{j}.geometry);
        if isempty(centre), continue; end
        passive_centres(end+1,:) = centre; %#ok<AGROW>
        passive_indices(end+1,1) = model.n_coils+k; %#ok<AGROW>
    end
end
passive_scatter = scatter(ax,passive_centres(:,1),passive_centres(:,2), ...
    14,zeros(size(passive_indices)),'filled','HandleVisibility','off');
half = 128;
current_colormap = [linspace(0,1,half)' linspace(0,1,half)' ones(half,1); ...
    ones(half,1) linspace(1,0,half)' linspace(1,0,half)'];
colormap(ax,current_colormap);
% Separate fixed scales avoid washing out coil colours with vessel peaks.
% Saturate at half the group peak to emphasize smaller current changes.
coil_limit = max(0.5*max(abs(Ic(frame_indices,1:model.n_coils)),[],'all'),eps);
passive_limit = max(0.5*max(abs(Ic(frame_indices,model.n_coils+1:end)),[],'all'),eps);
clim(ax,[-passive_limit passive_limit]);
cb = colorbar(ax); cb.Label.String = '\Delta I_e [A] (clipped)';
text(ax,0.02,0.03,sprintf('Coils: sqrt intensity, saturation \\pm %.3g A',coil_limit), ...
    'Units','normalized','BackgroundColor','white','FontSize',8);
levels = linspace(min(psi0(:)),max(psi0(:)),20);
[~,flux_contour] = contour(ax,RR,ZZ,psi0,levels,'Color',[0.6 0.6 0.6], ...
    'LineWidth',0.5,'HandleVisibility','off');
[~,initial_lcfs] = contour(ax,RR,ZZ,psi0,[psix psix], ...
    'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',1);
[~,lcfs] = contour(ax,RR,ZZ,psi0,[psix psix], ...
    'Color',[0.85 0.1 0.6],'LineWidth',2);
xpoint_marker = plot(ax,xpoint(1),xpoint(2),'x','Color',[0.85 0.1 0.6], ...
    'MarkerSize',9,'LineWidth',2,'HandleVisibility','off');
% Trapezoidal quadrature gives a current-weighted physical centroid.
r_weights = node_weights(RR(1,:));
z_weights = node_weights(ZZ(:,1));
weights = abs(model.Jtor(keepZ,keepR)).*(z_weights*r_weights');
centroid = [sum(RR.*weights,'all') sum(ZZ.*weights,'all')]/sum(weights,'all');
centroid_path = [repmat(centroid(1),numel(time),1),centroid(2)+zIp/model.Ip0];
trajectory = plot(ax,centroid(1),centroid(2),'-','Color',[0.1 0.4 0.85], ...
    'LineWidth',2);
centroid_marker = plot(ax,centroid(1),centroid(2),'o', ...
    'Color',[0.1 0.4 0.85],'MarkerFaceColor',[0.1 0.4 0.85]);
legend(ax,[initial_lcfs lcfs trajectory], ...
    {'Initial LCFS','Flux-map separatrix','Centroid trajectory'},'Location','best');
axis(ax,'equal'); box(ax,'on');
xlabel(ax,'R [m]'); ylabel(ax,'Z [m]');
title_handle = title(ax,'');
% A 5 cm trajectory is small on the full ITER cross section: show a zoom.
zoom_axes = axes(figures(2),'Position',[0.55 0.17 0.14 0.20]);
hold(zoom_axes,'on'); box(zoom_axes,'on'); grid(zoom_axes,'on');
zoom_trajectory = plot(zoom_axes,centroid(1),centroid(2),'-', ...
    'Color',[0.1 0.4 0.85],'LineWidth',2);
zoom_marker = plot(zoom_axes,centroid(1),centroid(2),'o', ...
    'Color',[0.1 0.4 0.85],'MarkerFaceColor',[0.1 0.4 0.85]);
xlim(zoom_axes,centroid(1)+[-0.025 0.025]);
zoom_z = [min(centroid_path(:,2)) max(centroid_path(:,2))];
zoom_padding = max(0.005,0.15*diff(zoom_z));
ylim(zoom_axes,zoom_z+[-zoom_padding zoom_padding]);
title(zoom_axes,'Centroid path','FontSize',9);
xlabel(zoom_axes,'R [m]'); ylabel(zoom_axes,'Z [m]');
zoom_axes.FontSize = 8;
animation.time = time(frame_indices);
animation.psix = zeros(Nframes,1);
animation.xpoint = zeros(Nframes,2);
animation.gradient_norm = zeros(Nframes,1);
animation.hessian_determinant = zeros(Nframes,1);
animation.z = zIp(frame_indices)/model.Ip0;
animation.initial_psix = psix0;
animation.frame_indices = frame_indices;
animation.centroid = centroid_path(frame_indices,:);
animation.n_active_currents = model.n_coils;
animation.n_passive_currents = model.n_loops;
animation.coil_colour_limit = coil_limit;
animation.passive_colour_limit = passive_limit;
animation.conductor_flux_amplitude = zeros(Nframes,1);
animation.plasma_flux_amplitude = zeros(Nframes,1);
for frame = 1:Nframes
    if ~isgraphics(figures(2)), break; end
    sample = frame_indices(frame);
    conductor_flux = reshape(Mxi*Ic(sample,:)',size(RR));
    plasma_flux = -(zIp(sample)/model.Ip0)*dpsi_p0_dZ;
    psi = psi0+conductor_flux+plasma_flux;
    animation.conductor_flux_amplitude(frame) = max(abs(conductor_flux(:)));
    animation.plasma_flux_amplitude(frame) = max(abs(plasma_flux(:)));
    [xpoint,psix,residual,determinant] = recompute_xpoint_flux( ...
        RR,ZZ,psi,xpoint,psix);
    animation.psix(frame) = psix;
    animation.xpoint(frame,:) = xpoint;
    animation.gradient_norm(frame) = residual;
    animation.hessian_determinant(frame) = determinant;
    flux_contour.ZData = psi;
    lcfs.ZData = psi;
    lcfs.LevelList = [psix psix];
    xpoint_marker.XData = xpoint(1);
    xpoint_marker.YData = xpoint(2);
    centroid_marker.YData = centroid_path(sample,2);
    trajectory.XData = centroid_path(1:sample,1);
    trajectory.YData = centroid_path(1:sample,2);
    zoom_trajectory.XData = trajectory.XData;
    zoom_trajectory.YData = trajectory.YData;
    zoom_marker.YData = centroid_path(sample,2);
    for k = 1:numel(coil_patches)
        intensity = max(-1,min(1,Ic(sample,coil_indices(k))/coil_limit));
        intensity = sign(intensity)*sqrt(abs(intensity));
        colour_index = 1+round((size(current_colormap,1)-1)* ...
            (1+intensity)/2);
        coil_patches(k).FaceColor = current_colormap(colour_index,:);
    end
    passive_scatter.CData = Ic(sample,passive_indices)';
    title_handle.String = sprintf('t = %.3f s; \\Delta z = %.3f mm; \\psi_X = %.5g Wb', ...
        time(sample),1000*zIp(sample)/model.Ip0,psix);
    drawnow limitrate;
    pause(0.001);
end
end

function check_flux_sign(RR,ZZ,psi0,psix,equilibrium)
% Ampere's law verifies that signed total flux matches the signed Ip.
r = RR(1,:); z = ZZ(:,1);
[dpsiR,dpsiZ] = gradient(psi0,r,z);
[dpsiRR,~] = gradient(dpsiR,r,z);
[~,dpsiZZ] = gradient(dpsiZ,r,z);
J = -(dpsiRR-dpsiR./RR+dpsiZZ)./(2*pi*4e-7*pi*RR);
inside = inpolygon(RR,ZZ,equilibrium.boundary.r,equilibrium.boundary.z);
J(~inside) = 0;
Ip_flux = trapz(r,trapz(z,J,1));
if sign(Ip_flux) ~= sign(equilibrium.Ip0)
    error('vertical:fluxSign','psi0 sign does not agree with equilibrium Ip.');
end
boundary_flux = interp2(RR,ZZ,psi0,equilibrium.boundary.r, ...
    equilibrium.boundary.z,'spline');
if max(abs(boundary_flux-psix)) > 0.01*(max(psi0(:))-min(psi0(:)))
    error('vertical:boundaryFlux','psi0 and psix signs/levels are inconsistent.');
end
end

function [point,psix,residual,determinant] = recompute_xpoint_flux(RR,ZZ,psi,seed,target_flux)
% Track a genuine saddle: grad(psi)=0 and det(Hessian)<0.
% The stationary point is found in a smooth interpolant, not at the
% nearest grid node or by evaluating the old boundary flux level.
r = RR(1,:); z = ZZ(:,1);
F = griddedInterpolant({z,r},psi,'spline','none');
h = min([diff(r(:));diff(z(:))])/20;
if ~isempty(seed)
    [point,residual,determinant,found] = refine_saddle(F,seed,r,z,h);
    if found
        psix = F(point(2),point(1));
        return
    end
end
% Initial localization: cells spanning zero in both gradient components.
[gR,gZ] = gradient(psi,r,z);
aR = cat(3,gR(1:end-1,1:end-1),gR(2:end,1:end-1), ...
    gR(1:end-1,2:end),gR(2:end,2:end));
aZ = cat(3,gZ(1:end-1,1:end-1),gZ(2:end,1:end-1), ...
    gZ(1:end-1,2:end),gZ(2:end,2:end));
[iz,ir] = find(min(aR,[],3)<=0 & max(aR,[],3)>=0 ...
    & min(aZ,[],3)<=0 & max(aZ,[],3)>=0);
best = Inf;
point = [NaN NaN]; residual = Inf; determinant = NaN;
for k = 1:numel(ir)
    candidate = [(r(ir(k))+r(ir(k)+1))/2 (z(iz(k))+z(iz(k)+1))/2];
    [candidate,cres,cdet,found] = refine_saddle(F,candidate,r,z,h);
    if ~found, continue; end
    candidate_flux = F(candidate(2),candidate(1));
    score = abs(candidate_flux-target_flux);
    if ~isempty(seed), score = score+norm(candidate-seed); end
    if score < best
        best = score; point = candidate; residual = cres; determinant = cdet;
    end
end
if ~isfinite(best)
    error('vertical:xpoint','No converged flux saddle was found in this frame.');
end
psix = F(point(2),point(1));
end

function [point,residual,determinant,found] = refine_saddle(F,point,r,z,h)
found = false; residual = Inf; determinant = NaN;
for iteration = 1:30
    if point(1)<=r(1)+2*h || point(1)>=r(end)-2*h ...
            || point(2)<=z(1)+2*h || point(2)>=z(end)-2*h
        return
    end
    R = point(1); Z = point(2); centre = F(Z,R);
    plusR = F(Z,R+h); minusR = F(Z,R-h);
    plusZ = F(Z+h,R); minusZ = F(Z-h,R);
    g = [(plusR-minusR);(plusZ-minusZ)]/(2*h);
    mixed = (F(Z+h,R+h)-F(Z+h,R-h)-F(Z-h,R+h)+F(Z-h,R-h))/(4*h^2);
    H = [(plusR-2*centre+minusR)/h^2 mixed; ...
        mixed (plusZ-2*centre+minusZ)/h^2];
    residual = norm(g); determinant = det(H);
    if residual < 1e-7
        found = determinant < 0;
        return
    end
    if any(~isfinite(H(:))) || rcond(H)<1e-12, return; end
    correction = H\g;
    % Limit Newton steps to avoid jumping to another magnetic critical point.
    max_step = 2*min([diff(r(:));diff(z(:))]);
    correction = correction*min(1,max_step/max(norm(correction),eps));
    point = point-correction';
end
end

function [vertices,centre] = conductor_outline(geometry)
vertices = []; centre = [];
switch geometry.geometry_type
    case 2 % Rectangle
        g = geometry.rectangle;
        centre = [g.r g.z];
        vertices = centre+[-g.width -g.height;g.width -g.height; ...
            g.width g.height;-g.width g.height]/2;
    case 3 % Oblique rectangle (scatter uses its centre)
        g = geometry.oblique;
        centre = [g.r g.z];
        vertices = centre;
    case 5 % Annular conductor cross section
        g = geometry.annulus;
        centre = [g.r g.z];
        theta = linspace(0,2*pi,40)';
        vertices = centre+g.radius_outer*[cos(theta) sin(theta)];
    otherwise
        if ~isempty(geometry.outline.r)
            vertices = [geometry.outline.r(:) geometry.outline.z(:)];
            centre = mean(vertices,1);
        end
end
if ~isempty(centre) && (any(~isfinite(centre)) || any(abs(centre)>1e10))
    vertices = []; centre = [];
end
end

function weights = node_weights(nodes)
% Trapezoidal quadrature on a tensor-product mesh.
nodes = nodes(:);
weights = [diff(nodes(1:2));nodes(3:end)-nodes(1:end-2); ...
    diff(nodes(end-1:end))]/2;
end
