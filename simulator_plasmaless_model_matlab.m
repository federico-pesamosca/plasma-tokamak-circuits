function out = simulator_plasmaless_model_matlab(TS,model,x0)
% run the ABCD model for an initial condition x0 and input timeseries TS
% to generate the output structure out, coherent with the simulink output
% of simulator_plasmaless_model_simulink, developed for the plasmaless model

A = model.A;
B = model.B;
C = model.C;
D = model.D;

y = zeros(numel(TS.Time),size(C,1));

for tt = 1:numel(TS.Time)
    if tt == 1
        x = x0;
    else
        x = xnew;
    end
    u = [TS.Data(tt,:)]'; % prepare input
    [y(tt,:),xnew] = plasmaless_timestep_forward(u,x,A,B,C,D);
end
out.tout = TS.Time;
out.simout = timeseries(y,TS.time);
