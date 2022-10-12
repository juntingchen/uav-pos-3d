function [Fmin, Xopt, stepCount] = finduavposStatK_mat(PosUE, PosBS, ...
                                                       U, fun, stepSizeMeter, ...
                                                       Maps, meterPerPixel, map_x0, ...
                                                       LosStat)
% Version 2.2: The input urban map is a set of matrices
% Version 2.1: K segment case
% Statistical optimization baseline

 
% PosUE3 = [PosUE, 0];
% PosBS3 = [PosBS, U.Hbs];

K = length(Maps) + 1;

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);

% Exhaustive algorithm (to be optimized)
L = norm(PosBS(1:2) - PosUE(1:2), 2);
rho_array = L - abs(stepSizeMeter) : - abs(stepSizeMeter): abs(stepSizeMeter);
Nsteps = length(rho_array);
Fvec = zeros(1, Nsteps);
for i = 1:Nsteps
    rho = rho_array(i);
    elev_angle = atan(U.Hdrone / rho);
    
    [~, I] = min(abs(LosStat.Angles - elev_angle));
    Prob = LosStat.LosFreq(:, I);
    Prob = Prob / sum(Prob);
    
    d1 = sqrt(rho ^ 2 + (U.Hdrone - U.Huser) ^ 2);
    gDrone = 0;
    for k = 1:K
        gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
    end
    
    d0 = sqrt((L - rho)^2 + (U.Hdrone - U.Hbs)^2);
    gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
    
    Fvec(i) = fun(gDrone, gBS); 
end
[~, I] = min(Fvec);
rho_opt = rho_array(I);

Xopt = PosUE + (rho_opt * Uvec(:)).';

% Evaluation on actual environment
los = IsLosK_discrete([PosUE, 0], [Xopt, U.Hdrone], Maps, meterPerPixel, map_x0);
ks = round((1 - los) * (K - 1) + 1);
d1 = sqrt(rho_opt ^ 2 + (U.Hdrone - U.Huser) ^ 2);
gDrone = 10 .^ ((log10(d1) * U.Alpha(ks) + U.Beta(ks)) / 10) / U.Noise;
d0 = sqrt((L - rho_opt)^2 + (U.Hdrone - U.Hbs)^2);
gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
Fmin = fun(gDrone, gBS); 

stepCount = Nsteps;
