function [Fmin, Xopt, stepCount] = finduavposStatK_mat3(PosUE, PosBS, ...
                                                       U, funs, stepSizeMeter, ...
                                                       Maps, meterPerPixel, map_x0, ...
                                                       LosStat)
% Version 3.2: The search should be done in 2D rather than as a line-search.
% Version 3.1b: Change all all U.Hdrone to Hsearch, and user U.Hmin as Hsearch                                                   
% [deleted] Version 3.1: Change all U.Hdrone to U.Hmin (use Hmin as the operational elevation) 
% Version 3: use funs for a set of segment specific objective functions
% Version 2.2: The input urban map is a set of matrices
% Version 2.1: K segment case
% Statistical optimization baseline

 
% PosUE3 = [PosUE, 0];
% PosBS3 = [PosBS, U.Hbs];

K = length(Maps) + 1;

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);

Hsearch = U.Hmin;

% Exhaustive algorithm (to be optimized)
L = norm(PosBS(1:2) - PosUE(1:2), 2);
rho_array = L - abs(stepSizeMeter) : - abs(stepSizeMeter): abs(stepSizeMeter);
Nsteps = length(rho_array);
Fvec = zeros(1, Nsteps);
for i = 1:Nsteps
    rho = rho_array(i);
    elev_angle = atan(Hsearch / rho);
    
    [~, I] = min(abs(LosStat.Angles - elev_angle));
    Prob = LosStat.LosFreq(:, I);
    Prob = Prob / sum(Prob);
    
    d1 = sqrt(rho ^ 2 + (Hsearch - U.Huser) ^ 2);
    d0 = sqrt((L - rho)^2 + (Hsearch - U.Hbs)^2);
    % gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
    % gDrone = 0;
    f = 0;
    for k = 1:K
        % gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
        fk = funs{k}(d1, d0);
        f = f + Prob(k) * fk;
    end
    
    % Fvec(i) = funs(gDrone, gBS); 
    Fvec(i) = f; 
    
end
[~, I] = min(Fvec);
rho_opt = rho_array(I);

Xopt = PosUE + (rho_opt * Uvec(:)).';

% Version 3.2: Improve the search by adjusting the heihgt as well.
% We use gradient search
Xopt0 = [Xopt(1:2), Hsearch];
Xopt = Xopt0;
DX = [1, 1, 1] * stepSizeMeter;
drho = stepSizeMeter / 2;
dz = stepSizeMeter / 2;
MAXLOOP = 200; cnt = 0;
F_array = zeros(1, MAXLOOP);
while norm(DX) > stepSizeMeter * 1e-3 && cnt < MAXLOOP
    cnt = cnt + 1;
    
    Xopt0 = Xopt;
    rho = norm(Xopt0(1:2) - PosUE(1:2)); z = Xopt0(3);
    
    % computer F0 = F(rho, z)
    elev_angle = atan(z / rho);
    [~, I] = min(abs(LosStat.Angles - elev_angle));
    Prob = LosStat.LosFreq(:, I); Prob = Prob / sum(Prob);
    d1 = sqrt(rho ^ 2 + (z - U.Huser) ^ 2);
    d0 = sqrt((L - rho)^2 + (z - U.Hbs)^2);
    F0 = 0;
    for k = 1:K
        % gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
        fk = funs{k}(d1, d0);
        F0 = F0 + Prob(k) * fk;
    end
    
    % computer F1 = F(rho + drho, z);
    elev_angle = atan(z / (rho + drho));
    [~, I] = min(abs(LosStat.Angles - elev_angle));
    Prob = LosStat.LosFreq(:, I); Prob = Prob / sum(Prob);
    d1 = sqrt((rho + drho) ^ 2 + (z - U.Huser) ^ 2);
    d0 = sqrt((L - (rho + drho))^2 + (z - U.Hbs)^2);
    F1 = 0;
    for k = 1:K
        % gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
        fk = funs{k}(d1, d0);
        F1 = F1 + Prob(k) * fk;
    end
    
    dF_drho = (F1 - F0) / drho;
    
    % Computer F2 = F(rho, z + dz)
    elev_angle = atan((z + dz) / rho);
    [~, I] = min(abs(LosStat.Angles - elev_angle));
    Prob = LosStat.LosFreq(:, I); Prob = Prob / sum(Prob);
    d1 = sqrt(rho ^ 2 + ((z + dz) - U.Huser) ^ 2);
    d0 = sqrt((L - rho)^2 + ((z + dz) - U.Hbs)^2);
    F2 = 0;
    for k = 1:K
        % gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
        fk = funs{k}(d1, d0);
        F2= F2 + Prob(k) * fk;
    end
    
    dF_dz = (F2 - F0) / dz;
    
    D = - [dF_drho, dF_dz];
    if norm(D) > stepSizeMeter
        D = D / norm(D) * stepSizeMeter;
    end
    
    % Compute new F and new Xopt
    dF = 1; st = 2;
    while dF > 0 
        st = st * 0.5;
        if st < 0.5^10
            st = 0;
        end
        rho1 = rho + D(1) * st; z1 = max(min(z + D(2) * st, U.Hmax), U.Hmin);
        elev_angle = atan(z1 / rho1);
        [~, I] = min(abs(LosStat.Angles - elev_angle));
        Prob = LosStat.LosFreq(:, I); Prob = Prob / sum(Prob);
        d1 = sqrt(rho1 ^ 2 + (z1 - U.Huser) ^ 2);
        d0 = sqrt((L - rho1)^2 + (z1 - U.Hbs)^2);
        F_array(cnt) = 0;
        for k = 1:K
            % gDrone = gDrone + Prob(k) * 10 .^ ((log10(d1) * U.Alpha(k) + U.Beta(k)) / 10) / U.Noise;
            fk = funs{k}(d1, d0);
            F_array(cnt) = F_array(cnt) + Prob(k) * fk;
        end
        Xopt = [PosUE + (rho1 * Uvec(:)).', z1];
        DX = Xopt  - Xopt0;
        dF = F_array(cnt) - F0 - 1e-10;
    end
end


% if cnt > 2
%     figure(10582), plot(F_array(1:cnt));
%     debug_here = 1;
% end


% Evaluation on actual environment
los = IsLosK_discrete([PosUE, 0], Xopt, Maps, meterPerPixel, map_x0);
ks = round((1 - los) * (K - 1) + 1);

rho_opt = norm(Xopt(1:2) - PosUE(1:2)); z = Xopt(3);

d1 = sqrt(rho_opt ^ 2 + (z - U.Huser) ^ 2);
% gDrone = 10 .^ ((log10(d1) * U.Alpha(ks) + U.Beta(ks)) / 10) / U.Noise;

d0 = sqrt((L - rho_opt)^2 + (z - U.Hbs)^2);
% gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;

% Fmin = funs(gDrone, gBS); 
Fmin = funs{ks}(d1, d0);

stepCount = Nsteps;
