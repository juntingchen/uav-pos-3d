function [fmin, lopt, Xopt3] = getcostf3d(PosUAV, PosUE, PosBS, los, U, funs)
% Version 3.1: Activate the UAV height constraint h < Hmax
% Version 3: use funs for a set of segment specific objective functions
% Version K, for K segmenet case, 
% where parameter los = 1, ..., 1/(K-1), 0
%
% Compute the partial optimal cost at a UAV position in 2D. 
% 
% INPUT 
%   PosUAV, PosUE, PosBS, positions in 3D
%   U, configuration structure, containing channel parameters
%   fun, objective function in terms of the channel gains
%
% OUTPUT
%   fmin,   the minimum cost over the ray with phi evalvation angle
%           depending on the UAV-user position
%   lopt,   corresponding coordinate in l.
%   Xopt    Optimal UAV position
%
% Constraint: l >= sqrt(rho^2 + h^2)
% NOTE: Assume user height = 0.

    L = norm(PosUE(1:2) - PosBS(1:2), 2);
    rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
    Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
    Posdelta = PosUAV(1:2) - PosUE(1:2);

    gamma= Posdelta(:).' * Uvec(:) / rho;
    theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

    h = PosUAV(3);
    Hb = PosBS(3);

    lmin = sqrt(rho^2 + h^2) * U.Hmin / h;
    lmax0 = sqrt(rho^2 + h^2) * U.Hmax / h;
    lmax = min(lmax0, (rho * L * cos(theta) + h * Hb) / sqrt(rho^2 + h^2));
    dl = (lmax - lmin) * 1e-6;
    tol = (lmax - lmin) * 1e-7;
    cnt = 0;
    MAXLOOP = 1000;
    if lmax - lmin <= tol
        l = lmin;
    end
    
    while lmax - lmin > tol && cnt < MAXLOOP
        cnt  = cnt + 1;

        l = (lmax + lmin) / 2;
        
        f_plus = fl(l + dl, rho, theta, L, h, Hb, U, funs, los);
        f_minus = fl(l - dl, rho, theta, L, h, Hb, U, funs, los);
        df = (f_plus - f_minus) / (2 * dl);

        if df <= 0
            lmin = l;
        else
            lmax = l;
        end

    end
    if cnt >= MAXLOOP
        error('failed to converge for l_star!');
    end

    lopt = l;
    fmin = fl(lopt, rho, theta, L, h, Hb, U, funs, los);

    % Compute the optimal UAV position in 3D (optimizing variable l only)
    h1 = lopt * h / sqrt(rho^2 + h^2);
    rho1 = lopt * rho / sqrt(rho^2 + h^2);
    
    M = [cos(theta), -sin(theta)
         sin(theta), cos(theta)];
    Xopt = PosUE(1:2) + (rho1 * M * Uvec(:)).';
    Xopt3 = [Xopt, h1];
end

function f = fl(l, rho, theta, L, h, Hb, U, funs, los)
    k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    % gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    % gainu = 10 ^ ((U.Alpha(k) * log10(du2u) + U.Beta(k)) / 10) / U.Noise;

    f = funs{k}(du2u, du2b);
    % f = funs(gainu, gainb);
end