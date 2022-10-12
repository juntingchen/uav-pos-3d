function D = uavSearchDirection_3d(PosUAV, PosUE, PosBS, ksegment, ksearch, U, funs)
% Version 3: Use funs for a set of segment specific objective functions.
% Built from uavSearchDreciton3DK.m
% 
% Version K, K segment case
% INPUT 
%   PosUAV, PosUE, PosBS, positions in 3D
%   U, configuration structure, containing channel parameters
%   fun, objective function in terms of the channel gains
%   stepSize, the step size

L = norm(PosUE(1:2) - PosBS(1:2), 2);
rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
Posdelta = PosUAV(1:2) - PosUE(1:2);

gamma= Posdelta(:).' * Uvec(:) / rho;
theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

h = PosUAV(3);
Hb = PosBS(3);

% ---- Numerical Computation of the Partial Derivatives ---- %
delta = 1.5; delta_t = 1.5e-1;

if rho < delta *2 
    D = 0;
    return
end

F_plus = F(rho + delta, theta, L, h, Hb, U, funs, ksearch);
F_minus = F(rho - delta, theta, L, h, Hb, U, funs, ksearch);
dF_drho = (F_plus - F_minus) / (2 * delta);

F_plus = F(rho, theta + delta_t, L, h, Hb, U, funs, ksearch);
F_minus = F(rho, theta - delta_t, L, h, Hb, U, funs, ksearch);
dF_dtheta = (F_plus - F_minus) / (2 * delta_t);

% Cond = rho^2 - L * cos(theta) * h * rho / U.Hmin + h^2 - L * cos(theta) * h^2 * Hb / U.Hmin;
Cond = rho - h / U.Hmin * L * cos(theta);

if dF_drho >= 0 || Cond > 0 % Check stopping criteria
    D = zeros(1, 2);
else
    if ksegment <= ksearch % LOS or virtual LOS
        D = Posdelta / norm(Posdelta, 2);
    else
        if abs(dF_dtheta) < 1e-10
            deltaX = Rot(theta) * Uvec(:);
        else
            deltaX = Rot(theta) * Uvec(:) + rho * dRot(theta) * Uvec(:) ...
                     * ( - dF_drho / dF_dtheta);
        end
        D = deltaX.' / norm(deltaX);
    end
        
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = Rot(ta)
    T = [cos(ta)    -sin(ta)
         sin(ta)    cos(ta)];
end

function T = dRot(ta)
    T = [-sin(ta)    -cos(ta)
         cos(ta)    -sin(ta)];
end



function f = F(rho, theta, L, h, Hb, U, funs, ksearch)
    % compute the function F(\rho,\theta) in the paper

    lmin = sqrt(rho^2 + h^2) * U.Hmin / h;
    lmax0 = sqrt(rho^2 + h^2) * U.Hmax / h;
    lmax = min(lmax0, (rho * L * cos(theta) + h * Hb) / sqrt(rho^2 + h^2));
    dl = (lmax - lmin) * 1e-6;
    tol = max(1, (lmax - lmin) * 1e-8);
    cnt = 0;
    MAXLOOP = 1000;
    if lmax - lmin <= tol
        l = lmin;
    end
    while lmax - lmin > tol && cnt < MAXLOOP
        cnt  = cnt + 1;

        l = (lmax + lmin) / 2;
        
        f_plus = fl(l + dl, rho, theta, L, h, Hb, U, funs, ksearch);
        f_minus = fl(l - dl, rho, theta, L, h, Hb, U, funs, ksearch);
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
    f = fl(lopt, rho, theta, L, h, Hb, U, funs, ksearch);
    
end

function f = fl(l, rho, theta, L, h, Hb, U, funs, ksearch)
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    % gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    % gainu = 10 ^ ((U.Alpha(ksearch) * log10(du2u) + U.Beta(ksearch)) / 10) / U.Noise;

    f = funs{ksearch}(du2u, du2b);
    % f = fun(gainu, gainb);
end