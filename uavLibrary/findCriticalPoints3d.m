function Rhos0 = findCriticalPoints3d(Rhos, ks, PosUAV, PosUE, PosBS, U, funs)
% Version 3: use funs for a set of segment specific objective functions
% Version K, for K segment case

    L = norm(PosUE(1:2) - PosBS(1:2), 2);
    rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
    Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
    Posdelta = PosUAV(1:2) - PosUE(1:2);

    gamma= Posdelta(:).' * Uvec(:) / rho;
    theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

    h = PosUAV(3);
    Hb = PosBS(3);

    if abs(theta) > 1e-6
        error('theta should be zero!');
    end

    Rhos0 = zeros(1, U.K);
    
    for k = 1:U.K
        I = find(ks == k);
        if isempty(I)
            Rhos0(k) = Rhos0(k - 1);
            continue
        end
        rhomin = min(Rhos(I));
        rhomax = max(Rhos(I));
        if k == 1
            rhomin = 0;
        elseif k == U.K
            rhomax = L;
        end 
        los = k; 
            
        delta = 1e-4;
        tol = L * 1e-7;
        cnt = 0;
        MAXLOOP = 1000;
        while rhomax - rhomin > tol && cnt < MAXLOOP
            cnt = cnt +1 ;
            rho = (rhomax + rhomin) / 2;

            F_plus = F(rho + delta, theta, L, h, Hb, U, funs, los);
            F_minus = F(rho - delta, theta, L, h, Hb, U, funs, los);
            dF_drho = (F_plus - F_minus) / (2 * delta);

            if dF_drho >= 0
                rhomax = rho;
            else
                rhomin = rho;
            end
        end

        if cnt >= MAXLOOP
            error('rho did not converge!');
            rho0 = 0;
        else
            rho0 = (rhomax + rhomin) / 2;
        end   
        Rhos0(k) = rho0;
        
    end

end

function f = F(rho, theta, L, h, Hb, U, funs, los)
    % compute the function F(\rho,\theta) in the paper

    lmin = sqrt(rho^2 + h^2) * U.Hmin / h;
    lmax0 = sqrt(rho^2 + h^2) * U.Hmax / h;
    lmax = min(lmax0, (rho * L * cos(theta) + h * Hb) / sqrt(rho^2 + h^2));
    tol = (lmax - lmin) * 1e-8;
    cnt = 0;
    MAXLOOP = 1000;
    if lmax - lmin <= tol
        l = lmin;
    end
    while lmax - lmin > tol && cnt < MAXLOOP
        cnt  = cnt + 1;

        l = (lmax + lmin) / 2;
        dl = (lmax - lmin) * 1e-6;
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
    f = fl(lopt, rho, theta, L, h, Hb, U, funs, los);
    
end

function f = fl(l, rho, theta, L, h, Hb, U, funs, los)
    k = los;
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    % gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    % gainu = 10 ^ ((U.Alpha(k) * log10(du2u) + U.Beta(k)) / 10) / U.Noise;
    
    f = funs{k}(du2u, du2b);
    % f = funs(gainu, gainb);
end
