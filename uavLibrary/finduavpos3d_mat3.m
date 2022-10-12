function [Fmin, Xopt, stepCount] = finduavpos3d_mat3(PosUE, PosBS, ...
                                                    U, funs, stepSizeMeter, ...
                                                    Maps, meterPerPixel, map_x0)
% Version 3: use funs for a set of segment specific objective functions
% Version 2: Matrix version, the input urban city map is a set of matrics

PosUE3 = [PosUE, 0];
PosBS3 = [PosBS, U.Hbs];

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
maxStep = ceil((2.4 * U.K - 1.4) * norm(PosBS(1:2) - PosUE(1:2), 2) / stepSizeMeter * U.Hdrone / U.Hmin);

Fmin = inf;
Xhat3 = [0, 0, 0];

XposArray = zeros(maxStep, 2);
FArray = zeros(maxStep, 1);

Xpos = PosBS;
Stage = 1;
cnt = 0;
Rhos = zeros(1, maxStep); % Record for Stage 1, rho value along the way
ks = zeros(1, maxStep);   % Record for stage 1, segment label along the way
ksearch = 1;
while Stage < 4 && cnt < maxStep
    cnt = cnt + 1;
    
    Xpos0 = Xpos;   % UAV Position
    
    try
    
        % los = IsLosK(PosUE, Xpos0, BldLines, BldHeight, U.Hdrone, BldTypes);
        los = IsLosK_discrete([PosUE, 0], [Xpos0, U.Hdrone], Maps, meterPerPixel, map_x0);
        ksegment = round((1 - los) * (U.K - 1) + 1);   % propagation segment index, k = 1,2,...,K
        [f, ~, Xpos3] = getcostf3d([Xpos0, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, funs);
        FArray(cnt) = f;
        XposArray(cnt, :) = Xpos0;
        if f < Fmin
            Fmin = f;
            Xhat3 = Xpos3;
        end

        if Stage == 1
            % Stage 1: Search on the User-BS axis
            if ksegment == 1 % LOS segment
                Rhos(cnt) = norm(PosUE(1:2) - Xpos0(1:2), 2);
                ks(cnt) = 1;
                Rhos0 = findCriticalPoints3d(Rhos, ks, [Xpos0, U.Hdrone], PosUE3, PosBS3, U, funs);

                rho0 = Rhos0(ksearch);
                theta0 = stepSizeMeter / rho0;
                M = [cos(theta0), -sin(theta0)
                     sin(theta0), cos(theta0)];
                Xpos = PosUE + (rho0 * M * Uvec(:)).';

                Stage = 2;

            elseif norm(PosUE - Xpos0) > 0.5  
                Rhos(cnt) = norm(PosUE(1:2) - Xpos0(1:2), 2);
                ks(cnt) = ksegment;
                searchDirection = (PosUE - PosBS) / norm(PosUE - PosBS);
                Xpos = Xpos0 + searchDirection * min(stepSizeMeter, norm(PosUE - Xpos0) - 0.25);

            else
                % Indoor case
                Stage = 4;
            end

        elseif Stage == 2
            % Stage 2: Search on the right branch
            searchDirection = uavSearchDirection_3d([Xpos, U.Hdrone], ...
                                                   PosUE3, PosBS3, ksegment, ksearch, U, funs);
            %
            Xpos = Xpos0 + searchDirection * stepSizeMeter;

            if norm(searchDirection) < 1e-10
                rho0 = Rhos0(ksearch);
                theta0 = - stepSizeMeter / rho0;
                M = [cos(theta0), -sin(theta0)
                     sin(theta0), cos(theta0)];
                Xpos = PosUE + (rho0 * M * Uvec(:)).';

                Stage = 3;
            end

        elseif Stage == 3
            % Stage 3: Search on the left branch
            searchDirection = uavSearchDirection_3d([Xpos, U.Hdrone], ...
                                                   PosUE3, PosBS3, ksegment, ksearch, U, funs);
            %
            Xpos = Xpos0 + searchDirection * stepSizeMeter;

            if norm(searchDirection) < 1e-10

                ksearch = ksearch + 1;
                if ksearch < U.K
                    rho0 = Rhos0(ksearch);
                    theta0 = stepSizeMeter / rho0;
                    M = [cos(theta0), -sin(theta0)
                         sin(theta0), cos(theta0)];
                    Xpos = PosUE + (rho0 * M * Uvec(:)).';
                    Stage = 2;
                else
                    Stage = 4; % Algorithm terminates
                end
            end

        else
            % The entire algorithm terminates
        end

    catch
        error('finduavpos3d_mat3: expected algorithm behavior!');
        debug_here = 1;
    end
    % figure(4), hold on, plot3(Xpos0(1), Xpos0(2), U.Hdrone, 'bx'); hold off
end
Xopt = Xhat3;
stepCount = cnt;
