%% Optimal UAV Placement in 3D Space
% Version with Shadow Fading components
% Built from Demo_sim_uav3d.m (verions 3??)
% Date: March 6th. 2020.
%
% Version 3: 
% Modified from sim_capacityK3_discrete.m
%   - Use the segment specific service model, where the objective function
%   is written as a set of objective functions {f_k(d1, d2)}, in which d1
%   and d2 are distances from the UE and BS, respectively. 
%
%   - Use matrix for (discrete) urban city map.
% 
%   - For the journal paper of 3D uav positioning
%
% Date: July 7th, 2018.
%
% Version 2: Build and use discrete elevation map
% (Compare the results to Demo_uavPosSearch3DK.m)
% Date: July 5th, 2018.
%
% Demonstration of hte UAV search trajectory in 3D space under the genernal 
% K segment propagation model. Here, we consider K = 3 segments, simulated as 
% a propagation consequence from K-1 = 2 types of obstacles.
% 
% Author: Junting Chen (juntingc@usc.edu)
% 
% Date: December 20, 2017.
% 
% The propagation parameters are obtained from LTE TR 36.814
%% Basic Configuration
% Here, I simply load the parameters from a pre-generated file.

clear
addpath ../uavLibrary
addpath ../..

DATA = load('../uavDataset/urbanMapSingleUser.mat');
U = DATA.U; 
% PosUE = [342, 237];
% PosUE = [94.0212  739.6973];
% PosUE = [498.8438  504.1995];
% PosUE = [  220.6620  743.5349];
% PosUE = rand(1, 2) * 500 + 100;
% PosUE = [578, 342];
PosUE = [496, 580];
% PosUE = [497, 600];
PosUE =[571, 524];
% PosUE =[566, 527];

PosBS = [150, 770];

% DATA = load('../uavDataset/urbanMapK_1m.mat');
DATA = load('../uavDataset/urbanMap_DC.mat');
Maps = DATA.Maps;
[MapWidth, MapHeight] = size(Maps{1});
Blds = DATA.BldArea;
stepsize = DATA.stepsize;
map_x0 = DATA.x0;
map_x0 = [0, 0];
Xvec = DATA.Xvec;
Yvec = DATA.Yvec;
U.K = length(Maps) + 1; K = U.K;
urbanMap = zeros(MapWidth, MapHeight);
for k = 1:U.K - 1
	urbanMap = urbanMap + Maps{k};
end

% Radio Configuration and initialization
U.A0 = -22; U.B0 = -40; % UAV-BS channel
U.A1 = -22.7; U.B1 = -28; U.S1 = 1; % Propagation segment 1 (LOS)
U.A2 = -28.4; U.B2 = -24; U.S2 = 3; % Propagation segment 2
U.A3 = -36.4; U.B3 = -22; U.S3 = 6; % Propagation segment 3 (NLOS)
% U.K = 3;
U.Alpha = [-22.7, -28.4, -36.5];
U.Beta= [-28, -24, -22];
U.Noise = 1e-10; U.Pb = 2; U.Pd = 2;
U.Hbs = 45;
U.Hmin = 45; % Minimum UAV operation height
U.Hmax = 200; % Maximum UAV operation height 
U.Hdrone = 50;
% Rmap.Alpha = [-20, -40]; Rmap.Beta = [-40, -40]; Rmap.Sigma2 = [3, 8]; Rmap.Pi = [0.3, 0.7];
% Chann.A0 = U.A0; Chann.B0 = U.B0; Chann.Noise = U.Noise;
dcor =  5; % Decorrelation distance in meter: h(x) = exp(- abs(x) / dcor);

stepSizeMeter = 3;

%% User Location 

% PosUE = rand(1, 2) * 500 + 400;
% PosUE = [660.8249  448.3650];
% PosUE = [590, 740]; % Case where optimal solution is in k=2 segment.
% % Random USER location
% inBld = 1;
% while inBld
%     pos = rand(1, 2) * 500 + 400;
%     inBld = 0;
%     for ib = 1:Nbld
%         in = inpoly(pos, Blds{ib, 1});
%         if in
%             inBld = 1;
%             break
%         end
%     end
% end
% PosUE = pos;
%% Objective function 
% Alternatively, you can also create a fun.m file to specify the objective function 
% fun(x, y), where the argument x is the UAV-USER SNR and y is the UAV-BS SNR
% 
% Supported functions: $f(x,y)=\max \{-\log_2(1+P_dx),-\log_2(1+P_b(y)\}$,$f(x,y) 
% = a_1f_1(x)+a_2f_x(y)$


fun = @(x,y) max(-log2(1 + U.Pd * real(x)), -log2(1 + U.Pb * real(y)));
% fun = @(x,y) log(1 ./ (U.Pd * real(x)) + 1 ./ (U.Pb * real(y) * real(x) * U.Pd));

%% Urban Map
uid = round(PosUE / stepsize);
show_map(Xvec, Yvec, urbanMap, Blds, 1);hold on
plot3(PosBS(1), PosBS(2), 50, 'r^', 'linewidth', 2, 'markersize', 9);
plot3(PosUE(1), PosUE(2), 50, 'ro', 'linewidth', 2, 'markersize', 9);hold off


%% Radio Map

Npt = 81;
% Xrange = Xvec(1) + (0:1/(Npt - 1):1) * MapWidth;
% Yrange = Yvec(1) + (0:1/(Npt - 1):1) * MapHeight;

PosMid = (PosUE + PosBS)/2;
WinLen = norm(PosUE - PosBS, 2)*0.5;
Xrange = PosMid(1) + (0:1/(Npt - 1):1) * WinLen;
Yrange = PosMid(2) - WinLen + (0:1/(Npt - 1):1) * WinLen;
Zrange = U.Hmin : stepSizeMeter : U.Hmax;
Nh = length(Zrange);

[Xmat, Ymat] = meshgrid(Xrange, Yrange);
Pmat = zeros(Npt);Cmat = zeros(Npt);Cmat3d = zeros(Npt); % Hmat = zeros(Npt);
Pmat3 = zeros(Npt, Npt, Nh); Cmat3d3 = zeros(Npt, Npt, Nh); % Hmat3 = zeros(Npt, Npt, Nh);
LosMat3 = zeros(Npt, Npt, Nh); SFnoise3 = zeros(Npt, Npt, Nh); Dmat3 = zeros(Npt, Npt, Nh);
for i = 1:Npt
    for j = 1:Npt
        for k = 1:Nh
            Hdrone = Zrange(k);
            
            Upos = [Xrange(i), Yrange(j)];
            du2u = norm([Upos, Hdrone] - [PosUE, 0]);
            % los = IsLosK(PosUE, Upos, BldLines, BldHeight, U.Hdrone, BldTypes);
            los = IsLosK_discrete([PosUE, 0], [Upos, Hdrone], Maps, stepsize, map_x0);
            if abs(los - 1) < 1e-9 % propagation segment 1
                gain_dB = log10(du2u) * U.A1 + U.B1;
                sf_noise = randn * U.S1;
            elseif abs(los - (1-1/(K-1))) < 1e-9     % propagation segment 2
                gain_dB = log10(du2u) * U.A2 + U.B2;
                sf_noise = randn * U.S2;
            else
                gain_dB = log10(du2u) * U.A3 + U.B3;
                sf_noise = randn * U.S3;
            end
            LosMat3(i, j, k) = los; SFnoise3(i, j, k) = sf_noise;
            Pmat3(i, j, k) = gain_dB;

            Dmat3(i, j, k) = norm([Upos, Hdrone] - [PosUE, 0]);
%             f2 = getcostf2DK([Upos, Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
%             Cmat(i, j) = - f2;
% 
%             [f3, lopt] = getcostfK([Upos, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
%             Cmat3d(i, j) = - f3;

%             Hmat(i, j) = sqrt(lopt^2 - norm(Upos(1:2) - PosUE(1:2), 2)^2);
        end
        
    end
end

Pmat3sf = zeros(Npt, Npt, Nh);
% Add correlated shadow fading
%  - Filter design
stepSizeZ = 3;
dx0 = [0: (1/(Npt - 1) * WinLen): dcor * 2]; dx1 = [-dx0(end:-1:1), dx0(2:end)];
dz0 = [0: stepSizeZ: dcor * 2]; dz1 = [-dz0(end:-1:1), dz0(2:end)];
dX = repmat(dx1(:), 1, length(dx1));
dY = repmat(dx1(:).', length(dx1), 1);
Hcor = zeros(length(dx1), length(dx1), length(dz1));
for i = 1:length(dz1)
    Hcor(:, :, i) = exp(- sqrt(dX.^2 + dY.^2 + dz1(i)^2) / dcor);
end
Hcor = Hcor / norm(Hcor(:));

SF3 = imfilter(SFnoise3, Hcor);
PmatSF3 = Pmat3 + SF3;

ShadowFading.SF3 = SF3;
ShadowFading.Xrange = Xrange;
ShadowFading.Yrange = Yrange;
ShadowFading.Zrange = Zrange;

% --
[~, zid] = min(abs(Zrange - U.Hdrone));
Pmat = Pmat3(:, :, zid) + SF3(:, :, zid);

for i = 1:Npt
    for j = 1:Npt
        Upos = [Xrange(i), Yrange(j)]; Hdrone = Zrange(zid); los = LosMat3(i, j, zid);
        
        f2 = getcostf2DK_sf([Upos, Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun, SF3(i, j, zid));
        Cmat(i, j) = - f2;
        
        % [f3, lopt] = getcostfK_sf([Upos, Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun, ShadowFading);
        % Cmat3d(i, j) = - f3; 
        % However, the above is not the correct way to compute the optimal
        % capacity under shadowing (exhaustive search baseline)
        
        C = 0;
        for k = 1:Nh
            gainu = 10^(PmatSF3(i, j, k)/10) / U.Noise;
            
            du2b = norm([Upos, Zrange(k)] - [PosBS, U.Hbs], 2);
            gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
            f3 = - fun(gainu, gainb);
            if f3 > C
                C = f3;
            end
        end
        Cmat3d(i, j) = C;
    end
end

show_map(Xrange, Yrange, Pmat, [], 2); title('Power map');

show_map(Xrange, Yrange, Cmat3d, [], 3); clim = caxis;
title(sprintf('Capacity map (3d), maxCap = %3.2f', max(Cmat3d(:))));

show_map(Xrange, Yrange, Cmat, [], 4); 
title(sprintf('Capacity map (2d), maxCap = %3.2f', max(Cmat(:))));caxis(clim);

%%
figure(103), 
I1 = find(abs(LosMat3 - 1) < 1e-9); 
    hsf1 = cdfplot(SF3(I1)); hold on, set(hsf1, 'linestyle', '-');
I2 = find(abs(LosMat3 - (1-1/(K-1))) < 1e-9); 
    hsf2 = cdfplot(SF3(I2)); set(hsf2, 'linestyle', '-.');
I3 = find(abs(LosMat3 - (1-2/(K-1))) < 1e-9); 
    hsf3 = cdfplot(SF3(I3)); set(hsf3, 'linestyle', '--');
xi = randn(1, 10000) * 6; hxi = cdfplot(xi); set(hxi, 'color', 'red'); hold off
legend('empirical LOS SF', 'empirical OLOS SF', 'empirical NLOS SF', '6 dB Gaussian', 'location', 'southeast');

figure(104), 
acor_sf = xcorr(SF3(1, :, zid)); 
M_acor_sf = (length(acor_sf) - 1) / 2; acor_x = [- M_acor_sf:0, 1:M_acor_sf ] * WinLen / (Npt - 1);
plot(acor_x, acor_sf); xlim([-dcor, dcor] * 10);

figure(105), plot(Dmat3(:), PmatSF3(:), 'bx');
%% Search Algorithm for Optimal UAV Position in 3D
% 
PosUE3 = [PosUE, 0];
PosBS3 = [PosBS, U.Hbs];

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
maxStep = ceil((2.4 * U.K - 1.4) * norm(PosBS(1:2) - PosUE(1:2), 2) / stepSizeMeter);

Fmin = inf;
Xhat3 = [0, 0, 0];
Xhat2 = [0, 0];

XposArray = zeros(maxStep, 2);
FArray = zeros(maxStep, 1);
FminArray = zeros(maxStep, 1);
Fcritical = 0;

Xpos = PosBS;
Stage = 1;
cnt = 0;
Rhos = zeros(1, maxStep); % Record for Stage 1, rho value along the way
ks = zeros(1, maxStep);   % Record for stage 1, segment label along the way
ksearch = 1;
while Stage < 4 && cnt < maxStep
    cnt = cnt + 1;
    
    Xpos0 = Xpos;   % UAV Position
    
    % los = IsLosK(PosUE, Xpos0, BldLines, BldHeight, U.Hdrone, BldTypes);
    los = IsLosK_discrete([PosUE, 0], [Xpos0, U.Hdrone], Maps, stepsize, map_x0);
    ksegment = round((1 - los) * (U.K - 1) + 1);   % propagation segment index, k = 1,2,...,K
    [f, ~, Xpos3] = getcostfK([Xpos0, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
    FArray(cnt) = f;
    XposArray(cnt, :) = Xpos0;
    if f < Fmin
        Fmin = f;
        Xhat3 = Xpos3;
        Xhat2 = Xpos0;
    end
    
    if cnt == 1
        FminArray(cnt) = Fmin;
    else
        if f < FminArray(cnt - 1)
            FminArray(cnt) = Fmin;
        else
            FminArray(cnt) = FminArray(cnt - 1);
        end
    end
    
    if Stage == 1
        % Stage 1: Search on the User-BS axis
        if ksegment == 1 % LOS segment
            Rhos(cnt) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt) = 1;
            Rhos0 = findCriticalPointsK(Rhos, ks, [Xpos0, U.Hdrone], PosUE3, PosBS3, U, fun);
            
            rho0 = Rhos0(ksearch);
            theta0 = stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos = PosUE + (rho0 * M * Uvec(:)).';
            
            Fcritical = Fmin;
            
            Stage = 2;
            
        elseif norm(PosUE - Xpos0) > stepSizeMeter 
            Rhos(cnt) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt) = ksegment;
            searchDirection = (PosUE - PosBS) / norm(PosUE - PosBS);
            Xpos = Xpos0 + searchDirection * stepSizeMeter;
            uavColor = [0.7, 0.7, 0.7];
            
        else
            % Indoor case
            Stage = 4;
        end
        
    elseif Stage == 2
        % Stage 2: Search on the right branch
        searchDirection = uavSearchDirection3DK([Xpos, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos = Xpos0 + searchDirection * stepSizeMeter;
        
        if ksegment <= ksearch % LOS or virtual LOS
            uavColor = [1, 0, 0];
        else
            uavColor = [0, 1, 0];
        end
        
        if norm(searchDirection) < 1e-10
            rho0 = Rhos0(ksearch);
            theta0 = - stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos = PosUE + (rho0 * M * Uvec(:)).';
            
            Stage = 3;
            FminArray(cnt) = Fcritical;
        end
        
    elseif Stage == 3
        % Stage 3: Search on the left branch
        searchDirection = uavSearchDirection3DK([Xpos, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos = Xpos0 + searchDirection * stepSizeMeter;
        
        if ksegment <= ksearch % LOS or virtual LOS
            uavColor = [1, 0, 0];
        else
            uavColor = [0, 1, 0];
        end
        
        if norm(searchDirection) < 1e-10
            
            ksearch = ksearch + 1;
            if ksearch < U.K && 0
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

end
Xopt = Xhat3;
%% This is for generating a figure in the journal paper:
show_map(Xrange, Yrange, Cmat3d, [], 30); title('Capacity map (3d)');clim = caxis; hold on
show_contour(Xrange, Yrange, Cmat3d, [], 30); title('Capacity map (3d)');caxis(clim); hold on
figure(30), hold on
% title(sprintf('3D Capacity Map, obj = %f, pos = (%d,%d,%d)', Fmin, round(Xopt)));
for t  = 1:cnt
    % plot3(XposArray(t, 1), XposArray(t, 2), abs(FArray(t) - 1e-3), 's', 'linewidth', 1, 'markersize', 4, 'color', [0,1,0]); 
    plot3(XposArray(t, 1), XposArray(t, 2), abs(FminArray(t) - 1e-3), 's', 'linewidth', 2, 'markersize', 2, 'color', [0.2,0.2,0.2]); 
    plot3(XposArray(t, 1), XposArray(t, 2), 0, '.', 'linewidth', 1, 'markersize', 1, 'color', [0.5,0.5,0.5]); 
end

plot3(Xhat2(1), Xhat2(2), abs(Fmin - 1e-3), 'xr', 'linewidth', 4, 'markersize', 16); hold off

% %
xlim([400, 600]);
ylim([450, 600]);
view([-33, 25]);

figure(31), 
show_contour(Xrange, Yrange, Cmat3d, [], 31); title('Capacity map (3d)');caxis(clim); hold on
title(sprintf('3D Capacity Map, obj = %f, pos = (%d,%d,%d)', Fmin, round(Xopt)));
for t  = 1:cnt
    % plot3(XposArray(t, 1), XposArray(t, 2), abs(FArray(t) - 1e-3), 's', 'linewidth', 1, 'markersize', 4, 'color', [0,1,0]); 
    plot3(XposArray(t, 1), XposArray(t, 2), abs(FminArray(t) - 1e-3), 's', 'linewidth', 2, 'markersize', 2, 'color', [0.2,0.2,0.2]); hold on
end
plot3(Xhat2(1), Xhat2(2), abs(Fmin - 1e-3), 'xr', 'linewidth', 4, 'markersize', 16); hold off
hold off
% %
xlim([400, 600]);
ylim([450, 600]);

%% Search Algorithm for Optimal UAV Position in 2D
% 
Fmin2d = inf;
Xhat2d = [0, 0];

XposArray2d = zeros(maxStep, 2);
FArray2d = zeros(maxStep, 1);

Xpos2d = PosBS;
Stage = 1;
cnt2d = 0;
Rhos2d = zeros(1, maxStep); % Record for Stage 1, rho value along the way
ks = zeros(1, maxStep);   % Record for stage 1, segment label along the way
ksearch = 1;
while Stage < 4 && cnt2d < maxStep
    cnt2d = cnt2d + 1;
    
    Xpos0 = Xpos2d;   % UAV Position
    
    % los = IsLosK(PosUE, Xpos0, BldLines, BldHeight, U.Hdrone, BldTypes);
    los = IsLosK_discrete([PosUE, 0], [Xpos0, U.Hdrone], Maps, stepsize, map_x0);
    ksegment = round((1 - los) * (U.K - 1) + 1);   % propagation segment index, k = 1,2,...,K
    f = getcostf2DK([Xpos0, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
    FArray2d(cnt2d) = f;
    XposArray2d(cnt2d, :) = Xpos0;
    if f < Fmin2d
        Fmin2d = f;
        Xhat2d = Xpos0;
    end
    
    if Stage == 1
        % Stage 1: Search on the User-BS axis
        if ksegment == 1 % LOS segment
            Rhos2d(cnt2d) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt2d) = 1;
            Rhos02d = findCriticalPoints2DK(Rhos2d, ks, [Xpos0, U.Hdrone], PosUE3, PosBS3, U, fun);
            
            rho0 = Rhos02d(ksearch);
            theta0 = stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
            
            Stage = 2;
            
        elseif norm(PosUE - Xpos0) > stepSizeMeter 
            Rhos2d(cnt2d) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt2d) = ksegment;
            searchDirection = (PosUE - PosBS) / norm(PosUE - PosBS);
            Xpos2d = Xpos0 + searchDirection * stepSizeMeter;
            uavColor = [0.7, 0.7, 0.7];
            
        else
            % Indoor case
            Stage = 4;
        end
        
    elseif Stage == 2
        % Stage 2: Search on the right branch
        searchDirection = uavSearchDirection2DK([Xpos2d, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos2d = Xpos0 + searchDirection * stepSizeMeter;
        
        if ksegment <= ksearch % LOS or virtual LOS
            uavColor = [1, 0, 0];
        else
            uavColor = [0, 1, 0];
        end
        
        if norm(searchDirection) < 1e-10
            rho0 = Rhos02d(ksearch);
            theta0 = - stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
            
            Stage = 3;
        end
        
    elseif Stage == 3
        % Stage 3: Search on the left branch
        searchDirection = uavSearchDirection2DK([Xpos2d, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos2d = Xpos0 + searchDirection * stepSizeMeter;
        
        if ksegment <= ksearch % LOS or virtual LOS
            uavColor = [1, 0, 0];
        else
            uavColor = [0, 1, 0];
        end
        
        if norm(searchDirection) < 1e-10
            
            ksearch = ksearch + 1;
            if ksearch < U.K
                rho0 = Rhos02d(ksearch);
                theta0 = stepSizeMeter / rho0;
                M = [cos(theta0), -sin(theta0)
                     sin(theta0), cos(theta0)];
                Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
                Stage = 2;
            else
                Stage = 4; % Algorithm terminates
            end
        end
        
    else
        % The entire algorithm terminates
    end
    

end
Xopt = Xhat2d;
%
figure(4), hold on
title(sprintf('2D Capacity Map, obj = %f, pos = (%d, %d)', Fmin2d, round(Xopt)));
for t  = 1:cnt2d
    plot3(XposArray2d(t, 1), XposArray2d(t, 2), abs(FArray2d(t) - 1e-3), 's', 'linewidth', 1, 'markersize', 4, 'color', [0,1,0]); 
end
plot3(Xhat2d(1), Xhat2d(2), abs(Fmin2d - 1e-3), 'xk', 'linewidth', 2, 'markersize', 11); hold off
% %