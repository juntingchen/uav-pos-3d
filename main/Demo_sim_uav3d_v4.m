%% Optimal UAV Placement in 3D Space
% Version 4: Change the objective function to depends on the distances
% (rather than the channel gains). Change the objective function
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
U.Hmax = 120;
% PosUE = [342, 237];
% PosUE = [94.0212  739.6973];
% PosUE = [498.8438  504.1995];
% PosUE = [99, 400];
% PosUE = [171, 258];
PosUE = [512, 644];
% PosUE = [  220.6620  743.5349];
PosBS = [150, 770];

% DATA = load('../uavDataset/urbanMapK_1m.mat');
DATA = load('../uavDataset/urbanMap_DC.mat');
Maps = DATA.Maps;
[MapWidth, MapHeight] = size(Maps{1});
Blds = DATA.BldArea;
stepsize = DATA.stepsize;
meterPerPixel = stepsize;
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
U.A3 = -36.4; U.B3 = -22; U.S3 = 3; % Propagation segment 3 (NLOS)
% U.K = 3;
U.Alpha = [-22.7, -28.4, -36.5];
U.Beta= [-28, -24, -22];
U.Noise = 1e-10; U.Pb = 2; U.Pd = 2;
U.Hbs = 45;
U.Hmin = 45; % Minimum UAV operation height
U.Hdrone = 50;
% Rmap.Alpha = [-20, -40]; Rmap.Beta = [-40, -40]; Rmap.Sigma2 = [3, 8]; Rmap.Pi = [0.3, 0.7];
% Chann.A0 = U.A0; Chann.B0 = U.B0; Chann.Noise = U.Noise;

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

% fun = @(x,y) max(-log2(1 + U.Pd * real(x)), -log2(1 + U.Pb * real(y)));
% fun = @(x,y) log(1 ./ (U.Pd * real(x)) + 1 ./ (U.Pb * real(y) * real(x) * U.Pd));

% General configurations
funs = cell(K, 1);
funs0 = cell(K, 1);
% 3GPP UMi model:
fc_3gpp = 2.5;    % carrier frequency = 2.5 GHz
floss_3gpp = 10;
bstxpower_dBm = 20;
uavtxpower_dBm = 17;
noise_figure = 7;   % 7 dB
bw_3gpp = 20e6;     % 20 MHz
snr_backoff = 3;    % 3 dB SNR backoff for rate estimation compared from AWGN channel
modulation_cutoff = 9.84;
pl3gpp = cell(3, 1); r1_3gpp = cell(3, 1); r0_3gpp = cell(3, 1); 
snr0_3gpp = cell(3, 1); snr1_3gpp = cell(3, 1);
snr0_mmw = cell(3, 1); snr1_mmw = cell(3, 1);
pl3gpp{1} = @(d) 22.0 + 28.0 * log10(d) + 20 * log10(fc_3gpp);
pl3gpp{2} = @(d) 22.0 + 28.0 * log10(d) + 20 * log10(fc_3gpp) + floss_3gpp;
pl3gpp{3} = @(d) 22.7 + 36.7 * log10(d) + 26 * log10(fc_3gpp);
for k = 1:3
    snr0 = @(d) 174 + bstxpower_dBm - pl3gpp{k}(d) - 10*log10(bw_3gpp) - noise_figure;
    r0_3gpp{k} = @(d) min(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    snr0_3gpp{k} = @(d) 10.^(snr0(d) / 10);
    
    snr1 = @(d) 174 + uavtxpower_dBm - pl3gpp{k}(d) - 10*log10(bw_3gpp) - noise_figure;
    r1_3gpp{k} = @(d) min(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    snr1_3gpp{k} = @(d) 10.^(snr1(d) / 10);
end

% mmWave model at 28 GHz
floss_mmw = 20; 
GainBf = 25; % 25 dB beamforming gain
bw_mmw = 1000e6; % 1 GHz
plmmw = cell(3, 1); r0_mmw = cell(3, 1); r1_mmw = cell(3, 1);
plmmw{1} = @(d) 61.4 + 20 * log10(d);
plmmw{2} = @(d) 61.4 + 20 * log10(d) + floss_mmw;
plmmw{3} = @(d) 72.0 + 29.2 * log10(d);
for k = 1:3
    snr0 = @(d) 174 + bstxpower_dBm - plmmw{k}(d) - 10*log10(bw_mmw) - noise_figure + GainBf;
    r0_mmw{k} = @(d) min(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    snr0_mmw{k} = @(d) 10.^(snr0(d) / 10);
    
    snr1 = @(d) 174 + uavtxpower_dBm - plmmw{k}(d) - 10*log10(bw_mmw) - noise_figure + GainBf;
    r1_mmw{k} = @(d) min(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    snr1_mmw{k} = @(d) 10.^(snr1(d) / 10);
end

% % -------------------------------------------------------------------------
% % Scenario Ia: Rate maximization problem I (K = 3 case). Similar to the 2D
% % paper 
% SCENARIO = 11;
% for k = 1:K
%     gainu = @(x) 10 ^ ((U.Alpha(k) * log10(x) + U.Beta(k)) / 10) / U.Noise;
%     gainb = @(y) 10 ^ ((U.A0 * log10(y) + U.B0) / 10) / U.Noise;
%     % f = fun(gainu, gainb);
%     fun = @(x,y) max(-log2(1 + U.Pd * real(gainu(x))), ...
%                      -log2(1 + U.Pb * real(gainb(y)))) * 1e-6;
%     funs{k} = fun;
%     
%     gainb2 = @(x) 10 ^ ((U.Alpha(k) * log10(x) + U.Beta(k)) / 10) / U.Noise;
%     fun0 = @(x) -log2(1 + U.Pb * gainb2(x)) * 1e-6;
%     funs0{k} = fun0;
% end
% obj_sign = - 1;

% -------------------------------------------------------------------------
% % Scenario Ib: Rate maximization problem I (K = 3 case). Similar to the 2D
% % paper 
% SCENARIO = 12;
% for k = 1:3
%     funs{k} = @(x,y) - min(r1_3gpp{k}(x), r0_3gpp{1}(y));
%     funs0{k} = @(x) - r0_3gpp{k}(x);
% end
% obj_sign = - 1;

% -------------------------------------------------------------------------
% % Scenario II: Outage probability minimization. Similar to the 2D paper
% SCENARIO = 2;
% target_rate = 2;
% for k = 1:K
%     funs{k} = @(x,y) (1/snr1_mmw{k}(x) + 1/snr0_mmw{1}(y));% * (2^(2 * target_rate) - 1)^2;
%     funs0{k} = @(x) (1/snr0_mmw{k}(x));% * (2^(target_rate) - 1)^2;
% end
% obj_sign = 1;

% % -------------------------------------------------------------------------
% Scenario III: mmWave & 3GPP UMi hybrid system
% SCENARIO = 3;
% funs{1} = @(x,y) - 0.5 * min(r1_mmw{1}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
% funs{2} = @(x,y) - 0.5 * min(r1_mmw{2}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
% funs{3} = @(x,y) - 0.5 * min(r1_3gpp{3}(x) * bw_3gpp, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
% funs0{1} = @(x) - r0_mmw{1}(x) * bw_mmw / 1e9;   % Tbit
% funs0{2} = @(x) - r0_mmw{2}(x) * bw_mmw / 1e9;   % Tbit
% funs0{3} = @(x) - r0_3gpp{3}(x) * bw_3gpp / 1e9;   % Tbit
% 
% obj_sign = - 1;

% % -------------------------------------------------------------------------
% % Scenario IV: Energy constrained deliveray
% SCENARIO = 4;
% Battery_Wh = 80;     % Watt.hour = 3.6kJ, battery capacity
% Pu_hover = 200;      % Watt, hovering power
% Pu_cruise = 200;     % Watt, cruise power
% Pt_circuit = 2;      % watt, Circuit power
% cruise_velocity = 5; % m/s
% 
% txtime = @(d) (Battery_Wh * 3600 - Pu_cruise * d * 2 / cruise_velocity) / (Pu_hover + Pt_circuit + 10^(uavtxpower_dBm / 10)/1000);
% funs{1} = @(x,y) - min(r1_mmw{1}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
% funs{2} = @(x,y) - min(r1_mmw{2}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
% funs{3} = @(x,y) - min(r1_3gpp{3}(x) * bw_3gpp, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
% funs0{1} = @(x) - r0_mmw{1}(x) * bw_mmw / 1e12 * txtime(0);   % Tbit
% funs0{2} = @(x) - r0_mmw{2}(x) * bw_mmw / 1e12 * txtime(0);   % Tbit
% funs0{3} = @(x) - r0_3gpp{3}(x) * bw_3gpp / 1e12 * txtime(0);   % Tbit
% 
% obj_sign = - 1;

% % -------------------------------------------------------------------------
% % Scenario V: Energy minimization for message delivery
SCENARIO = 5;
Pu_hover = 200;      % Watt, hovering power
Pu_cruise = 200;     % Watt, cruise power
Pt_circuit = 2;      % watt, Circuit power
cruise_velocity = 5; % m/s
Bits_Gb = 10;    % Giga bits

txtime = cell(3, 1);
for k = 1:3
    txtime{k} = @(d) Bits_Gb * 1e9 / (bw_mmw * 0.5 * r1_mmw{k}(d));
    funs{k} = @(d1, d0) Pu_cruise * d0 / cruise_velocity ...
        + (Pu_hover + Pt_circuit + 10^(uavtxpower_dBm / 10)/1000) * txtime{k}(d1);
end
obj_sign = 1;



%% Urban Map
uid = round(PosUE / stepsize);
show_map(Xvec, Yvec, urbanMap, Blds, 1);hold on
plot3(PosBS(1), PosBS(2), 50, 'r^', 'linewidth', 2, 'markersize', 9);
plot3(PosUE(1), PosUE(2), 50, 'ro', 'linewidth', 2, 'markersize', 9);hold off


%% Radio Map

Npt = 70;

% Xrange = Xvec(1) + (0:1/(Npt - 1):1) * MapWidth;
% Yrange = Yvec(1) + (0:1/(Npt - 1):1) * MapHeight;

% Lshow = 50;
% Xrange = PosUE(1) + ((0:1/(Npt - 1):1) - 0.5) * Lshow;
% Yrange = PosUE(2) + ((0:1/(Npt - 1):1) - 0.5) * Lshow;

midpos = (PosUE(1:2) + PosBS(1:2)) / 2;
L = norm(PosUE - PosBS);
Xrange = midpos(1) + ((0:1/(Npt - 1):1) - 1/2) * L;
Yrange = midpos(2) + ((0:1/(Npt - 1):1) - 1/2) * L;

[Xmat, Ymat] = meshgrid(Xrange, Yrange);
Pmat = zeros(Npt);Cmat = zeros(Npt);Cmat3d = zeros(Npt);Hmat = zeros(Npt);
for i = 1:Npt
    for j = 1:Npt
        Upos = [Xrange(i), Yrange(j)];
        du2u = norm([Upos, U.Hdrone] - [PosUE, 0]);
        % los = IsLosK(PosUE, Upos, BldLines, BldHeight, U.Hdrone, BldTypes);
        los = IsLosK_discrete([PosUE, 0], [Upos, U.Hdrone], Maps, stepsize, map_x0);
        if abs(los - 1) < 1e-9 % propagation segmen 1
            gain_dB = log10(du2u) * U.A1 + U.B1;
            
        elseif abs(los - (1-1/(K-1))) < 1e-9     % propagation segment 2
            gain_dB = log10(du2u) * U.A2 + U.B2;
        else
            gain_dB = log10(du2u) * U.A3 + U.B3;
        end
        Pmat(i, j) = gain_dB;
        
        % f2 = getcostf2d([Upos, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, funs);
        % Cmat(i, j) = f2;
        
        [f3, lopt] = getcostf3d([Upos, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, funs);
        Cmat3d(i, j) = f3;
        
        Hmat(i, j) = sqrt(lopt^2 - norm(Upos(1:2) - PosUE(1:2), 2)^2);
        
    end
end

show_map(Xrange, Yrange, Pmat, [], 2); title('Power map');
show_map(Xrange, Yrange, Cmat3d, [], 3); title('Capacity map (3d)');clim = caxis;
% show_map(Xrange, Yrange, Cmat, [], 4); title('Capacity map (2d)');caxis(clim);

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

        elseif norm(PosUE - Xpos0) > stepSizeMeter 
            Rhos(cnt) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt) = ksegment;
            searchDirection = (PosUE - PosBS) / norm(PosUE - PosBS);
            Xpos = Xpos0 + searchDirection * stepSizeMeter;

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

end
Xopt = Xhat3;
Fmin = obj_sign * Fmin;
%
figure(3), hold on
title(sprintf('3D Capacity Map, obj = %f, pos = (%d,%d,%d)', Fmin, round(Xopt)));
for t  = 1:cnt
    plot3(XposArray(t, 1), XposArray(t, 2), abs(FArray(t) - 1e-3), 's', 'linewidth', 1, 'markersize', 4, 'color', [0,1,0]); 
end
plot3(Xhat2(1), Xhat2(2), abs(Fmin - 1e-3), 'xk', 'linewidth', 2, 'markersize', 11); hold off
% %