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
% ------------------
% Massive simulation
% Version: Build and use a discrete urban city map.
% Date: July 5th, 2018.

% clear
addpath ../uavLibrary
addpath(genpath('../../../Library')),

simCaseID = 0;
simFileNameToken = 'uavpos3d_mmwrate';

Nue = 10000; % 10000 = 5.8 hours (all schemes)
PosBS = [150, 770];
SCENARIO = 3;

% Load city map
DATA = load('../uavDataset/urbanMap_DC.mat');
Maps = DATA.Maps;
[MapWidth, MapHeight] = size(Maps{1});
% Blds = DATA.BldArea;
meterPerPixel = DATA.stepsize;
map_x0 = DATA.x0;
Xvec = DATA.Xvec;
Yvec = DATA.Yvec;
K = length(Maps) + 1; U.K = K;
urbanMap = zeros(MapWidth, MapHeight);
for k = 1:K - 1
	urbanMap = urbanMap + Maps{k}; % urbanMap = FoliageMap + BldMap
end

% Channel configuration
% U.A0 = -22; U.B0 = -40; % UAV-BS channel
% U.A1 = -22.7; U.B1 = -28; U.S1 = 1; % Propagation segment 1 (LOS)
% U.A2 = -28.4; U.B2 = -24; U.S2 = 3; % Propagation segment 2
% U.A3 = -36.4; U.B3 = -22; U.S3 = 3; % Propagation segment 3 (NLOS)
% U.Alpha = [U.A1, U.A2, U.A3];
% U.Beta= [U.B1, U.B2, U.B3];
% U.Noise = 1e-10; U.Pb = 2; U.Pd = 2;

U.Huser = 0;    % meter, ground level
U.Hbs = 45;     % meter, BS height
U.Hmin = 45;    % meter, minimum UAV operation height
U.Hmax = 120;
U.Hdrone = 50;  % meter, UAV search height
U.minDistance = 50;    % Minimum distance for a LOS user to the BS
stepSizeMeter = 3;  % UAV search step size

%% Specify the objective function(s)
% NOTE THAT we redefine the objective function according to the 3D UAV
% positioning paper, where x and y in f(x,y) represent the distances
% (rather than gains in the 2D UAV paper).
mymin = @(x,y) (x + y - abs(x - y)) / 2;
% -- 
% General configurations
funs = cell(K, 1);
funs0 = cell(K, 1);
% 3GPP UMi model:
fc_3gpp = 2.5;    % carrier frequency = 2.5 GHz
floss_3gpp = 10;
bstxpower_dBm = 23;
uavtxpower_dBm = 20;
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
    % r0_3gpp{k} = @(d) min(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    r0_3gpp{k} = @(d) mymin(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    snr0_3gpp{k} = @(d) 10.^(snr0(d) / 10);
    
    snr1 = @(d) 174 + uavtxpower_dBm - pl3gpp{k}(d) - 10*log10(bw_3gpp) - noise_figure;
    % r1_3gpp{k} = @(d) min(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    r1_3gpp{k} = @(d) mymin(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    snr1_3gpp{k} = @(d) 10.^(snr1(d) / 10);
end

% mmWave model at 28 GHz
floss_mmw = 20; 
GainBf = 20; % 20 dB beamforming gain
bw_mmw = 500e6; % 500 MHz
plmmw = cell(3, 1); r0_mmw = cell(3, 1); r1_mmw = cell(3, 1);
plmmw{1} = @(d) 61.4 + 20 * log10(d);
plmmw{2} = @(d) 61.4 + 20 * log10(d) + floss_mmw;
plmmw{3} = @(d) 72.0 + 29.2 * log10(d);
for k = 1:3
    snr0 = @(d) 174 + bstxpower_dBm - plmmw{k}(d) - 10*log10(bw_mmw) - noise_figure + GainBf;
    % r0_mmw{k} = @(d) min(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    r0_mmw{k} = @(d) mymin(log2(1 + 10.^((snr0(d) - snr_backoff)/10)), modulation_cutoff);
    snr0_mmw{k} = @(d) 10.^(snr0(d) / 10);
    
    snr1 = @(d) 174 + uavtxpower_dBm - plmmw{k}(d) - 10*log10(bw_mmw) - noise_figure + GainBf;
    % r1_mmw{k} = @(d) min(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    r1_mmw{k} = @(d) mymin(log2(1 + 10.^((snr1(d) - snr_backoff)/10)), modulation_cutoff);
    snr1_mmw{k} = @(d) 10.^(snr1(d) / 10);
end

% % -------------------------------------------------------------------------
% Scenario Ia: Rate maximization problem I (K = 3 case). Similar to the 2D
% paper 
if SCENARIO == 11
    for k = 1:K
        gainu = @(x) 10 ^ ((U.Alpha(k) * log10(x) + U.Beta(k)) / 10) / U.Noise;
        gainb = @(y) 10 ^ ((U.A0 * log10(y) + U.B0) / 10) / U.Noise;
        % f = fun(gainu, gainb);
        fun = @(x,y) max(-log2(1 + U.Pd * real(gainu(x))), ...
                         -log2(1 + U.Pb * real(gainb(y)))) * 1e-6;
        funs{k} = fun;

        gainb2 = @(x) 10 ^ ((U.Alpha(k) * log10(x) + U.Beta(k)) / 10) / U.Noise;
        fun0 = @(x) -log2(1 + U.Pb * gainb2(x)) * 1e-6;
        funs0{k} = fun0;
    end
    obj_sign = - 1;
end

% -------------------------------------------------------------------------
% Scenario Ib: Rate maximization problem I (K = 3 case). Similar to the 2D
% paper 
if SCENARIO == 12
    for k = 1:3
        funs{k} = @(x,y) - 0.5 * min(r1_3gpp{k}(x), r0_3gpp{1}(y));
        funs0{k} = @(x) - r0_3gpp{k}(x);
    end
    obj_sign = - 1;
    
    sim_scheme = [1 1 1 0 1 0];
end

% -------------------------------------------------------------------------
% % % Scenario II: Outage probability minimization. Similar to the 2D paper
if SCENARIO == 2
    % target_rate = 2;

    for k = 1:K
        funs{k} = @(x,y) (1/snr1_3gpp{k}(x) + 1/snr0_3gpp{1}(y));% * (2^(2 * target_rate) - 1)^2;
        funs0{k} = @(x) (1/snr0_3gpp{k}(x));% * (2^(target_rate) - 1)^2;
    end

%     funs{1} = @(x,y) (1/snr1_mmw{1}(x) + 1/snr0_mmw{1}(y));
%     funs{2} = @(x,y) (1/snr1_mmw{2}(x) + 1/snr0_mmw{1}(y));
%     funs{3} = @(x,y) (1/snr1_3gpp{3}(x) + 1/snr0_mmw{1}(y));
%     funs0{1} = @(x) (1/snr0_mmw{1}(x));
%     funs0{2} = @(x) (1/snr0_mmw{2}(x));
%     funs0{3} = @(x) (1/snr0_3gpp{3}(x));

    obj_sign = 1;
    
    sim_scheme = [1 1 1 0 1 1];
end

% % -------------------------------------------------------------------------
% Scenario III: mmWave & 3GPP UMi hybrid system
if SCENARIO == 3
    % funs{1} = @(x,y) - 0.5 * min(r1_mmw{1}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    % funs{2} = @(x,y) - 0.5 * min(r1_mmw{2}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    % funs{3} = @(x,y) - 0.5 * min(r1_3gpp{3}(x) * bw_3gpp, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    funs{1} = @(x,y) - 0.5 * mymin(r1_mmw{1}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    funs{2} = @(x,y) - 0.5 * mymin(r1_mmw{2}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    funs{3} = @(x,y) - 0.5 * mymin(r1_3gpp{3}(x) * bw_3gpp, r0_mmw{1}(y) * bw_mmw) / 1e9;   % Tbit
    
    funs0{1} = @(x) - r0_mmw{1}(x) * bw_mmw / 1e9;   % Tbit
    funs0{2} = @(x) - r0_mmw{2}(x) * bw_mmw / 1e9;   % Tbit
    funs0{3} = @(x) - r0_3gpp{3}(x) * bw_3gpp / 1e9;   % Tbit

    obj_sign = - 1;
    
    sim_scheme = [1 1 0 0 1 0];
end

% % -------------------------------------------------------------------------
% % Scenario IV: Energy constrained deliveray
if SCENARIO == 4
    Battery_Wh = 80;     % Watt.hour = 3.6kJ, battery capacity
    Pu_hover = 200;      % Watt, hovering power
    Pu_cruise = 200;     % Watt, cruise power
    Pt_circuit = 2;      % watt, Circuit power
    cruise_velocity = 5; % m/s

    txtime = @(d) (Battery_Wh * 3600 - Pu_cruise * d * 2 / cruise_velocity) / (Pu_hover + Pt_circuit + 10^(uavtxpower_dBm / 10)/1000);
    funs{1} = @(x,y) - 0.5 * mymin(r1_mmw{1}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
    funs{2} = @(x,y) - 0.5 * mymin(r1_mmw{2}(x) * bw_mmw, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
    funs{3} = @(x,y) - 0.5 * mymin(r1_3gpp{3}(x) * bw_3gpp, r0_mmw{1}(y) * bw_mmw) / 1e12 * txtime(y);   % Tbit
    funs0{1} = @(x) - r0_mmw{1}(x) * bw_mmw / 1e12 * txtime(0);   % Tbit
    funs0{2} = @(x) - r0_mmw{2}(x) * bw_mmw / 1e12 * txtime(0);   % Tbit
    funs0{3} = @(x) - r0_3gpp{3}(x) * bw_3gpp / 1e12 * txtime(0);   % Tbit

    obj_sign = - 1;
    
    sim_scheme = [1 1 1 0 1 0];
end

% % -------------------------------------------------------------------------
% % Scenario V: Energy minimization for message delivery
if SCENARIO == 5
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
    
    sim_scheme = [0 1 1 0 1 1];
end

% % -------------------------------------------------------------------------
% Pre-calculation for speed-up
syms d1 d2
for i =1:3
    F = matlabFunction(feval(funs{1}, d1, d2));
    funs{1} = F;
    
    F = matlabFunction(feval(funs0{1}, d1));
    funs0{1} = F;
end

%% Generate user topology and air-to-ground propagation statistics
% Nuser = 5000;
% [PosArray, Angles, LosFreq] = sim_capacityK3_discrete_genstat(Nuser, Maps, map_x0, meterPerPixel);
% show_map(Xvec, Yvec, urbanMap, Blds, 1); hold on
% for i = 1:Nuser
%     plot3(PosArray(i, 1),PosArray(i, 2), 50, 'ro', 'markersize', 4);
% end
% hold off
% figure(2),
% plot(Angles / pi * 180, LosFreq);
% legend('LOS', 'OLOS', 'NLOS');
% title('Marginal distribution under given elevation angles');
% save topology_DC PosArray Angles LosFreq

% Nota that: LosFreq is the marginal distribution under a given elevation
% angle
DATA = load('../uavDataSet/topology_DC.mat');
PosArray = DATA.PosArray; Angles = DATA.Angles; LosFreq = DATA.LosFreq;
LosStat = struct(); LosStat.Angles = Angles; LosStat.LosFreq = LosFreq;
%%
N_scheme = 6;
Alg_scheme_name = {
    'Direct BS-User link'   % 1
    'Statistical Method'
    'Simple Search'
    '2D Optimal'
    '3D Optimal (proposed)'
    '3D Optimal (exhaustive)'
};

% sim_scheme = [0 1 1 0 1 0];

tic

Nue = min(size(PosArray, 1), Nue);
Rates0 = zeros(Nue, N_scheme);
strongUserIds = zeros(Nue, 1);
failIds = ones(Nue, 1);
parfor i = 1:Nue

    PosUE = PosArray(i, :);
    
    % los = IsLosK(PosUE, [PosBS, U.Hbs], BldLines, BldHeight, U.Hdrone, BldTypes);
    los = IsLosK_discrete([PosUE, 0], [PosBS, U.Hbs], Maps, meterPerPixel, map_x0);
    if los == 1
        strongUserIds(i) = 1;
        % continue    % We are only interested in the case where the direct BS-user link is blocked
    end
    if norm(PosUE - PosBS, 2) < U.minDistance % && los == 1
        % failIds(i) = 1;
        continue      % We ignore the strong user (close to BS with LOS condition)
    end
    
    try
        if sim_scheme(2)
            [FminStat, XhatStat] = ...
                finduavposStatK_mat3(PosUE, PosBS, U, funs, stepSizeMeter, Maps, meterPerPixel, map_x0, LosStat);
        else
            FminStat = 0;
        end
            
        if sim_scheme(3)
            [Fmin1, Xhat1] = ...
                finduavpos1d_mat3(PosUE, PosBS, U, funs, stepSizeMeter, Maps, meterPerPixel, map_x0);
        else
            Fmin1 = 0;
        end
        
        if sim_scheme(4)
            [Fmin2, Xhat2] = ...
                finduavpos_mat3(PosUE, PosBS, U, funs, stepSizeMeter, Maps, meterPerPixel, map_x0);
        else
            Fmin2 = 0;
        end
        
        if sim_scheme(5)
            [Fmin3, Xhat3] = ...
                finduavpos3d_mat3(PosUE, PosBS, U, funs, stepSizeMeter, Maps, meterPerPixel, map_x0);
        else 
            Fmin3 = 0;
        end
        
        if sim_scheme(6)
            [Fmin_exhst, Xhat_exhst] = ...
                finduavpos3d_matexhst3(PosUE, PosBS, U, funs, stepSizeMeter, Maps, meterPerPixel, map_x0);
        else
            Fmin_exhst = 0;
        end
        

        if sim_scheme(1)
            % Direct BS-user link
            k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
            d = norm([PosBS, U.Hbs] - [PosUE, 0], 2);
            % snr = 10 ^ ((U.Alpha(k) * log10(d) + U.Beta(k)) / 10) / U.Noise;
            F0 = funs0{k}(d);
        else
            F0 = 0;
        end
        
        failIds(i) = 0;
        
    catch
        F0 = 0; Fmin1 = 0; Fmin2 = 0; Fmin3 = 0; FminStat = 0; Fmin_exhst = 0;
        failIds(i) = 1;
        fprintf('error found at user_id = %d\n', i);
    end

    Rates0(i, :) = obj_sign * [F0, FminStat, Fmin1, Fmin2, Fmin3, Fmin_exhst];
end

toc

%% Save results
clear DATA
[~, fileHead, ~, dataFileName] = save_result(simCaseID, simFileNameToken);


%% Plot results
my_line_styles = {'-',...   % proposed
                  '--', ... % exhaustive 
                  '-.', ... % statistical, BS-user
                  ':'...    % 1d
                  }.';
co = [   0    0.4470    0.7410 % blue -> proposed
    0.8500    0.3250    0.0980 % orange -> 1d baseline
    0.9290    0.6940    0.1250 % yellow -> statistics baseline
    0.4940    0.1840    0.5560 % purple -> exhaustive baseline
    0.4660    0.6740    0.1880 % green -> BS baseline
    0.3010    0.7450    0.9330 % Tiffany blue
    0.6350    0.0780    0.1840 % red 
    ];

disp_line_styles = my_line_styles([3 4 3 1 2]);
schemes_to_show = [1 2 3 4 5 6];
schemes_to_show = schemes_to_show(schemes_to_show .* sim_scheme > 0);

validUserId = failIds < 1;
Rates = Rates0(validUserId, :);

Nue = size(Rates, 1);
n_schemes_to_show = length(schemes_to_show);

% X_data = zeros(100, n_schemes_to_show);
% F_data = zeros(100, n_schemes_to_show);

Npoints = 1000;
eX_data = zeros(Npoints, n_schemes_to_show);
eF_data = zeros(Npoints, n_schemes_to_show);

eX1_data = zeros(Npoints, n_schemes_to_show);
eF1_data = zeros(Npoints, n_schemes_to_show);

for i = 1:n_schemes_to_show
    n = schemes_to_show(i);
%     r_vec = Rates(:, n);
%     [F1,X1] = ksdensity(r_vec, 'Support', 'positive');
%     F1_cumsum = cumsum(F1);
%     
%     X_data(:, i) = X1(:);
%     F_data(:, i) = F1_cumsum(:) / max(F1_cumsum);
   
    [eF, eX] = ecdf(Rates(:, n));   % CDF of the original data
    [eF1, eX1] = ecdf(Rates(:, n) ./ Rates(:, 1));  % CDF of the gain w.r.t. baseline
    
    Nef = length(eF);
    I = 1 : (Nef - 1) / Npoints : Nef - (Nef - 1) / (2 * Npoints);
    I = round(I);
    
    eX_data(:, i) = eX(I(:));
    eF_data(:, i) = eF(I(:));

    Nef = length(eF1);
    I = 1 : (Nef - 1) / Npoints : Nef - (Nef - 1) / (2 * Npoints);
    I = round(I);
    eX1_data(:, i) = eX1(I(:));
    eF1_data(:, i) = eF1(I(:));
    
end

% CDF figure ---------------------------------
hf = figure(2);
switch(SCENARIO)
    case {11, 12, 3}
        set(groot,'defaultAxesColorOrder',co([5 3 2 1 4],:))
        % set(groot,'defaultAxesLineStyleOrder','-.|:|-.|-|--')
        p_handle = plot(eX_data, eF_data,'linewidth',2);
        set(p_handle, {'LineStyle'}, disp_line_styles(1:size(eF_data,2)));
        set(gca, 'FontSize', 14);
        legend(Alg_scheme_name{schemes_to_show}, 'location', 'southeast');
        set(gca, 'YTick', 0:0.2:1);
        ylabel('CDF');
        
        if SCENARIO == 3
            xlabel('Tera bit');
            xlim([0 2]);
        else
            xlabel('bps/Hz');
        end
        
    case 5
        set(groot,'defaultAxesColorOrder',co([1 2 3],:))
        p_handle = plot(eX_data(:, end:-1:1) / 1000, eF_data(:, end:-1:1),'linewidth',2);
        set(p_handle, {'LineStyle'}, my_line_styles([1 4 3]));
        set(gca, 'FontSize', 14);
        legend(Alg_scheme_name{schemes_to_show(end:-1:1)}, 'location', 'southeast');
        set(gca, 'YTick', 0:0.2:1);
        xlabel('Energy consumption [kJ]');
        ylabel('CDF');
        
    case 2
        plot(eX_data(:, end:-1:2), eF_data(:, end:-1:2),'linewidth',2);
        xlabel('BER(UAV)/BER(BS-USER)');
        ylabel('CDF');
        set(gca, 'YTick', 0:0.2:1);
        set(gca, 'xscale', 'log');
        legend(Alg_scheme_name{schemes_to_show(end:-1:2)}, 'location', 'southeast');
        %xlim([-20 10]);
        xlim([1e-3, 1e-0]);
end
    

% tune_figure,
% %
% figureFileName = 'ecdf'; 
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, figureFileName, {'fig','jpg', 'eps'});
% end


% BAR figure ---------------------------------

switch(SCENARIO)
    case {3, 4}
        hf = figure(3);
        % rateNoUav = Rates(:, 1);
        rate_to_sort = Rates(:, 2);
        [~, sortedIndex] = sort(rate_to_sort, 'ascend');

        low20percentileIndex = sortedIndex(1:round(Nue * 0.2));
        high20percentileIndex = sortedIndex(round(Nue * 0.8): end);

        RateLow = mean(Rates(low20percentileIndex, :), 1);
        RateMean = mean(Rates, 1);
        RateHigh = mean(Rates(high20percentileIndex, :), 1);
        bar([RateLow(schemes_to_show)
                  RateMean(schemes_to_show)
                  RateHigh(schemes_to_show)]);
        set(gca, 'FontSize', 14);
        legend(Alg_scheme_name{schemes_to_show}, 'location', 'northwest');
        set(gca, 'XTickLabel', {'20th percentile', 'Mean', 'Top 20th percentile'});

        % set(gca, 'YTick', 0:1:4);
        ylabel('Throughput [Terabit]');
        % label the bars
        Xdata = [RateLow(schemes_to_show)
                 RateMean(schemes_to_show)
                 RateHigh(schemes_to_show)];
        bartext = [];
        for i = 1:size(Xdata, 1)
            for j = 1:size(Xdata, 2)
                bartext(i, j) = text(i + (j - 3.0) * 0.150, Xdata(i, j) + 0.01, ...
                    sprintf('%1.2f', Xdata(i, j)), 'fontsize', 11);
            end
        end
        
        % Use the handles TH to modify some properties
        set(bartext,'Horizontalalignment','center',...
        'verticalalignment','bottom') ;
%         tune_figure,
        %
%         figureFileName = 'bar'; 
%         if exist('fileHead', 'var') == 1
%             save_figure(hf, fileHead, figureFileName, {'fig','jpg', 'eps'});
%         end
end

