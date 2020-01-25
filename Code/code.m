clc 
close all
clear all

%% turn monteCarloSimulation for distribution of the steady state on/off
monteCarloSimulation = true;
if(monteCarloSimulation)
    monteCarloMaxIter = 500;
else
    monteCarloMaxIter = 1;
end
monteCarloSteadyState = zeros(monteCarloMaxIter,1);

%% Parameters & utility functions

% set realdaten to true or false
realdaten = false;

% Define which real data you want to use. Has to correspond exactly to
% names from dataset
useRealDataFrom = 'Barisal';

% set use of real Evaporation data to true or false
if (realdaten)
    real_Evap = true;
else
    real_Evap = false;
end

% Preferences
a1 = 1;
a2 = 2;
b1 = -2;
b2 = -3;

M = 7; % Capacity Reservior
beta = 0.9; % Discount factor

% utility functions
utilFar = @(x) (a1 / (1 + b1)) * x .^ (1 + b1);
utilRec = @(s,x) (a2 / (1 + b2)) * (s - x) .^ (1 + b2);

%% Create grid and iteration parameters

% We need a grid to evaluate the value function for different water levels
% and different periods.
% We define dimWL different water levels
dimWL = 1000;
waterLevel = linspace(0,M,dimWL);

% Iteration parameters
maxIter = 10000;
tol = 1e-6;

% Number of periods for forward simulation (will be overwritten if
% using realdaten
T = 1000;

%% Rainfall simulation
% Rainfall parameters for simulation
mu = 0;
sigma = 1;

% Realdaten       
data = readtable('bangladesh_weather_formatted.csv');
% Find Column indices of data you want to use
indData = find(strcmp(table2cell(data(1,:)), useRealDataFrom));

% Extract the wanted rainfall data
avgRain = str2double(table2array(data(5:end,indData(1,4)))); % avg per Month

avgRain = avgRain * 12; % avg per Year

% Exclude NaNs
avgRain = avgRain(~isnan(avgRain));

% scale down to be sensible in our example
avgRain = avgRain / 900;

% calculate mu and sigma from the real dataset
% mu = mean(data);
% sigma = std(data);

%% Potentially read evaporation data
if (real_Evap)
    % Extract needed Data
    altitude = str2double(table2array(data(2,indData(1,1))));
    latitude = str2double(table2array(data(3,indData(1,1))));
    avgMaxTemp = str2double(table2array(data(5:end,indData(1,1))));
    avgMinTemp = str2double(table2array(data(5:end,indData(1,2))));
    avgMeanTemp = (avgMaxTemp + avgMinTemp) / 2;
    avgHumidity = str2double(table2array(data(5:end,indData(1,3))));

    % Calculate the Evaporation by formulas from paper
    e = (700 .* (avgMeanTemp + 0.006 .* altitude) ./...
    (100 - latitude) + 15 .* ((100 - avgHumidity) / 5))...
    ./ (80 - avgMeanTemp); % Avg Evaporation in mm per day

    e = e * 365.25; % Avg Evaporation in mm per year

    % Exclude NaNs
    e = e(~isnan(e));

    % scale down to be sensible in our example
    e = e / 900;

    % Expected value of evaporation
    E_Evap = mean(e);

    % Discretization of expected value of rain and rounding to fit our grid 
    dimE_Evap = round(E_Evap/(M/dimWL));
else 
    e = zeros(T,1);
    dimE_Evap = 0;
end
%% Gauss-Hermite to calculate the expected value of the rain distribution
if (realdaten)
    E_rain = mean(avgRain);
else
    n = 10;
    [x_i,w] = GaussHermite(n);

    x_trans = sqrt(2)*sigma.*x_i + mu;

    % Expected value of rain
    E_rain = (1/sqrt(pi)) * w' * exp(x_trans);
end

% Discretization of expected value of rain and rounding to fit our grid 
dimE_Rain = round(E_rain/(M/dimWL));

%% Compute Value function

% Create value function V for the combined value of farmers and 
% recreational users for each water level, respecting the discounted value
% for the next period and respecting the expected rainfall
V = zeros(dimWL,1);

tic
for j = 1:maxIter
    V_old = V;

    % We need an auxiliary matrix "aux" for each water level and depending
    % on that each possible amount of water used for irrigation. Then we'll
    % use aux to find the maximum over all amounts of water used for
    % irrigation for each water level in the reservoir
    aux = zeros(dimWL, dimWL) + NaN;
    aux_farm = zeros(dimWL, dimWL) + NaN;
    aux_rec = zeros(dimWL, dimWL) + NaN;
    for iWL = 1:dimWL
        for iIrrigation = 1:iWL
            aux(iWL, iIrrigation) = utilFar(waterLevel(iIrrigation)) +...
                utilRec(waterLevel(iWL), waterLevel(iIrrigation)) +...
                beta*V_old(min(max(iWL-iIrrigation+1+dimE_Rain-dimE_Evap,1),dimWL));
            aux_farm(iWL, iIrrigation) = utilFar(waterLevel(iIrrigation));
            aux_rec(iWL, iIrrigation) = utilRec(waterLevel(iWL), waterLevel(iIrrigation));
        end
    end
    [V, optIrrigation_ind]= max(aux,[],2);  

    % Termination check: Break if norm is smaller then tolerance for all
    % values that are not -inf
    if norm(V_old(V ~= -inf) - V(V ~= -inf)) < tol
        break;
    end
toc
end

% don't print all the diagrams if running the Monte Carlo Simulation
if(monteCarloSimulation==false)
    figure(6)
    hold on
    surf(aux_farm);
    colormap(jet);
    xlabel('Amount of Water Used For Irrigation');
    ylabel('Amount of Water in the Reservoir');
    zlabel('Value');
    hold off

    figure(7)
    hold on
    surf(aux_rec);
    colormap(jet);
    xlabel('Amount of Water Used For Irrigation');
    ylabel('Amount of Water in the Reservoir');
    zlabel('Value');
    hold off

    % Plot Value function
    figure(1)
    hold on
    plot(waterLevel,V(:));
    title('Value function');
    xlabel('Water Level');
    ylabel('maximum value V');
    axis([0 M -20 0]);
    hold off
end

%% Start of Monte Carlo Simulation for Steady Sta
for monteInd=1:monteCarloMaxIter
    %% Actually simulate the rainfall or use real data
    % rainfall: lognormal distribution with parameters mu and sigma
    if (realdaten)
        r = avgRain;
    else
        r = exp(mu + sigma.*randn(T,1));
    end
    
    %% forward iteration
    
    % Overwrite T to match available data points if using realdaten
    if (realdaten)
        T = size(r,1);
    end

    % Index for current water level; set to 1 (=empty) in period 1
    waterInd = zeros(1, T+1);
    waterInd(1) = 1;

    % Index for water that is used for irrigation
    irrigationInd = zeros(1, T);
    steadyStateLvls = zeros(2, T);
    steadyStatePeriod = T;
    periodsForMean = 200;
    for i=1:T

        irrigationInd(i) = optIrrigation_ind(waterInd(i));

        waterInd(i+1) = min(max(waterInd(i) - irrigationInd(i) + round(r(i)/(M/dimWL)) - round(e(i)/(M/dimWL)),1),dimWL);
        
        if (i>periodsForMean)
            steadyStateLvls(1,i) = mean(waterLevel(waterInd(i-periodsForMean:i)));
        else
            steadyStateLvls(1,i) = mean(waterLevel(waterInd(1:i)));
        end
        if (i~=1)
            % the differences to the previous mean
            steadyStateLvls(2,i) = steadyStateLvls(1,i)-steadyStateLvls(1,i-1);
        end
    end

    for i=2:T
        if (abs(steadyStateLvls(2,i)) < 0.0001)
            steadyStatePeriod = i;
            steadyStateLvl = steadyStateLvls(1,i);
            fprintf('Steady state found in period %s\n',num2str(i));
            break;
        end
    end
    %steadyStatePeriod = 20; % Change this parameter to determine the period
                               % after which the system is expected to be in a 
                               % steady state (has to be smaller then T!)

    %% TODO better threshold calculation by using more that one period to evaluate the difference
    % steadyStateLvl = mean(waterLevel(waterInd(steadyStatePeriod:end)));
    monteCarloSteadyState(monteInd) = steadyStateLvl;
    fprintf('Iteration %s ended with steady State: %s\n',num2str(monteInd), num2str(steadyStateLvl));
end

if(monteCarloSimulation)
    edges= [0:0.1:7];
    figure(1)
    hold on
    histogram(monteCarloSteadyState, edges, 'Normalization', 'probability')
    title('Probability of Water Level Steady States');
    ylabel('probability');
    xlabel('water level in reservoir');
    hold off
else
    figure(8)
    hold on
    plot(waterLevel(waterInd));
    plot(steadyStateLvls(1,:));
    plot(steadyStateLvls(2,:));
    plot([1 T],[steadyStateLvl steadyStateLvl],'--g');
    area = patch([steadyStatePeriod-periodsForMean steadyStatePeriod steadyStatePeriod steadyStatePeriod-periodsForMean],[0 0 7 7],'r');
    alpha(area,.2)
    xlim([1 T]);
    ylim([-1 11]);
    legend('water level in reservoir',append('Mean last ', num2str(periodsForMean), ' periods'), ...
        'Difference in mean to previous', 'steady state level');
    title('Steady State of Water Level');
    xlabel('period');
    ylabel('amount');
    set(gca,'FontSize',12)
    set(gcf,'Units','Centimeters','position',[0,0,16,12]);
    hold off
    
    % Plot water level, amount of water used for irrigation and steady state
    % value
    figure(2)
    hold on
    plot(waterLevel(waterInd));    
    plot(waterLevel(irrigationInd));
    plot(r);
    xlim([1 T]);
    legend('water level in reservoir','water used for irrigation','rain');
    title('Optimal Irrigation Policy');
    xlabel('period');
    ylabel('amount');
    hold off

    % plot histogram
    figure(3)
    hold on
    histogram(waterLevel(waterInd(steadyStatePeriod:end)),30, 'Normalization','probability');
    title('Distribution of Water Level');
    ylabel('probability');
    xlabel('water level in reservoir');
    hold off

    % plot optimal irrigation policy
    figure (4)
    hold on
    plot(linspace(0,M,dimWL), waterLevel(optIrrigation_ind))
    title('Optimal Irrigation Policy')
    xlabel('water level in reservoir')
    ylabel('irrigation amount')
    hold off

    % Plot water level, amount of water used for irrigation and steady state
    % value
    figure(5)
    hold on
    plot(waterLevel(waterInd));
    plot(steadyStateLvls(1,:));
    plot(steadyStateLvls(2,:));
    plot([1 T],[steadyStateLvl steadyStateLvl],'--g');
    xlim([1 T]);
    legend('water level in reservoir',strcat('Mean last ', periodsForMean, ' periods'), ...
        'Difference in mean to previous', 'steady state level');
    title('Steady State of Water Level');
    xlabel('period');
    ylabel('amount');
    hold off
end
