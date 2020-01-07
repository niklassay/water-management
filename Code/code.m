clc 
close all
clear all
%% Parameters & utility functions

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

% Number of periods for forward simulation
T = 100;

%% Rainfall simulation
% Rainfall parameters (might be overwritten by real data)
mu = 0;
sigma = 1;

% % Realdaten
% c = 4; % Column: 1,2,3,4 = POONDI,CHOLAVARAM,REDHILLS,CHEMBARAMBAKKAM
% data = csvread('chennai_reservoir_rainfall_formatted.csv',1,c,[1 c 16 c]); % Syntax csvread(filename,R1,C1,[R1 C1 R2 C2])
% mu = mean(data);
% sigma = std(data);

% rainfall: lognormal distribution with parameters mu and sigma
% r = lognrnd(mu, sigma, T, 1); % Alternative
r = exp(mu + sigma.*randn(T,1)); 

%% Gauss-Hermite to calculate the expected value of the rain distribution
n = 10;
[x_i,w] = GaussHermite(n);

x_trans = sqrt(2)*sigma.*x_i + mu;

% Expected value of rain
E_rain = (1/sqrt(pi)) * w' * exp(x_trans);
% E_rain = lognstat(mu, sigma) % Alternative

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
    for iWL = 1:dimWL
        for iIrrigation = 1:iWL
            aux(iWL, iIrrigation) = utilFar(waterLevel(iIrrigation)) +...
                utilRec(waterLevel(iWL), waterLevel(iIrrigation)) +...
                beta*V_old(min(max(iWL-iIrrigation+1+dimE_Rain,1),dimWL));
        end
    end
    [V, optIrrigation_ind]= max(aux,[],2);  
    
    % Termination check: Break if norm is smaller then tolerance for all values that are not -inf        
    if norm(V_old(V ~= -inf) - V(V ~= -inf)) < tol
        break;
    end
toc
end

% Plot Value function
figure(1)
hold on
plot(waterLevel,V(:));
title('Value function');
xlabel('Water Level');
ylabel('maximum value V');
axis([0 M -20 0]);
hold off

%% forward iteration

% Index for current water level; set to 1 (=empty) in period 1
waterInd = zeros(1, T+1);
waterInd(1) = 1;

% Index for water that is used for irrigation
irrigationInd = zeros(1, T);

for i=1:T
    
    irrigationInd(i) = optIrrigation_ind(waterInd(i));
    
    waterInd(i+1) = min(waterInd(i) - irrigationInd(i) + round(r(i)/(M/dimWL)),dimWL);
    
end

steadyStateThreshold = 20; % Change this parameter to determine the period 
                           % after which the system is expected to be in a 
                           % steady state (has to be smaller then T!)
                       
steadyStateLvl = mean(waterLevel(waterInd(steadyStateThreshold:end)));

% Plot water level, amount of water used for irrigation and steady state
% value
figure(2)
hold on
plot(waterLevel(waterInd));    
plot(waterLevel(irrigationInd));
plot([1 T],[steadyStateLvl steadyStateLvl],'--g');
xlim([1 T]);
legend('water left','water used for irrigation','steady state level');
title('Optimal Irrigation Policy');
xlabel('period');
ylabel('water level in reservoir');
hold off

% plot histogram
figure(3)
hold on
histogram(waterLevel(waterInd(steadyStateThreshold:end)),30, 'Normalization','probability');
title('Steady-State distribution of water level');
ylabel('probability');
xlabel('water level in reservoir');
hold off