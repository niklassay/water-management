%% Reset the Script Environment
clc;
close all;
clear all;

%% Simulation Parameters 

% Switch to use a monte Carlo simulation to determine the steady state
useMonteCarloSimulation = false;

% The count of iterations to be used in the monte carlo simulation, if
% enabled
monteCarloMaxIter = 500;

% Switch to use real rainfall and evaporation data 
useRealData = false;

% Set the file name of the weather dataset
realDataSource = 'bangladesh_weather_formatted.csv';

% Set the name of the location in the weather dataset to be used to
% simulate the evaporation
locationName = 'Barisal';

% Reservoir surface to volume ratio
% This exemplary value is calculated for the Rurstausee with data from
% https://de.wikipedia.org/wiki/Rurtalsperre
surface_to_volume_ratio = 7.83 / 202.6;

% Parameters of the value function
a1 = 1;
a2 = 2;
b1 = -2;
b2 = -3;

% Total capacity of the reservoir
M = 7;

% Discount factor
beta = 0.9;

% Number of periods for the forward simulation, will be overwritten when
% using real data with the actually available amount of period data
T = 1000;

% Parameters for rainfall distribution, if simulated
mu = 0;
sigma = 1;


% Count of different water levels
dimWL = 1000;

% The simulation stops if either the iteration limit or
% the tolerance limit is reached
% Iteration limit
valueFunctionMaxIter = 10000;
% Tolerance limit
valueFunctionTolerance = 1e-6;

% The count of periods to use to calculate the mean steady state levels
steadyStateMeanPeriods = 200;
% The steady state calculation stops if tolerance imit is reached
steadyStateTolerance = 0.0001;

%% Initialization

if(useRealData)
    % Load weather dataset
    data = readtable(realDataSource);
end

% Initialize the monte carlo simulation
if (~useMonteCarloSimulation)
    monteCarloMaxIter = 1;
end
monteCarloSteadyState = zeros(monteCarloMaxIter,1);

% Utility function of the the farmers
utilFar = @(x) (a1 / (1 + b1)) * x .^ (1 + b1);
% Utility function of the recreational users
utilRec = @(s,x) (a2 / (1 + b2)) * (s - x) .^ (1 + b2);

% Create a grid of the water levels to later evaluate the value function 
% for different water levels and periods
waterLevel = linspace(0,M,dimWL);

%% Rainfall simulation
if (useRealData)
    % Find column indices of the necessary data
    indData = find(strcmp(table2cell(data(1,:)), locationName));
    % Extract the rainfall data and calculate the average per year
    avgRain_mm = str2double(table2array(data(5:end,indData(1,4)))) * 12;
    
    % Exclude missing values
    avgRain_mm = avgRain_mm(~isnan(avgRain_mm));

    % Norm the amount of rain in mm to the volume of the reservoir
    avgRain = M*surface_to_volume_ratio*avgRain_mm/1000;

    % Note that because of the weather conditions of the location chosen in 
    % bangladesh, the amount of rainfall is not high enough to provide a
    % significant amount of value through irrigatio for the farmers. For 
    % this demo however, we want the reservoir to be refilled in order to 
    % give the farmers at least some utility. Thus we will
    % artificially increase the amount of rainfall here. 
    avgRain = avgRain * 4;

    % Expected value of the rainfall
    E_rain = mean(avgRain);
else
    % Use gauss-hermite to calculate the expected value of the rainfall
    % distribution
    
    % Number of nodes
    n = 10;
    
    [x_i,w] = GaussHermite(n);
    x_trans = sqrt(2)*sigma.*x_i + mu;
    
    % Expected value of the rainfall
    E_rain = (1/sqrt(pi)) * w' * exp(x_trans);
end

% Discretize and round the expected value of the rain to fit to
% the grid
dimE_Rain = round(E_rain/(M/dimWL));


%% Evaporation Simulation

if (useRealData)
    % Extract the data necessary for the evaporation simulation
    altitude = str2double(table2array(data(2,indData(1,1))));
    latitude_n = str2double(table2array(data(3,indData(1,1))));
    avgMaxTemp = str2double(table2array(data(5:end,indData(1,1))));
    avgMinTemp = str2double(table2array(data(5:end,indData(1,2))));
    avgMeanTemp = (avgMaxTemp + avgMinTemp) / 2;
    avgHumidity = str2double(table2array(data(5:end,indData(1,3))));

    % Calculate the average dew temperature using (lawrence2005relationship)
    % Inputs: avgHumidity in degrees celsius
    %         avgHumidity isurface_to_volume_ration percent
    % Output: avgDewTemp in degrees celsius
    avgDewTemp = (100 - avgHumidity) / 5;

    % Calculate the altitude in dergees south
    latitude_s = -latitude_n;

    % Calculate the average evaporation using (linacre1977simple)
    % Inputs: avgMeanTemp in degrees celsius
    %         avgDewTemp in degrees celsius
    %         latitude in degrees south    %CORRECT?
    %         altitude in metres
    % Output: avgEvap in mm per year
    avgEvap_mm = (700 .* (avgMeanTemp + 0.006 .* altitude)./(100 - latitude_s) ...
                + 15 .* (avgMeanTemp - avgDewTemp))./ (80 - avgMeanTemp) * 365.25; 

    % Exclude missing values
    avgEvap_mm = avgEvap_mm(~isnan(avgEvap_mm));

    % Note that because of the weather in the chosen location in
    % bangladesh, the evaporation is higher than the 

    % Norm the amount of rain in mm to the volume of the reservoir
    avgEvap = M*surface_to_volume_ratio*avgEvap_mm/1000;

    % Expected value of the evaporation
    E_Evap = mean(avgEvap);

    % Discretize and round the expected value of the evaporation to fit to
    % the grid
    dimE_Evap = round(E_Evap/(M/dimWL));
else 
    % Use no evaporation
    avgEvap = zeros(T,1);
    dimE_Evap = 0;
end

[V, optIrrigation_ind, aux_farm, aux_rec] = ValueFunction(dimWL, valueFunctionMaxIter, waterLevel, dimE_Rain, ...
                                                dimE_Evap, utilFar, utilRec, beta, valueFunctionTolerance);



%% (Monte Carlo) Simulation for the Steady State

for monteInd=1:monteCarloMaxIter
    if (useRealData)
        % Use real rainfall data
        r = avgRain;
    else
        % Use lognormal distributed rainfall
        r = exp(mu + sigma.*randn(T,1));
    end
    
    % Overwrite T to match available data points if using real data
    if (useRealData)
        T = size(r,1);
    end

    % Index for the current water level
    waterInd = zeros(1, T+1);
    % Set the index to 1 (empty) in period 1
    waterInd(1) = 1;

    % Amount of water (as index) over time to be used for irrigation
    irrigationInd = zeros(1, T);
    
    % Steady state levels of water in the reservoir over time
    steadyStateLvls = zeros(2, T);
    % Initialize the steady state period with the amount of
    % periods of the forward iteration in case no steady state will be
    % found
    steadyStatePeriod = T;
    
    % Perform a forward iteration
    for i=1:T
        % Irrigate with the optimal amount of water concerning the current
        % water level
        irrigationInd(i) = optIrrigation_ind(waterInd(i));
        
        % Calculate the water level in the following period
        waterInd(i+1) = min(max(waterInd(i) - irrigationInd(i) + round(r(i)/(M/dimWL)) ...
                            - round(avgEvap(i)/(M/dimWL)),1),dimWL);
        
        % Calculate the steady state level of the water reservoir
        if (i > steadyStateMeanPeriods)
            steadyStateLvls(1,i) = mean(waterLevel(waterInd(i - steadyStateMeanPeriods:i)));
        else
            steadyStateLvls(1,i) = mean(waterLevel(waterInd(1:i)));
        end
        
        % Calculate the change in the steady state levels
        if (i > 1)
            steadyStateLvls(2,i) = steadyStateLvls(1,i)...
                                   - steadyStateLvls(1,i-1);
        end
    end
    
    % Try to find a steady state
    if(~useMonteCarloSimulation)
        steadyStateLvl = steadyStateLvls(1,i);
    else
        for i=2:T
            % A steady state has been found if the changes in the steady state
            % levels are below the steady state tolerance
            if (abs(steadyStateLvls(2,i)) < steadyStateTolerance)
                steadyStatePeriod = i;
                steadyStateLvl = steadyStateLvls(1,i);
                fprintf('Steady State found in period %s\n', num2str(i));
                break;
            end
        end
    end
    % Save the steady state of the current monte carlo iteration
    monteCarloSteadyState(monteInd) = steadyStateLvl;
    fprintf('Iteration %s ended with Steady State %s\n', ...
             num2str(monteInd), num2str(steadyStateLvl));
end


%% plots

if(useMonteCarloSimulation)
    edges = [0:0.1:7];
    % Plot monte carlo steady states
    figure(1)
    hold on
    histogram(monteCarloSteadyState, edges, 'Normalization', 'Probability')
    title('Probability of Water Level Steady States');
    ylabel('Probability');
    xlabel('Water Level in of the Reservoir');
    hold off
else
    % Plot the value function
    figure(1)
    hold on
    plot(waterLevel,V(:));
    title('Value Function');
    xlabel('Water Level');
    ylabel('Maximum Value V');
    axis([0 M -20 0]);
    hold off
    
    % Plot water level, amount of water used for irrigation, rainfall and
    % evaporation
    figure(2)
    hold on
    plot(waterLevel(waterInd));    
    plot(waterLevel(irrigationInd));
    plot(r);
    plot(avgEvap);
    xlim([1 T]);
    legend('Water Level of the Reservoir','Water used for Irrigation',...
    'Rain','Evaporation');
    title('Optimal Irrigation Policy');
    xlabel('Period');
    ylabel('Amount');
    hold off

    % Plot water level histogram
    figure(3)
    hold on
    histogram(waterLevel(waterInd(steadyStatePeriod:end)),30, ... 
        'Normalization','Probability');
    title('Distribution of Water Level');
    ylabel('Probability');
    xlabel('Water Level of the Reservoir');
    hold off

    % Plot optimal irrigation policy
    figure(4)
    hold on
    plot(linspace(0,M,dimWL), waterLevel(optIrrigation_ind))
    title('Optimal Irrigation Policy')
    xlabel('Water Level of the Reservoir')
    ylabel('Irrigation Amount')
    hold off

    % Plot water level, steady state mean, mean-delta steady state level
    figure(5)
    hold on
    plot(waterLevel(waterInd));
    plot(steadyStateLvls(1,:));
    plot(steadyStateLvls(2,:));
    plot([1 T],[steadyStateLvl steadyStateLvl],'--g');
    xlim([1 T]);
    legend('Water Level of the Reservoir',strcat('Mean last ', num2str(steadyStateMeanPeriods), ' Periods'), ...
        'Difference in Mean to previous', 'Steady State Level');
    title('Steady State of Water Level');
    xlabel('Period');
    ylabel('Amount');
    hold off
    
    % Note: figures 6 and 7 require manual panning and zooming in order to
    % see the details

    % Plot the auxiliary matrix of the farmers
    figure(6)
    hold on
    surf(aux_farm);
    colormap(jet);
    title('Auxiliary Matrix of the Farmers');
    xlabel('Amount of Water Used For Irrigation');
    ylabel('Amount of Water in the Reservoir');
    zlabel('Value');
    hold off
    
    % Plot the auxiliary matrix of the recreational user
    figure(7)
    hold on
    surf(aux_rec);
    colormap(jet);
    title('Auxiliary Matrix of the Recreational Users');
    xlabel('Amount of Water Used For Irrigation');
    ylabel('Amount of Water in the Reservoir');
    zlabel('Value');
    hold off

    % Plot simulation steady state
    figure(8)
    hold on
    plot(waterLevel(waterInd));
    plot(steadyStateLvls(1,:));
    plot(steadyStateLvls(2,:));
    plot([1 T],[steadyStateLvl steadyStateLvl],'--g');
    area = patch([steadyStatePeriod-steadyStateMeanPeriods steadyStatePeriod steadyStatePeriod steadyStatePeriod-steadyStateMeanPeriods],[0 0 7 7],'r');
    alpha(area,.2)
    xlim([1 T]);
    ylim([-1 11]);
    legend('Water Level of the Reservoir',append('Mean last ', ... 
        num2str(steadyStateMeanPeriods), ' Periods'), ...
        'Difference in Mean to previous', 'Steady State Level');
    title('Steady State of Water Level');
    xlabel('Period');
    ylabel('Amount');
    set(gca,'FontSize',12)
    set(gcf,'Units','Centimeters','position',[0,0,16,12]);
    hold off
    
    
    % Plot optimal irrigation policy
    figure(9)
    hold on
    colormap_jet = colormap(parula);
    for beta=0.95:-0.05:0.95
        % Computation of the Value Function
        [V, optIrrigation_ind, aux_farm, aux_rec] = ValueFunction(dimWL, valueFunctionMaxIter, waterLevel, dimE_Rain, dimE_Evap, utilFar, utilRec, beta, valueFunctionTolerance);
        plot(linspace(0,M,dimWL), waterLevel(optIrrigation_ind), 'color', colormap_jet(round(256*beta),:))
    end
    title('Optimal Irrigation Policy')
    xlabel('Water Level of the Reservoir')
    ylabel('Irrigation Amount')
    hold off
    
    % Plot the aggregated use for farmers and recreational users
    figure(10)
    hold on
    set(gca,'FontSize',12)
    set(gcf,'Units','Centimeters','position',[0,0,16,8]);
    %title('Aggreated Use for Farmers and Recreational Users');
    yyaxis left;
    plot(utilRec(waterLevel(waterInd(:, 1:T)), waterLevel(irrigationInd)) + utilFar(waterLevel(irrigationInd)),'color','#A2142F');
    ylabel('Use');
    ylim([-100 0]);
    yyaxis right;
    plot(waterLevel(waterInd),'-', 'color', '#0072BD');
    plot(waterLevel(irrigationInd),'-','color','#D95319');
    plot(r,'-','color','#EDB120');
    ylabel('Amount');
    ylim([0 10]);
    legend('Aggregated Use','Water Level of the Reservoir','Water used for Irrigation',...
    'Rain');
    xlim([1 100]);
    xlabel('Period');
    hold off
    
    figure(11)
    set(gcf,'Units','Centimeters','position',[0,0,8,8]);
    edges = [0:0.05:7];
    hold on
    histogram(r, edges, 'Normalization', 'Probability')
    title('Probability of Rainfall');
    ylabel('Probability');
    xlabel('Amount of Rain');
    hold off
end



%% Computation of the Value Function
function [V, optIrrigation_ind, aux_farm, aux_rec] = ValueFunction(dimWL, valueFunctionMaxIter, waterLevel, dimE_Rain, ...
                                                dimE_Evap, utilFar, utilRec, beta, valueFunctionTolerance)

    % Create a value function for the combined value of farmers and
    % recreational users for each water level, respecting the discounted value
    % for the next period and the expected rainfall and evaporation
    V = zeros(dimWL,1);

    tic
    for j = 1:valueFunctionMaxIter
        V_old = V;

        % Create an auxiliary matrix for each water level and depending on 
        % that each possible amount of water to be used for irrigation. 
        % Then use the auxiliary matrix to find the maximum over all 
        % amounts of water used for irrigation for each water level in the 
        % reservoir
        aux = zeros(dimWL, dimWL) + NaN;
        aux_farm = zeros(dimWL, dimWL) + NaN;
        aux_rec = zeros(dimWL, dimWL) + NaN;
        for iWL = 1:dimWL
            for iIrrigation = 1:iWL
                aux(iWL, iIrrigation) = utilFar(waterLevel(iIrrigation)) + utilRec(waterLevel(iWL), waterLevel(iIrrigation)) ...
                    + beta*V_old(min(max(iWL-iIrrigation +1+dimE_Rain-dimE_Evap,1),dimWL));
                % Generate aux_farm and aux_rec for plotting in a later 
                % step 
                aux_farm(iWL, iIrrigation) = utilFar(waterLevel( iIrrigation));
                aux_rec(iWL, iIrrigation) = utilRec(waterLevel(iWL), waterLevel(iIrrigation));
            end
        end
        % Compute V and the optimum irrigation index (for each water level
        % step) using the bellman equation

        [V, optIrrigation_ind]= max(aux,[],2);

        % Termination check: break if the norm is smaller than the 
        % tolerance (excluding -inf values)
        if norm(V_old(V ~= -inf) - V(V ~= -inf)) < valueFunctionTolerance
            break;
        end
    toc
    end
end