%% parameters
% preference parameters
a1 = 1;
a2 = 2;
b1 = -2;
b2 = -3;

% maximum capacity of the reservoir
totalReservoirCapacity = 7;

% initial fill level of the reservoir
initialReservoirCapacity = totalReservoirCapacity/2;

% discount factor
beta = 0.9;

% rainfall distribution
mu = 0; 
sigma = 1;

% simulation periods
periods = 365;


%% utility functions
% farmers utility
uFarmer = @(x) a1/(1+b1)*x.^(1+b1);

% recreational users utility
uRecr = @(s,x) a2/(1+b2)*(s-x).^(1+b2);

% total utility
uTotal = @(s,x) uFarmer(x) + uRecr(s,x);


%% other functions
% rainfall
r = lognrnd(mu,sigma, 1, periods);

%% simulation

reservoirFillLevels = zeros(periods,1);
reservoirFillLevels(1) = initialReservoirCapacity;

for p=2:periods
    reservoirFillLevels(p) = reservoirFillLevels(p-1)+r(p);
    reservoirFillLevels(p) = min(totalReservoirCapacity, reservoirFillLevels(p));
end
reservoirFillLevels

