%% Assignment 2 Portfolio Management 
% Long-term strategic asset allocation
clear
clc
% initialisation
Data = xlsread('DataAssignment2.xls');
Data(:,1) = [];
% Brandt parameters
N = 10000;
K = 18;
o=size(Data,1);
%% Exercise 1 ( Lecture 6 )

% a : VAR(1) model 
% Data parsing
r =  log(1 + Data(:,1));
s = Data(:,2);
y = [r, s];

x = [ones(o,1) s];
x = x(1:o-1,:);
y = y(2:end,:);
b = inv(x'*x)*x'*y;
residuals = x*b-y;
Sigma = cov(residuals);

% Make the restricted VAR(1) model
% b : Assume the new parameters of Brandt
a = [0.2049;-0.1694];
a_vol = [0.0839;0.0845];
AR_1 = [0.0568;0.9514];
AR_1_vol = [0.0249;0.0251];
Q = [0.006225 -0.006044 ; -0.006044 0.006316];
s0 = -3.5;


%% simulate N paths of K periods by the plug-in method

% loop over all paths
% simulate per period
for k = 1:K
    % loop over all paths
    for n = 1:N
    	eps = mvnrnd([0;0],Q);
       	if k == 1
           r_PI_sim(n,k) = a(1) + s0*AR_1(1) + eps(1);
           s_PI_sim(n,k) = a(2) + s0*AR_1(2) + eps(2);
       	else
           r_PI_sim(n,k) = a(1) + s_PI_sim(n,k-1)*AR_1(1) + eps(1);
           s_PI_sim(n,k) = a(2) + s_PI_sim(n,k-1)*AR_1(2) + eps(2);
       	end
    end
    Cumu_PI(k) = var(sum(r_PI_sim,2))/k;
end



%% c    NEW MODEL rt+1 
mu_r = 0.0072;
SSE = 0.3828;
obs = 60;

%% simulate N paths of K periods by the Decision-theory method

% verschil met chris: Hier wordt geen wortel genomen, not sure if must.
for n = 1:N
    sig(n) = sqrt(igammarnd(SSE,obs-1));
end
mu = normrnd(mu_r,sig/obs);

r_DT_sim = [];
for k = 1:K
   eps = normrnd(0,sig);
   r_DT_sim(:,k) = mu + eps;
   cumu_DT(k) = var(sum(r_DT_sim,2))/k;
end




%% Exercise 2 ( Lecture 7 & 8 )
% a calculate variances
%used to calculate 0.5 * var(r_t+1, r_t+2)
var1 = 6.225*10^-3;
var2 = 5.892*10^-3;
covar = -6.044*10^-3;
beta1r = 0.0568;

var_fut = beta1r^2 * var2 + var1;
temp_covar = beta1r*covar; 

var_scaled = 0.5*(var1+var_fut + 2*temp_covar);

% b discuss (no code needed), see pdf File (Report)

% c discuss (no code needed), see pdf File (Report)

%% Exercise 3 ( Lecture 9 )
% Numerically calculate optimal portfolio weights

%% a Myopic Strategy
% initialisation
gamma = 5;                  % CRRA coefficient 
path = r_PI_sim(:,1);       % simulated path from 1b
alphaG = linspace(0,1,100); % portfolio weights grid
rf = 0;
CU_Best = -inf;
alphaBest = -inf;

% Calculate optimal portfolio
for i = 1:length(alphaG)
    alpha = alphaG(i);
    RU = ((alpha*exp(path)+(1-alpha)*exp(rf)).^(1-gamma))/(1-gamma);
    CU = mean(RU);
    if (CU > CU_Best)
        CU_Best = CU;
        alphaBest = alpha;
    end
end

% Printing the results
disp(CU_Best)
disp(alphaBest)

%% b Dynamic Strategy
% initialisation
gamma = 5;                  % CRRA coefficient 
path = r_PI_sim;            % simulated path from 1b
rf = 0;
alphaG = linspace(0,1,100); % portfolio weights grid
U_fut = ones(N,1)*0.1; % initialise future utility
CU_Best = ones(N,1)*-inf; % intialise optimal utility
alphaBest = zeros(N,1); % initialise portfolio weights 

% calculate the optimal portfolio
for i=(K-1):-1:0
    if i > 0 % for each period except T
        for j = 1:length(alphaG)
            alpha = alphaG(j);
            RU = (1\(1-gamma))*((alpha*exp(path(:,i+1))+ (1-alpha)*exp(rf)).^(1-gamma)).*U_fut;
            x = [ones(N,1) s_PI_sim(:,i)];
            b = (x'*x)\x'*RU;
            CU = x*b;
            for n = 1:N
                if CU(n)>CU_Best(n)
                    CU_Best(n) = CU(n);
                    alphaBest(n) = alpha;            
                end
            end
        end
        U_fut = (alphaBest.*exp(path(:,i+1))+(1-alphaBest).*exp(rf)).^(1-gamma) .* U_fut;
        CU_Best = ones(N,1)*-inf; % Re-intialise 
        alphaBest = zeros(N,1); % Re-initialise 
    else % period T
         for j = 1:length(alphaG)
            alpha = alphaG(j);
            RU = (1\(1-gamma))*((alpha*exp(path(:,i+1))+ (1-alpha)*exp(rf)).^(1-gamma)).*U_fut;
            x = ones(N,1);
            b = (x'*x)\x'*RU;
            CU = x*b;
            for n = 1:N
                if CU(n)>CU_Best(n)
                    CU_Best(n) = CU(n);
                    alphaBest(n) = alpha;            
                end
            end
         end
    end
end

% Printing the results
disp(CU_Best(1))
disp(alphaBest(1))

%% c Buy-and-hold Strategy
% initialisation
gammar = 5;                 % CRRA coefficient 
path_PI = r_PI_sim;         % simulated path from 1b
path_DT = r_DT_sim;         % simulated path from 1c
rf = 0;                     % risk-free rate
alphaG = linspace(0,1,100); % portfolio weights grid

%% plug-in method paths
% initialisation
path = sum(path_PI,2);  % cumulative returns
CU_Best = -inf;         % initial values
alphaBest = -inf;       % initial values

% Calculate optimal portfolio
for i = 1:length(alphaG)
    alpha = alphaG(i);
    RU = ((1\(1-gamma))*(alpha*exp(path)+(1-alpha)*exp(rf)).^(1-gamma));
    CU = mean(RU);
    if (CU > CU_Best)
        CU_Best = CU;
        alphaBest = alpha;
    end
end

% Printing the results
disp(CU_Best)
disp(alphaBest)

%% decision-theory method paths
% initialisation
path = sum(path_DT,2);  % cumulative returns
CU_Best = -inf;         % initial values
alphaBest = -inf;       % initial values

% Calculate optimal portfolio
for i = 1:length(alphaG)
    alpha = alphaG(i);
    RU = (1\(1-gamma))*(alpha*exp(path)+(1-alpha)*exp(rf)).^(1-gamma);
    CU = mean(RU);
    if (CU > CU_Best)
        CU_Best = CU;
        alphaBest = alpha;
    end
end

% Printing the results
disp(CU_Best)
disp(alphaBest)





