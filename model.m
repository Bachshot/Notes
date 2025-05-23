%% File Info.

%{

    model.m
    -------
    This code sets up the model.

%}

%% Model class.

classdef model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
        function par = setup()            
            %% Structure array for model parameters.
            
            par = struct();
            
            %% Preferences.
            
            par.beta = 0.96; % Discount factor: Lower values of this mean that consumers are impatient and consume more today.
            par.sigma = 2.00; % CRRA: Higher values of this mean that consumers are risk averse and do not want to consume too much today.
            
            par.gamma = 1.00; % Weight on leisure: Higher values mean that leisure has a higher weight in the utility function.
            par.nu = 0.04; % Frisch Elasticity: Higher values of this mean that the labor choice becomes more sensitive to productivity shocks.
            % Tax Parameters
            par.lambda = 0.2; % Tax level
            par.tau = 0.2; % Progressivity
            par.chi = 0.5;

            par.slen = 2;                               %% updated
            par.prob_s = [0.5, 0.5];                    %% updated

            wl = 2;                                   %% updated
            wh = 5;                                   %% updated
            par.w = [wl, wh];                           %% updated
            par.skill_names = {'Low skill', 'High skill'}; %% optional for display

            assert(par.beta > 0 && par.beta < 1.00,'Discount factor should be between 0 and 1.\n')
            assert(par.sigma > 0,'CRRA should be at least 0.\n')
            assert(par.gamma > 0,'The weight on leisure should be at least 0.\n')
            assert(par.nu > 0,'The sensitivity to productivity shocks should be positive.\n')

            %% Technology.

            par.alpha = 0.33; % Capital's share of income.
            par.delta = 0.05; % Depreciation rate of physical capital.

            assert(par.alpha > 0 && par.alpha < 1.00,'Capital share of income should be between 0 and 1.\n')
            assert(par.delta >= 0 && par.delta <= 1.00,'The depreciation rate should be from 0 to 1.\n')

            par.sigma_eps = 0.07; % Std. dev of productivity.
            par.rho = 0.90; % Persistence of AR(1) process.
            par.mu = 0.0; % Intercept of AR(1) process.

            assert(par.sigma_eps > 0,'The standard deviation of the shock must be positive.\n')
            assert(abs(par.rho) < 1,'The persistence must be less than 1 in absolute value so that the series is stationary.\n')

            %% Simulation parameters.

            par.seed = 2025; % Seed for simulation.
            par.T = 100; % Number of time periods.

        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Capital grid.

            par.kss = (par.alpha/((1/par.beta)-(1-par.delta)))^(1/(1-par.alpha)); % Steady state capital in the deterministic case.
             
            par.klen = 300; % Grid size for k.
            par.kmax = 0.9*par.kss; % Upper bound for k.
            par.kmin = 0.01*par.kss; % Minimum k.
            
            assert(par.klen > 5,'Grid size for k should be positive and greater than 5.\n')
            assert(par.kmax > par.kmin,'Minimum k should be less than maximum value.\n')
            
            par.kgrid = linspace(par.kmin,par.kmax,par.klen)'; % Equally spaced, linear grid for k and k'.

            %% Discretized productivity process.
                  
            par.Alen = 7; % Grid size for A.
            par.m = 3; % Scaling parameter for Tauchen.
            
            assert(par.Alen > 3,'Grid size for A should be positive and greater than 3.\n')
            assert(par.m > 0,'Scaling parameter for Tauchen should be positive.\n')
            
            [Agrid,pmat] = model.tauchen(par.mu,par.rho,par.sigma_eps,par.Alen,par.m); % Tauchen's Method to discretize the AR(1) process for log productivity.
            par.Agrid = exp(Agrid); % The AR(1) is in logs so exponentiate it to get A.
            par.pmat = pmat; % Transition matrix.

        end
        
        %% Tauchen's Method
        
        function [y,pi] = tauchen(mu,rho,sigma,N,m)
            %% Construct equally spaced grid.
        
            ar_mean = mu/(1-rho); % The mean of a stationary AR(1) process is mu/(1-rho).
            ar_sd = sigma/((1-rho^2)^(1/2)); % The std. dev of a stationary AR(1) process is sigma/sqrt(1-rho^2)
            
            y1 = ar_mean-(m*ar_sd); % Smallest grid point is the mean of the AR(1) process minus m*std.dev of AR(1) process.
            yn = ar_mean+(m*ar_sd); % Largest grid point is the mean of the AR(1) process plus m*std.dev of AR(1) process.
            
	        y = linspace(y1,yn,N); % Equally spaced grid.
            d = y(2)-y(1); % Step size.
	        
	        %% Compute transition probability matrix from state j (row) to k (column).
        
            ymatk = repmat(y,N,1); % States next period.
            ymatj = mu+rho*ymatk'; % States this period.
        
	        pi = normcdf(ymatk,ymatj-(d/2),sigma) - normcdf(ymatk,ymatj+(d/2),sigma); % Transition probabilities to state 2, ..., N-1.
	        pi(:,1) = normcdf(y(1),mu+rho*y-(d/2),sigma); % Transition probabilities to state 1.
	        pi(:,N) = 1 - normcdf(y(N),mu+rho*y+(d/2),sigma); % Transition probabilities to state N.
	        
        end
        
        %% Utility function.
        
             function u = utility(c, n, par, g)
                        %% CRRA utility with public goods
            
                        if nargin < 4
                            g = 1;
                        end
            
                        if any(c <= 0, 'all') || any(g <= 0, 'all')
                            u = -Inf;
                            return
                        end
            
                        un = ((1 - n).^(1 + 1 / par.nu)) ./ (1 + 1 / par.nu);
            
                        if par.sigma == 1
                            uc = log(c);
                            ug = log(g);
                        else
                            uc = (c.^(1 - par.sigma)) ./ (1 - par.sigma);
                            ug = (g.^(1 - par.sigma)) ./ (1 - par.sigma);
                        end
            
                        u = uc + par.gamma * un + par.chi * ug;
                    end
            %% Tax function: T(y) = y - lambda * y^(1 - tau)

            function T = tax(y, par) %% new
                T = y - par.lambda * y.^(1 - par.tau);
        end
        
    end
end