function main

clear all
clc

%Define 
l = 4;                                  
dx = 0.25;
Ts = 0.01;
T = 0.1;
X = meshgrid(-l+dx/2:dx:l-dx/2);
Y = X.';
N = length(X);

rho0 = zeros(N,N);
nu0 = zeros(N,N);
rho0(X<3 & X>-3 & Y<3 & Y>-3) = 1/4;

% vectorize the quantities
rho0 = rho0(:);
nu0 = nu0(:);

rhoR = discrete_evol (rho0, nu0, Ts);
% imagesc(reshape(rho0,N,N));
% imagesc(reshape(rhoR,N,N));

% mpc controller
nx = length(rho0);
ny=nx; nu=nx;
nlobj = nlmpc(nx,ny,nu); % as I understand we have the same number of state variables as output variables and control variables
nlobj.PredictionHorizon = 10;
nlobj.ControlHorizon = 5;
nlobj.Model.StateFcn = "discrete_evol";
nlobj.Model.IsContinuousTime = false;
nlobj.Model.NumberOfParameters = 1;
nlobj.Model.OutputFcn = 'OC_transportOutputFcn';

validateFcns(nlobj,rho0,nu0,[],{Ts});

nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};

% initialize the states
x = rho0;
y = x;

% reference to track, i.e. move the initial square by one unit on x-axis
yref = zeros(nx,nx);
yref(X<4 & X>-2 & Y<3 & Y>-3) = 1/4;
mv = zeros(nx,1);

% Run the simulation for |20| seconds.
Duration = 20;
hbar = waitbar(0,'Simulation Progress');
xHistory = x;
for ct = 1:(20/Ts)
    [mv,nloptions,info] = nlmpcmove(nlobj,x,mv,yref,[],nloptions);
    % Implement first optimal control move and update plant states.
    x = discrete_evol(x,mv,Ts);
    % Generate sensor data with some white noise.
    % y = x([1 3]) + randn(nx,1)*0.01; 
    % Save plant states for display.
    xHistory = [xHistory x]; %#ok<*AGROW>
    waitbar(ct*Ts/20,hbar);
end
close(hbar)

end