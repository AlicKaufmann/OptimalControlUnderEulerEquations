
function [rhoR] = discrete_evol (rho0,nu0, dt)
%SINGLE2D: numeric simulation of the dynamics of one single species in 2D
%   r = single2d(r0,l, W, dt, T) 
%   Input:
%       r0:         initial density, NxN matrix
%       l:          domain[-l,l]x[-l,l]
%       W:          W(x) interacting potentials
%       dt:         time step
%       T:          simulation time. Total #iterations = T/dt
%   Output: 
%       r:          density at t = T, NxN matrix
%
%   Optional parm:
%   r = single2d(.. ,H): H is a symbolic function for internal energy as a
%   function of the density. Default H(r) = 0
%
%   r = single2d(.. ,V) optionally sets the environmental confinement
%   potential V, which is a NxN matrix. Default: 0
%
%   r = single2d(.. ,e) sets the diffusion coefficient for some e > 0.
%   Default e = 0
%
%   r = single2d(.. ,'v') or r = single2d(.. ,'V') enables visual display
%   during the simulation. Default disabled.
%   
%   r = single2d(.. ,'solver') where 'solver' sets the numeric method used
%   for ODE. Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.

% Matrixify the data
rho0 = reshape(rho0,int64(sqrt(length(rho0))),int64(sqrt(length(rho0))));
nu0 = reshape(nu0,int64(sqrt(length(nu0))),int64(sqrt(length(nu0))));

N = length(rho0);


% Optional input defaults
e = 0;                      % No diffusion term
v = false;                  % No visible update during simulation
V = zeros(N, N);            % No confinement potential
s = 's';

% mesh parameters
l=4;
dx = 2*l/N;

% self and cross-interaction kernels
syms K(x,y) H(x,y)
K(x,y) = -exp(- x^2 - y^2)/pi; % potential describing self-interaction population-population
H(x,y) = -exp(- x^2 - y^2)/pi; % potential describing cross interaction population-agents
K = getKernal2(K,dx,l);
H = getKernal2(H,dx,l);

L = @(X) drho_dt(X,nu0,K,H,dx, e); % X is the current rho and Y the current nu
rhoR = rho0;
X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';

%compute the next step with time sample dt
if s=='s'
    rhoR = SSPRK3(rhoR,L,dt);
else
    rhoR = Euler(rhoR,L,dt);
end 
% draw(X,Y,rhoR);
% title('after one time step');
% drawnow;

% vectorize rhoR
rhoR=rhoR(:);
end

function [K] = getKernal2 (kernel, dx, l)
X = meshgrid(-l+dx/2:dx:l-dx/2);
Y = X.';
K=double(kernel(X,Y));
end







