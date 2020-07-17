% ------------------------------------------------------------------
% SINDy method for discovering slow timescale dynamics 
% ------------------------------------------------------------------
% Application to a signal of the form
%
%           x(t) = eps*P(t) + C(eps*t)     
%
% Here eps > 0 is taken to be small. The periodic component P(t)
% is given by a Fourier series with randomly generated coefficients. The
% chaotic slow component C(t) is a trajectory on the circularly symmetric
% attract of Thomas governed by the equations:
%
%           C1' = sin(C2) - b*C1
%           C2' = sin(C3) - b*C2
%           C3' = sin(C1) - b*C3
%
% Here b > 0 is a system parameter.
%
% This code is associated with the paper 
% "Sparse Identification of Slow Timescale Dynamics" by Jason J. 
% Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020).This script is 
% used to obtain the results in Section IV C.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Number of variables
n = 3;

%Model parameters for Thomas system
b = 0.2;

%Slow timescale parameter
eps = 0.1;

%ODE generation parameters
dt = 0.01;
tspan = (0:100000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

%Thomas trajectory
x0(1,:) = [0.3; 0.2; 0.1]; 
[t, xdat] = ode45(@(t,x) Thomas(x,b,eps),tspan,x0(1,:),options);
l = length(t);

% Generate Signal
T = 0.25; %Period of fast scale
N = 10; %number of Fourier modes
eta = 0;% %strength of noise to add to signal
periodic = [0, 0, 0].*ones(l,3) + (2*rand(1,3)-ones(1,3)).*ones(l,3);
for k = 1:N
    Q = 2*rand(2,3)-ones(2,3);
    periodic = periodic + [sin(k*2*pi*t/T), cos(k*2*pi*t/T)]*Q.*(ones(l,3) + eta*(2*rand(l,3) - ones(l,3)));
end
xdat = eps*periodic + xdat;

%% Coarsened Data

xt(1:n,1) = xdat(1,:)'; %Start measuring at t = 0 
count = 2;

for j = 1:l-1
   if (mod(t(j+1),T) >= 0) && (mod(t(j),T) >= T - 1*dt)  
        xt(1:n,count) = xdat(j,:)';
        count = count+1;
   end
end

xtnext = xt(:,2:end)';
xt = xt(:,1:end-1)';

%% SINDy for Slow Discovery

% Access SINDy directory
addpath Util

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 1; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = eps^2;      % lambda is our sparsification knob.

% apply iterative least squares/sparse regression
Xi = sparsifyDynamics(Theta,xtnext,lambda,n);
if n == 4
[yout, newout] = poolDataLIST({'x','y','z','w'},Xi,n,polyorder,usesine);
elseif n == 3
[yout, newout] = poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
elseif n == 2
 [yout, newout] = poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
elseif n == 1 
  [yout, newout] = poolDataLIST({'x'},Xi,n,polyorder,usesine);
end 

fprintf('SINDy model: \n ')
for k = 2:size(newout,2) 
    SINDy_eq = newout{1,k}; 
    SINDy_eq = [SINDy_eq  ' = '];
    new = 1;
   for j = 2:size(newout, 1)
       if newout{j,k} ~= 0 
           if new == 1 
             SINDy_eq = [SINDy_eq  num2str(newout{j,k}) newout{j,1} ];  
             new = 0;
           else 
             SINDy_eq = [SINDy_eq  ' + ' num2str(newout{j,k}) newout{j,1} ' '];
           end 
       end
   end
  fprintf(SINDy_eq)
  fprintf('\n ')
end 


%% Thomas system right-hand-side

function dx = Thomas(x,b,eps)

    dx = eps*[sin(x(2)) - b*x(1); sin(x(3)) - b*x(2); sin(x(1)) - b*x(3)];

end





