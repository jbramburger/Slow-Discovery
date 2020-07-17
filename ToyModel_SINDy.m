% ------------------------------------------------------------------
% SINDy method for model discovery 
% ------------------------------------------------------------------
% Application to the toy model from the paper
% 'Dynamic mode decomposition for multiscale nonlinear physics' by 
% D. Dylewsky, M. Tao, and J.N Kutz (PRE, 2019). Training data is loaded 
% from the file toy_model_data.mat.
%
% This code is associated with the paper 
% "Sparse Identification of Slow Timescale Dynamics" by Jason J. 
% Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020).This script is 
% used to obtain the results in Section III.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

% Load data
load toy_model_data.mat
tspan = toy_model_data_t;
xdat = toy_model_data';

% Initializations
n = 4; % Number of components
l = length(tspan);
dt = tspan(2) - tspan(1);

%% SINDy for Continuous-Time Model Discovery

% Access SINDy directory
addpath Util

cut = 2.5e5; %same sample size as the forecast model
xt = xdat(:,2:cut-1);
dxt = (xt(:,2:end)-xt(:,1:end-1))/dt;

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 0.001;      % lambda is our sparsification knob.

% apply iterative least squares/sparse regression
Xi = sparsifyDynamics(Theta,dxt,lambda,n);
if n == 4
[yout, newout] = poolDataLIST({'x','y','z','w'},Xi,n,polyorder,usesine);
elseif n == 3
[yout, newout] = poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
elseif n == 2
 [yout, newout] = poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
elseif n == 1 
  [yout, newout] = poolDataLIST({'x'},Xi,n,polyorder,usesine);
end 
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

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







