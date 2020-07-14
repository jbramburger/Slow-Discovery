% ------------------------------------------------------------------
% SINDy method for discovering slow timescale dynamics 
% ------------------------------------------------------------------
% Application to the motion of Saturn as governed by a 3-body system
% with the Sun, Jupiter, and Saturn. Data is restricted to the orbital 
% plane and loaded from the file three_body_data.mat.
%
% This code is associated with the paper 
% "Sparse Identification of Slow Timescale Dynamics" by Jason J. 
% Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020).This script is 
% used to obtain the results in Section III B.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

% Load Saturn data
load saturn_data;
data = saturn_data; 

% Initializations
n = 2; %number of components
l = length(tspan);
dt = tspan(2) - tspan(1);
T = 29.5; % Saturn fast period
sample = [1 2]; % restrict to orbital plane
count = 1;

%% Coarsened Data

for j = 1:l-1
   if mod(tspan(j+1),T) == 0 
        xt(1:n,count) = data(sample,j);
        count = count+1;
   end
end

cut = 1500; % Truncate size of input data 
xtnext = xt(:,2:cut)';
xt = xt(:,1:cut-1)';

%% SINDy for Slow Discovery

% Access SINDy directory
addpath Util

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 1; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 0.001;      % lambda is our sparsification knob.

% apply iterative least squares/sparse regression
Xi = sparsifyDynamics2(Theta,xtnext,lambda,n);
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
