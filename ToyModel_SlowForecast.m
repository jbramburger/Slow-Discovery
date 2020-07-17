% ------------------------------------------------------------------
% SINDy method for forecasting slow timescale dynamics 
% ------------------------------------------------------------------
% Application to the toy model from the paper
% 'Dynamic mode decomposition for multiscale nonlinear physics' by 
% D. Dylewsky, M. Tao, and J.N Kutz (PRE, 2019). Data is loaded from
% the file toy_model_data.mat.
%
%
% This code is associated with the paper 
% "Sparse Identification of Slow Timescale Dynamics" by Jason J. 
% Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020).This script is 
% used to obtain the results in Section II.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

% Load data
load toy_model_data.mat
tspan = toy_model_data_t;
data = toy_model_data';

% Initializations
n = 2; % Number of components
l = length(tspan);
dt = tspan(2) - tspan(1);
T = 2*pi/10; %Fast scale period
if n == 2
    sample = [1 2]; % Modes to sample
else
    sample = 1:n;
end
start = 50;
xt(1:n,1) = data(sample,start); %Start measuring at t = 0 
count = 2;

%% Coarsened Data

for j = start:l-1
   if (mod(tspan(j+1),T) >= 0) && (mod(tspan(j),T) >= T - dt)  
        xt(1:n,count) = data(sample,j);
        count = count+1;
   end
end

cut = 200; %Beyond about 300 numerical error shows up in the data
xtnext = xt(:,2:cut)';
xt = xt(:,1:cut-1)';

%% SINDy for Forecasting Slow Dynamics

% Access SINDy directory
addpath Util

addpath Util

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = .001;      % lambda is our sparsification knob.

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

%% Simulate Poincare Map

sample_end = 5000;

a = zeros(sample_end,1); %SINDy map solution
b = zeros(sample_end,1);
a(1) = xt(1,1);
b(1) = xt(1,2);

for k = 1:sample_end-1
    
    % Constant terms
    a(k+1) = Xi(1,1);
    b(k+1) = Xi(1,2);
    
    %Polynomial terms
   for p = 1:polyorder
       for j = 0:p
           a(k+1) = a(k+1) + Xi(1 + j + p*(p+1)/2,1)*(a(k)^(p-j))*(b(k)^j);
           b(k+1) = b(k+1) + Xi(1 + j + p*(p+1)/2,2)*(a(k)^(p-j))*(b(k)^j);
       end
   end
   
   if usesine == 1
        a(k+1) = a(k+1) + Xi((p+1)*p/2+p+2,1)*sin(a(k)) + Xi((p+1)*p/2+p+3,1)*sin(b(k))+ Xi((p+1)*p/2+p+4,1)*cos(a(k)) + Xi((p+1)*p/2+p+5,1)*cos(b(k));
        b(k+1) = b(k+1) + Xi((p+1)*p/2+p+2,2)*sin(a(k)) + Xi((p+1)*p/2+p+3,2)*sin(b(k))+ Xi((p+1)*p/2+p+4,2)*cos(a(k)) + Xi((p+1)*p/2+p+5,2)*cos(b(k));
   end
   
   if abs(a(k)) >= 1e2
        break
   elseif abs(b(k)) >= 1e2
        break
   end
end

%% Plot Results

% Figure 1: Simulation vs training data - Component 1
figure(1)
plot(1:500,a(1:500),'b.','MarkerSize',10)
hold on
plot(xt(:,1),'k.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('Iterations','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Forecasted Dynamics and Training Data for Component 1','Interpreter','latex','FontSize',20,'FontWeight','Bold')
legend({'Forecasted Dynamics','Training Data'}, 'Interpreter','latex','FontSize',16,'Location','best')

% Figure 2: Simulation vs training data - Component 2
figure(2)
plot(1:500,b(1:500),'r.','MarkerSize',10)
hold on
plot(xt(:,2),'k.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('Iterations','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Forecasted Dynamics and Training Data for Component 2','Interpreter','latex','FontSize',20,'FontWeight','Bold')
legend({'Forecasted Dynamics','Training Data'}, 'Interpreter','latex','FontSize',16,'Location','best')



