% ------------------------------------------------------------------
% Script for generating the toy model data 
% ------------------------------------------------------------------
% Model under investigation is the toy model from the paper
% 'Dynamic mode decomposition for multiscale nonlinear physics' by 
% D. Dylewsky, M. Tao, and J.N Kutz (PRE, 2019). Out put is the file
% toy_model_data.mat to be used by the scripts ToyModel_SINDy.m and
% ToyModel_Forecast.m.
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

%% Initializations

% Initial Conditions
x1_0 = 0;
x2_0 = 0.5;
y1_0 = 0;
y2_0 = 0.5;
x0 = [x1_0; x2_0; y1_0; y2_0];

% System parameters
epsilon=0.01;
delta = 4;

%% RK4 integration of the mixed-scale system

% ODE parameters
T=2^7;
h=epsilon/100;
TimeSpan = 0:h:T;   
TimeSteps=length(TimeSpan)-1;

x = zeros(4,TimeSteps+1);
x(:,1) = [x1_0; x2_0; y1_0; y2_0];

% Numerical integration
[t, x] = ode23(@(t,x) rhs(x,delta,epsilon),TimeSpan,x0);

tStep = mean(diff(t))*4;
nSteps = ceil(T/tStep);
tN = 0:tStep:T;
tN = tN(1:nSteps); %match to nSteps

xOld = x;
x = interp1(t,xOld,tN); %enforce evenly spaced time steps
TimeSpan = tN;

%% Randomized Linear Mixing

rng(111); %seed random
n = size(x,2);
 
% Generate a random unitary matrix M
X = rand(n)/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
M = real(Q*R);
x = x*M;

%% Save generated data

save('toy_model_data.mat','x','TimeSpan')

%% Toy model right-hand-side

function dx = rhs(x,delta,epsilon)
    x1 = x(1); x2 = x(2); y1 = x(3); y2 = x(4);
    dx=[x2; -y1^2 * x1^3; y2; -epsilon^(-1)*y1 - delta*y1^3];
end


