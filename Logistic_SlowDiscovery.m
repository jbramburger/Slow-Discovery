% ------------------------------------------------------------------
% SINDy method for discovering slow timescale dynamics 
% ------------------------------------------------------------------
% Application to the singularly perturbed non-autonomous Logistic
% equation
%
%           x' = eps*x*(1 + sin(T*t) - x)     
%
% Here eps > 0 is taken to be small and T > is real-valued.
%
% This code is associated with the paper 
% "Sparse Identification of Slow Timescale Dynamics" by Jason J. 
% Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020).This script is 
% used to obtain the results in Section IV A.
% ------------------------------------------------------------------

% Clean workspace
clear all
close all 
clc

format long

%Model parameters
eps = 0.01; % epsilon value
T = 2*pi; % period of forcing terms

%Inializations for generating trajectories
m = 2; %Dimension of ODE
n = m-1; %Dimension of Poincaré section
dt = 0.005;
tspan = (0:20000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
kfinal = 5;
if kfinal >= 1
    for k = 1:kfinal
        x0(k,:) = [4*rand; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) Logistic(x,eps,T),tspan,x0(k,:),options);
    end
end

%% Coarsened Data

%Counting parameter
count = 1;

%Known equilibrium entry
Psec(1) = 0;
PsecNext(1) = 0;

%Create Poincare section data
for i = 1:kfinal
    for j = 1:length(xdat(i,:,1))-1 
        if  (mod(xdat(i,j,2),1) >= 0.95 && mod(xdat(i,j+1,2),1) <= 0.05) %&& j >= length(xdat(i,:,1))/2) 
            temp(count) = xdat(i,j+1,1); %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1)'];
    PsecNext = [PsecNext; temp(2:length(temp))'];
   	count = 1;
    temp = [];
end

%% SINDy for Slow Discovery

% Access SINDy directory
addpath Util

% Create the recurrence data
xt = Psec;
xtnext = PsecNext;

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

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

%% Simulate SINDy Map

a = zeros(10000,1); %SINDy map solution
b = zeros(10000,1);
a(1) = 0.1;
b(1) = x0(1,1);

for k = 1:9999
   for j = 1:length(Xi)
      a(k+1) = a(k+1) + Xi(j)*a(k)^(j-1);
      b(k+1) = b(k+1) + Xi(j)*b(k)^(j-1);
   end
end

%% Plot Results

% Figure 1: Continuous-time solution
figure(1)
plot(tspan,xdat(1,:,1),'b','LineWidth',2)
set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Solution of the ODE','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Figure 2: Simulations of the discovered Poincaré map
figure(2)
plot(1:3000,a(1:3000),'b.','MarkerSize',10)
hold on
plot(1:3000,b(1:3000),'r.','MarkerSize',10)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
title('Iterates of the Discovered Poincaré Mapping','Interpreter','latex','FontSize',20,'FontWeight','Bold')

%% Logistic right-hand-side
function dx = Logistic(x,eps,T)

    dx = [eps*(x(1)*(1- x(1) + sin(T*x(2)))); 1];

end

















