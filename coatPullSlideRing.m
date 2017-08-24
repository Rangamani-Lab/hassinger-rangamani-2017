%% Solve the ODEs for membrane with curvature-generating coat and pulling force
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%%

% Solves system of ODEs for axisymmetric membrane vesicle with set pole
% height => use to recover applied force

% t = mesh points in terms of area
% Sol = solution to equations as a matrix of values at each mesh point
%   x = Sol(1,:), y = Sol(2,:), psi = Sol(3,:), h = Sol(4,:), l = Sol(5,:)
% alpha = dimensionless patch area
% mesh = initial solver mesh (set values between 0 and 1), i.e. mesh = 0:1/1000:1
% initSol = initial guess, from memTubeInit.m
% pole = pole position
function [t,Sol,f] = coatPullSlideRing(alpha, mesh, lambda, alpha0, rF, rIn, rOut, zp, f0, k0, dk, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

% declare and assign global variables to be used in nested functions
global a gt a0 iSol lam c0 g delta p aF aIn aOut yp

a = alpha;
a0 = alpha0;
gt = t;
lam = lambda*R0^2/k0;
c0 = R0*C0;
g = gamma;
delta = dk;
p = P*R0^3/k0;
aF = ((rF/R0)^2)/2;
aIn = ((rIn/R0)^2)/2;
aOut = ((rOut/R0)^2)/2;
yp = zp/R0;
iSol = initSol;

beta = f0*R0^3/k0;

% initial guess structure
solinit = bvpinit(t,@mat4init, beta);

% solver options; increasing maximum number of mesh points solver will use
options = bvpset('NMax', 100*length(t), 'RelTol', 1e-3);

% solve the boundary value problem
sol = bvp4c(@mat4ode,@mat4bc,solinit,options);

% extract force
f = sol.parameters*k0/R0^3*2*pi*R0^2;

% evaluate the solution on the original mesh points
Sol = deval(sol,t);

% plot the resultant profile of the membrane
coatArea = [0 alpha0];
actArea = [0 aF; aIn aOut];
plotTitle = sprintf('Membrane profile, f = %0.3f pN', f);
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
yLim = [-500 0];
plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], xLim, yLim, plotTitle, 0)


%%Define the variables
% X(1)=x;
% X(2)=y;
% X(3)=phi;
% X(4)=H;
% X(5)=L;

%%%%%the differential equations
%------------------------------------
function dXdt = mat4ode(t, X, beta)
%parameters
global c0 a0 g delta p aF aIn aOut

% spontaneous curvature
c = 0.5*c0*(1 - tanh(g*(t - a0))); 

% derivative of spontaneous curvature
dc = 0.5*c0*g*(tanh(g*(t - a0))^2 - 1);

% bending modulus
b = 1 + 0.5*(delta-1)*(1 - tanh(g*(t - a0)));

% opposing pressure from rigid wall
dp = -p*(tanh(g*(X(2))));

% no wall to oppose pressure
%dp = p; 

% applied force
fbar = beta*(0.5*((1 - tanh(g*(t - aF)))/aF - (tanh(g*(t - aIn)) - tanh(g*(t - aOut)))/(aOut-aIn)));

dXdt = [cos(X(3))/X(1)
        sin(X(3))/X(1)
        (2*X(1)*X(4) - sin(X(3)))/X(1)^2
        X(5)/X(1)^2 + dc
        dp/b + fbar*cos(X(3)) + 2*X(4)*((X(4)-c)^2 + X(6)/b) - 2*(X(4)-c)*(X(4)^2 + (X(4)-sin(X(3))/X(1))^2)
        2*b*(X(4) - c)*dc - fbar*sin(X(3))/X(1)
        ];
            

%-------------------------boundary conditions-------------

function res = mat4bc(Xa,Xb,beta) % y position at pole specified
global lam yp
ds = 1e-4;
   
    res = [ Xa(1) - ds
            Xa(2) - yp
            Xb(2)         
            Xa(3)               
            Xb(3)             
            Xa(5) 
            Xb(6) - lam
            ];
        
        
%-----------------------------------Initial guesses------------


function Xinit = mat4init(t)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Xinit = iSol(:,find(gt==t));