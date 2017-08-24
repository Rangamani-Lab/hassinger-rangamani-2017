%% Solve the ODEs for membrane shape with curvature-generating coat, Gaussian 
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2016
%
%   Last Edited: 8/27/2016
%
%%

% Inputs:
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   alpha0 - dimensionless coat area
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the solution

% Outputs:
%   t - area mesh points
%   Sol - solution array


function [t, Sol] = endoClathrin(alpha, mesh, lambda, alpha0, k0, dk, dkG, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

% declare and assign global variables to be used in nested functions
global a gt a0 iSol lam c0 g delta delkG p pw wd

a = alpha;
a0 = alpha0;
gt = t;
lam = lambda*R0^2/k0;   % dimensionless membrane tension
c0 = R0*C0;             % dimensionless preferred curvature
g = gamma;
delta = dk;
delkG = dkG;
p = P*R0^3/k0;          % dimensionless pressure
pw = 1*R0^3/k0;         % dimensionless pressure from rigid wall
wd = 10/R0;             % dimensionless wall "distance"
iSol = initSol;

% initial guess structure
solinit = bvpinit(t,@mat4init);

% solver options; increasing maximum number of mesh points solver will use
options = bvpset('NMax', 100*length(t), 'RelTol', 1e-3);

% solve the boundary value problem
sol = bvp4c(@mat4ode,@mat4bc,solinit,options);

% evaluate the solution on the original mesh points
Sol = deval(sol,t);

% plot the resultant profile of the membrane
coatArea = [0 alpha0];
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
%xLim = [-350 350];
yLim = [-300 100];
%yLim = [-100 300];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f', alpha0, C0);
plotMemProfileArea(Sol, t, R0, coatArea, [], [], xLim, yLim, plotTitle, 0);


%%Define the variables
% X(1)=x;
% X(2)=y;
% X(3)=phi;
% X(4)=H;
% X(5)=L;

%%%%%the differential equations
%------------------------------------
function dXdt = mat4ode(t, X)
%parameters
global c0 a0 g delta delkG p pw wd

% spontaneous curvature
c = 0.5*c0*(1 - tanh(g*(t - a0))); 

% derivative of spontaneous curvature
dc = 0.5*c0*g*(tanh(g*(t - a0))^2 - 1);

% bending modulus
b = 1 + 0.5*(delta-1)*(1 - tanh(g*(t - a0)));

% derivative of bending modulus
db = 0.5*(delta-1)*g*(tanh(g*(t - a0))^2 - 1);

% Derivative of Gaussian modulus
dkg = 0.5*delkG*g*(tanh(g*(t - a0))^2 - 1);

% Second derivative of Gaussian modulus
ddkg = delkG*g^2*tanh(g*(t - a0))*(1 - tanh(g*(t - a0))^2);

% opposing pressure from rigid wall
%dp = -p*(tanh(g*(X(2)))) - pw/2*(1 + tanh(g*(X(2) - wd)));

% no wall to oppose pressure
dp = p;

% ODEs - see Hassinger et al, 2016
dXdt = [cos(X(3))/X(1)
        sin(X(3))/X(1)
        (2*X(1)*X(4)-sin(X(3)))/X(1)^2
        X(5)/X(1)^2 + dc - db/b*(X(4) - c)
        dp/b + 2*X(4)*((X(4)-c)^2 + X(6)/b) - 2*(X(4)-c)*(X(4)^2 + (X(4)-sin(X(3))/X(1))^2) - db/b*X(5) - ddkg/b*X(1)*sin(X(3)) - dkg/b*cos(X(3))*(2*X(4) - sin(X(3))/X(1))
        -db*(X(4) - c)^2 + 2*b*(X(4) - c)*dc - dkg*(X(4)^2 - (X(4) - sin(X(3))/X(1))^2)
        ];
            

%-------------------------boundary conditions-------------

function res = mat4bc(Xa,Xb) 
global lam
ds = 1e-4;  % small offset to prevent division by 0
   
    % boundary conditions - see Hassinger et al, 2016
    res = [ Xa(1) - ds
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
    
    
    