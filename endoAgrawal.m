%% Solve the ODEs for membrane shape with curvature-generating coat, arc-length parametrization
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%   Last Edited: 8/26/2016
%
%%

% Inputs:
%   s - arc-length meshing for the domain, i.e. 0:0.01:20
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - spontaneous curvature of the coat
%   lambda - membrane tension at the boundary, in units of pN/nm
%   S0 - dimensionless coat arc-length
%   R0 - nondimensionalization length
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   initSol - initial guess for the solution

% Outputs:
%   a0 - dimensionless coat area
%   s - arc-length mesh points
%   Sol - solution array

function [a0,s,Sol] = endoAgrawal(s, gamma, C0, lambda, S0, R0, k0, dk, dkG, P, initSol)

% declare and assign global variables to be used in nested functions
global g b gs s0 iSol lam delta delkG p pw wd
g = gamma;
b = C0*R0;
s0 = S0;
lam = lambda*R0^2/k0;
gs = s;
delta = dk;
delkG = dkG;
p = P*R0^3/k0;
pw = 1*R0^3/k0;
wd = 10/R0;
iSol = initSol;

% spontaneous curvature
%c = 0.5*c0*R0*(1 - tanh(gamma*(t - alpha0)));         

% derivative of spontaneous curvature
%dc = 0.5*c0*R0*gamma*(tanh(gamma*(t - alpha0)).^2 - 1);

% bending modulus
%b = 1 + 0.5*dk*(1 - tanh(gamma*(t - alpha0)));

% initial guess structure
solinit = bvpinit(s,@mat4init);

% solver options; increasing maximum number of mesh points solver will use
options = bvpset('NMax', 100*length(s));

% solve the boundary value problem
sol = bvp4c(@mat4ode,@mat4bc,solinit,options);

% evaluate the solution on the original mesh points
Sol = deval(sol,s);

% calculate area covered by the coat
if S0 ~= 0
    a0 = trapz(s(s<=S0),Sol(1,s<=S0));
else
    a0 = 0;
end

% plot the resultant profile of the membrane
xLim = [-350 350];
yLim = [-300 100];
plotMemProfileArc(Sol, s, R0, [0 s0], [], [], xLim, yLim, sprintf('Membrane profile, \\lambda = %g pN/nm, s_0 = %0.2f', lambda, s0))
%ylim([-10*R0 0.25*R0]);


%%Define the variables
% X(1)=x;
% X(2)=y;
% X(3)=phi;
% X(4)=H;
% X(5)=L;

%%%%%the differential equations
%------------------------------------
function dXdt = mat4ode(s, X)
%parameters
global s0 g b delta p pw wd delkG

% spontaneous curvature
c = 0.5*b*(1 - tanh(g*(s - s0))); 

% derivative of spontaneous curvature
dc = 0.5*b*g*(tanh(g*(s - s0))^2 - 1);

% bending modulus
k = 1 + 0.5*(delta-1)*(1 - tanh(g*(s - s0)));

% derivative of bending modulus
dk = 0.5*(delta-1)*g*(tanh(g*(s - s0))^2 - 1);

% Derivative of Gaussian modulus
dkg = 0.5*delkG*g*(tanh(g*(s - s0))^2 - 1);

% Second derivative of Gaussian modulus
ddkg = delkG*g^2*tanh(g*(s - s0))*(1 - tanh(g*(s - s0))^2);

% opposing pressure from rigid wall
%dp = -p*(tanh(g*(X(2)))) - pw/2*(1 + tanh(g*(X(2) - wd)));

% no wall to oppose pressure
dp = p;

dXdt = [cos(X(3))
        sin(X(3))
        (2*X(1)*X(4)-sin(X(3)))/X(1)
        X(5)/X(1) + dc - dk/k*(X(4) - c)
        X(1)*((dp)/k + 2*(X(4)*((X(4)-c)^2 + X(6)/k) - (X(4)-c)*(X(4)^2 + (X(4)-sin(X(3))/X(1))^2))) - dk/k*X(5) - ddkg/k*sin(X(3)) - dkg/k*cos(X(3))*(2*X(4) - sin(X(3))/X(1))   
        -dk*(X(4) - c)^2 + 2*k*(X(4) - c)*dc - dkg*(X(4)^2 - (X(4) - sin(X(3))/X(1))^2)
        ];
            

%-------------------------boundary conditions-------------

function res = mat4bc(Xa,Xb) 
global lam
ds = 0.0001; % small offset to prevent division by 0
   
    % boundary conditions - see Hassinger et al, 2016
    res = [ Xa(1)-ds
            Xb(2)         
            Xa(3)               
            Xb(3)       
            Xa(5) 
            Xb(6) - lam
            ];
        
        
%-----------------------------------Initial guesses------------


function Xinit = mat4init(s)
 
 global gs iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Xinit = iSol(:,find(gs==s));
    
    
    