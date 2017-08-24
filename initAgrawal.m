%% Initialize flat patch
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%%

% Initializes a flat, circular patch of membrane for use with the
% axisymmetric Helfrich membrane model

% t = mesh points in terms of area
function initSol = initAgrawal(t, lambda, R0, k0)

ds=0.0001;    % Small deviation from x=0 (necessary to avoid division by 0)

lam = lambda*R0^2/k0;

initSol = [ds + t                   % x
           zeros(1, length(t))      % y
           zeros(1, length(t))      % psi
           zeros(1, length(t))      % h
           zeros(1, length(t))      % l  
           lam*ones(1, length(t))]; % lambda