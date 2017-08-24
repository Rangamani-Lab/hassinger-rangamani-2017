%% Solves ODEs for membrane shape over a range of spontaneous curvatures using endoAgrawal
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
%   s0Rng - range of dimensionless coat arc-lengths
%   R0 - nondimensionalization length
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   initSol - initial guess for the solution, can be input as empty array

% Outputs:
%   a0rng - returned range of coat areas in case of early termination
%   s0rng - returned range of coat arc-lengths in case of early termination
%   endoClathrinSol - solution array

function [a0rng, s0rng, endoAgrawalSol] = loopEndoAgrawal(s, gamma, C0, lambda, s0rng, R0, k0, dk, dkG, P, initSol)

tic

a0rng = zeros(1,length(s0rng));

mesh = 0:1/(length(s)-1):1;

if isempty(initSol)
    initSol = endoInit(s(end)^2/2, mesh, lambda, k0, R0);        % initial guess
end
%FvsYp = zeros(2, length(a0rng));    % initialize FvsYp matrix
endoAgrawalSol = zeros(6, length(s), length(s0rng));   % initialize solution matrix

% display a status bar for the calculation
h = waitbar(0,sprintf('\\lambda = %g pN/nm, Calculating... t_0 = %d/%0.2f', lambda, 0, max(s0rng)));

%figure; % open a figure for the intermediate solutions

% loop over the yprng vector
for ii = 1:length(s0rng)
    
    lastwarn('');
    
    try
    
    % solve for the iith value of yprng
    [a0,s,Sol] = endoAgrawal(s, gamma, C0, lambda, s0rng(ii), R0, k0, dk, dkG, P, initSol);
    
    [warnStr, warnID] = lastwarn;
    
    if strcmp(warnID, 'MATLAB:bvp4c:RelTolNotMet')
        
        error(warnStr)
        
    end
    
    catch ME
        
        display(ME.identifier);
        
        endoAgrawalSol = endoAgrawalSol(:,:,1:ii-1);
        
        s0rng = s0rng(1:ii-1);
        
        a0rng = a0rng(1:ii-1);
        
        display(sprintf('Final solution: t0 = %0.3f', s0rng(end)));
        
        break;
        
    end
    
    % assign iith solution
    endoAgrawalSol(:,:,ii) = Sol;
    
    % assign iith value of FvsYp
    %FvsYp(:,ii) = [Sol(2,1), 2*Sol(5,1)];
    
    a0rng(ii) = a0;
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
    % update the status bar
    waitbar(ii/length(s0rng), h, sprintf('\\lambda = %g pN/nm, Calculating... t_0 = %0.2f/%0.2f', lambda, s0rng(ii), max(s0rng)))
    
end

close(h)    % close status bar

toc

[y, Fs] = audioread('notify.wav');
sound(y, Fs);
clear y Fs

% plot the resultant profile of the membrane
%plotMemProfileArc(Sol, s, R0, [0 s0rng(end)], [], [], [], [], sprintf('Membrane profile, \\lambda = %0.3f pN/nm, P = %0.3f pN/nm^2, s_0 = %0.3f', lambda, P, s0rng(end)))
