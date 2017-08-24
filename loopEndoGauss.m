%% Solves ODEs for membrane shape over a range of spontaneous curvatures using endoClathrin
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%   Last Edited: 8/27/2016
%
%%

% Inputs:
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   a0 - dimensionless coat area
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference across the membrane, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0Rng - preferred curvatures of the coat to loop over, in units of nm^-1
%   R0 - nondimensionalization length
%   initSol - initial guess for the first solution, can be input as empty
%       array

% Outputs:
%   C0Rng - returned range of coat spontaneous curvatures in case of early termination 
%   endoClathrinSol - solution array

function [dkGrng, endoGaussSol] = loopEndoGauss(alpha, mesh, lambda, a0, k0, dk, dkGrng, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lambda, k0, R0);        % initial guess
end

endoGaussSol = zeros(6, length(mesh), length(dkGrng));   % initialize solution matrix

% display a status bar for the calculation
h = waitbar(0,sprintf('\\alpha = %0.0f, Calculating... dkG = %d/%0.4f', alpha, dkGrng(1), dkGrng(end)));

figure; % open a figure for the intermediate solutions

% loop over the C0Rng vector
for ii = 1:length(dkGrng)
   
    % update the status bar
    waitbar(ii/length(dkGrng), h, sprintf('\\alpha = %0.0f, Calculating... dkG = %0.4f/%0.4f', alpha, dkGrng(ii), dkGrng(end)))
    

    try
    
    % solve for the iith value of C0Rng
    [~,Sol] = endoClathrin(alpha, mesh, lambda, a0, k0, dk, dkGrng(ii), P, gamma, C0, R0, initSol);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        endoGaussSol = endoGaussSol(:,:,1:ii-1);
        
        dkGrng = dkGrng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    endoGaussSol(:,:,ii) = Sol;
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

display(sprintf('Final solution: dkG = %0.4f', dkGrng(end)));

% plot the resultant profile of the membrane
plotTitle = sprintf('Membrane profile, \\lambda = %0.3f pN/nm,  P = %0.3f pN/nm^2, dkG = %0.3f', lambda, P, dkGrng(end));
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
yLim = [-200 50];
plotMemProfileArea(Sol(:,:,end), t, R0, [0 a0], [], [], xLim, yLim, plotTitle, 0)
