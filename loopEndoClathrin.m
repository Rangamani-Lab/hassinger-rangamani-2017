%% Solves ODEs for membrane shape over a range of coat areas using endoClathrin
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
%   a0rng - range of dimensionless coat areas to loop over
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference across the membrane, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the first solution, can be input as empty
%       array

% Outputs:
%   a0rng - returned range of coat areas in case of early termination 
%   endoClathrinSol - solution array

function [a0rng, endoClathrinSol] = loopEndoClathrin(alpha, mesh, lambda, a0rng, k0, dk, dkG, P, gamma, C0, R0, initSol)

tic

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lambda(1), k0(1), R0);        % initial guess
end

endoClathrinSol = zeros(6, length(mesh), length(a0rng));   % initialize solution matrix

% display a status bar for the calculation
h = waitbar(0,sprintf('\\alpha = %0.0f, Calculating... \\alpha_0 = %d/%0.3f', alpha, 0, max(a0rng)));

figure; % open a figure for the intermediate solutions

% loop over the a0rng vector
for ii = 1:length(a0rng)
   
    % update the status bar
    waitbar(ii/length(a0rng), h, sprintf('\\alpha = %0.0f, Calculating... \\alpha_0 = %0.3f/%0.3f', alpha, a0rng(ii), max(a0rng)))
    

    try
    
    % solve for the iith value of a0rng
    [~,Sol] = endoClathrin(alpha, mesh, lambda, a0rng(ii), k0, dk, dkG, P, gamma, C0, R0, initSol);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        endoClathrinSol = endoClathrinSol(:,:,1:ii-1);
        
        a0rng = a0rng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    endoClathrinSol(:,:,ii) = Sol;
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

display(sprintf('Final solution: a0 = %0.3f', a0rng(end)));

% plot the resultant profile of the membrane
coatArea = [0 a0rng(end)];
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
yLim = [-300 100];
%yLim = [-100 300];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f', a0rng(end), C0);
plotMemProfileArea(endoClathrinSol(:,:,end), t, R0, coatArea, [], [], xLim, yLim, plotTitle, 0);

toc