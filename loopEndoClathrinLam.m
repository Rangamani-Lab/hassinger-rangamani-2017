%% Solves ODEs for membrane shape over a range of membrane tensions using endoClathrin
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%   Last Edited: 2/22/2016
%
%%

% Inputs:
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lamRng - membrane tensions to loop over, in units of pN/nm
%   a0 - dimensionless coat area
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference across the membrane, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat, in units of nm^-1
%   R0 - nondimensionalization length
%   initSol - initial guess for the first solution, can be input as empty
%       array

% Outputs:
%   C0Rng - returned range of coat spontaneous curvatures in case of early termination 
%   endoClathrinSol - solution array

function [lamRng, endoClathrinSol] = loopEndoClathrinLam(alpha, mesh, lamRng, a0, k0, dk, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lamRng(1), k0, R0);        % initial guess
end

endoClathrinSol = zeros(6, length(mesh), length(lamRng));   % initialize solution matrix

% display a status bar for the calculation
h = waitbar(0,sprintf('\\lambda = %0.4f, Calculating... %0.0f/%0.0f', lamRng(1), 1, length(lamRng)));

figure; % open a figure for the intermediate solutions

% loop over the lamRng vector
for ii = 1:length(lamRng)
   
    % update the status bar
    waitbar(ii/length(lamRng), h, sprintf('\\lambda = %0.4f, Calculating... %0.0f/%0.0f', lamRng(ii), ii, length(lamRng)))
    

    try
    
    % solve for the iith value of lamRng
    [~,Sol] = endoClathrin(alpha, mesh, lamRng(ii), a0, k0, dk, P, gamma, C0, R0, initSol);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        endoClathrinSol = endoClathrinSol(:,:,1:ii-1);
        
        lamRng = lamRng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    endoClathrinSol(:,:,ii) = Sol;
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

display(sprintf('Final solution: \\lambda = %0.3f', lamRng(end)));

coatArea = [0 a0];
xLim = [-sqrt(2*alpha)*R0 sqrt(2*alpha)*R0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f', a0, C0);
plotMemProfileArea(endoClathrinSol(:,:,end), t, R0, coatArea, [], [], xLim, [], plotTitle, 0);
