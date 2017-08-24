%% Solves ODEs for membrane shape over a range of tip mean curvatures using coatPinch
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2016
%
%   Last edited: 3/14/2016
%
%%

% Input(s):
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   acoat0 - dimensionless coat area
%   alphaFi - dimensionless minimium area of applied force
%   alphaFf - dimensionless maximum area of applied force
%   HpRng - range of pole mean curvatures, in units of nm^-1
%   F0 - Initial guess for the applied force, in units of pN
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the solution

% Output(s):
%   FvsHp - array of computed force vs. position
%   coatPullSol - solution array
%   compHpRng - computed range of pole z-positions, in units of nm

function [FvsHp, coatPinchSol, compHpRng] = loopCoatPinch(alpha, mesh, lambda, acoat, aFi, aFf, HpRng, F0, k0, dk, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lambda, k0, R0);        % initial guess
end

FvsHp = zeros(2, length(HpRng));    % initialize FvsHp matrix
coatPinchSol = zeros(6, length(mesh), length(HpRng));   % initialize solution matrix
compHpRng = HpRng;

% display a status bar for the calculation
h = waitbar(0,sprintf('Calculating %d/%d... Hp = %0.4f nm^{-1}, F_0 = %0.0f pN', 1, length(HpRng), HpRng(1), F0));

figure; % open a figure for the intermediate solutions

% loop over the HpRng vector
for ii = 1:length(HpRng)
    
    lastwarn('');
   
    % update the status bar
    waitbar(ii/length(HpRng), h, sprintf('Calculating %d/%d... Hp = %0.4f nm^{-1}, F_0 = %0.3f pN', ii, length(HpRng), HpRng(ii), F0))
    

    try
    
    % solve for the iith value of HpRng
    [t,Sol,F] = coatPinch(alpha, mesh, lambda, acoat, aFi, aFf, HpRng(ii), F0, k0, dk, P, gamma, C0, R0, initSol);
    
    [warnStr, warnID] = lastwarn;
    
    % catches errors from loopCoatPinch
    if strcmp(warnID, 'MATLAB:bvp4c:RelTolNotMet')
        
        error(warnStr)
        
    end
    
    catch ME
        
        display(ME.message);
        
        coatPinchSol = coatPinchSol(:,:,1:ii-1);
        
        FvsHp = FvsHp(:,1:ii-1);
        
        compHpRng = compHpRng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    coatPinchSol(:,:,ii) = Sol;
    
    % assign iith value of FvsHp
    FvsHp(:,ii) = [Sol(4,1)/R0, F];
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
    % set initial guess for force for next iteration
    F0 = F;
    
end

close(h)    % close status bar

display(sprintf('Final solution: H_p = %0.4f nm^{-1}, F = %.3f pN', compHpRng(end), F));

% plot the resultant profile of the membrane
coatArea = [0 acoat];
actArea = [aFi aFf];
plotTitle = sprintf('Membrane profile, \\lambda = %0.3f pN/nm,  P = %0.3f pN/nm^2, F = %0.2f pN', lambda, P, F);
%xLim = [-20 20];

plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], [], [], plotTitle, 0)