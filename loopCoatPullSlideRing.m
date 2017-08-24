%% Loops over a range of pole heights to extract force vs. displacement curve
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%%

% finalSol = solution for final pole position
%   x = Sol(1,:), y = Sol(2,:), psi = Sol(3,:), h = Sol(4,:), l = Sol(5,:)
% FvsYp = matrix of values for force vs. displacement plot
%   yp =FvsYp(1,:), f = FvsYp(2,:)
% alpha = dimensionless patch area
% meshPt = number of mesh points, 1000 typically sufficient
% yprng = vector of pole postion values to solve for, i.e. yprng = 0:0.1:20
function [FvsZp, coatPullSol, compZpRng] = loopCoatPullSlideRing(alpha, mesh, lambda, acoat, rF, rIn, rOut, zpRng, f0, k0, dk, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lambda, k0, R0);        % initial guess
end

FvsZp = zeros(2, length(zpRng));    % initialize FvsYp matrix
coatPullSol = zeros(6, length(mesh), length(zpRng));   % initialize solution matrix
compZpRng = zpRng;

% display a status bar for the calculation
h = waitbar(0,sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', zpRng(1), zpRng(end)));

figure; % open a figure for the intermediate solutions

% loop over the yprng vector
for ii = 1:length(zpRng)
    
    lastwarn('');
   
    % update the status bar
    waitbar(ii/length(zpRng), h, sprintf('Calculating... z_p = %0.2f nm/%0.2f nm', zpRng(ii), zpRng(end)))
    

    try
    
    % solve for the iith value of yprng
    [~,Sol,f] = coatPullSlideRing(alpha, mesh, lambda, acoat, rF, rIn, rOut, zpRng(ii), f0, k0, dk, P, gamma, C0, R0, initSol);
    
    [warnStr, warnID] = lastwarn;
    
    if strcmp(warnID, 'MATLAB:bvp4c:RelTolNotMet')
        
        error(warnStr)
        
    end
    
    catch ME
        
        display(ME.message);
        
        coatPullSol = coatPullSol(:,:,1:ii-1);
        
        FvsZp = FvsZp(:,1:ii-1);
        
        compZpRng = compZpRng(1:ii-1);
        
        break;
        
    end
    
    % assign iith solution
    coatPullSol(:,:,ii) = Sol;
    
    % assign iith value of FvsYp
    FvsZp(:,ii) = [Sol(2,1)*R0, f];
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

display(sprintf('Final solution: z_p = %0.2f nm', compZpRng(end)));

% plot the resultant profile of the membrane
coatArea = [0 acoat];
%actArea = [0 ((rF/R0)^2)/2; aIn aOut];
plotTitle = sprintf('Membrane profile, \\lambda = %0.3f pN/nm,  P = %0.3f pN/nm^2, A_{coat} = %0.3f nm^2', lambda, P, acoat*2*pi*R0^2);
%xLim = [-20 20];

%plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], [], [], plotTitle, 1, 1, 1, 1)
