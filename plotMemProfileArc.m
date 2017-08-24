%% Plot membrane profiles
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeley
%
%   Copyright 2015
%
%%

function plotMemProfileArc(Sol, s, R0, coatArc, actArc, barArc, xLim, yLim, plotTitle, varargin)

%figure

fontsize = 48;%32;  %54;  %70;
lineWidth = 12;  %6;
axesWidth = 6;

xLabelOn = 1;
yLabelOn = 1;
xTickLabelOn = 1;
yTickLabelOn = 1;

% handle arbitrary number of inputs
if (~isempty(varargin))
    for ii = 1:length(varargin)
        switch ii

            case 1
            if (~isempty(varargin{1}))
                xLabelOn = varargin{1};
            end
            
            case 2
            if (~isempty(varargin{2}))
                yLabelOn = varargin{2};
            end
            
            case 3
            if (~isempty(varargin{3}))
                xTickLabelOn = varargin{3};
            end
            
            case 4
            if (~isempty(varargin{4}))
                yTickLabelOn = varargin{4};
            end
            
        end
    end
end

%memColor = 'black'; 
memColor = [0.9139    0.7258    0.3063];
coatColor = 'blue';
%coatColor = [0    0.4470    0.7410];
actColor = 'red';
%actColor = [0.8500    0.3250    0.0980];
barColor = 'green';
%barColor = [0.5023    0.7491    0.4776];

plot(Sol(1,:)*R0, Sol(2,:)*R0, -Sol(1,:)*R0, Sol(2,:)*R0, 'Color', memColor, 'LineWidth', lineWidth);

hold on

if ~isempty(coatArc)
    for ii = 1:size(coatArc,1)

        sMin = coatArc(ii,1);
        sMax = coatArc(ii,2);

        plot(Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, ...
            -Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, 'Color', coatColor, 'LineWidth', lineWidth);

    end
end

if ~isempty(actArc)
    for ii = 1:size(actArc,1)

        sMin = actArc(ii,1);
        sMax = actArc(ii,2);

        plot(Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, ...
            -Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, 'Color', actColor, 'LineWidth', lineWidth);

    end
end

if ~isempty(barArc)
    for ii = 1:size(barArc,1)

        sMin = barArc(ii,1);
        sMax = barArc(ii,2);

        plot(Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, ...
            -Sol(1,s>=sMin & s<=sMax)*R0, Sol(2,s>=sMin & s<=sMax)*R0, 'Color', barColor, 'LineWidth', lineWidth);

    end
end


if xLabelOn == 1
    xlabel('R (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
end

if yLabelOn == 1
    ylabel({'Z (nm)'}, 'FontSize',fontsize, 'FontName', 'Helvetica');
end

set(gca,'FontSize',fontsize-2, 'FontName', 'Helvetica', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Linewidth', axesWidth);

if xTickLabelOn == 0
    set(gca, 'XTickLabel', []);
end

if yTickLabelOn == 0
    set(gca, 'YTickLabel', []);
end

if ~isempty(plotTitle)
    title(plotTitle, 'FontSize', fontsize+4, 'FontName', 'Helvetica');
end

if ~isempty(xLim)
    xlim(xLim);
end

if ~isempty(yLim)
    ylim(yLim);
end

daspect([1 1 1]);


hold off

end