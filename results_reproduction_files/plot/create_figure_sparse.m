function create_figure_sparse(YMatrix1, dataset)
%CREATEFIGURE(YMatrix1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 20-Sep-2022 18:00:36

% Create figure
figure1 = figure('OuterPosition',[852 410 600 413]);

% Create axes
axes1 = axes('Position',...
    [0.0535077288941736 0.0714285714285714 0.919095010831854 0.874404761904762]);
hold(axes1,'on');


% Create multiple lines using matrix input to plot
plot1 = stem(YMatrix1);
set(plot1(1),'LineWidth',1, 'Color',[0 0 0]);
set(plot1(2),'LineWidth',1,'LineStyle','-.',...
    'Color',[0 0.447058823529412 0.741176470588235], 'Marker', 'x');
% set(plot1(3),'LineWidth',2,'Color',[0 0 0]);
% set(plot1(4),'LineWidth',1,'Color',[0 0.447058823529412 0.741176470588235]);
switch dataset
    case 'A'
xlim([0,200])
ylim([0,25])
    case 'B'
       xlim([0,200])
ylim([0,30]) 
end
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 220]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 7]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create axes


