%% compute sensitivity indices
sensitivities;

%% FIGURES

addpath( '../helper_functions' );

xvalues = {'Glycolysis^{II}_L','Gluconeogenesis^{II}_L', 'Glycogenesis_L','Glycogenolysis_L', 'Lipolysis_A'};
yvalues = {'\uparrow 5 muscle mass','\downarrow 5 muscle mass', ...
    '\uparrow 5 fat mass','\downarrow 5 fat mass'};

fg = figure(1);
subplot(2,2,1)
colormap parula
h=heatmap(xvalues,yvalues,male_HiC_3);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-4e-3 4e-3];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,2)
colormap parula
h=heatmap(xvalues,yvalues,female_HiC_3);
%h.ColorScaling = 'scaledrows';
h.ColorLimits = [-4e-3 4e-3];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,3)
colormap parula
h=heatmap(xvalues,yvalues,male_HiF_3);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-4e-3 4e-3];
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,4)
h=heatmap(xvalues,yvalues,female_HiF_3);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-4e-3 4e-3];
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
colormap parula
set(gca,'FontSize',18,'FontName','Times New Roman')
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)'},...
    'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
    'FontSize', 18, 'FontWeight', 'normal')

fg = figure(2);
subplot(2,2,1)
colormap parula
h=heatmap(xvalues,yvalues,male_HiC_9);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.01 0.01];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,2)
colormap parula
h=heatmap(xvalues,yvalues,female_HiC_9);
%h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.01 0.01];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,3)
colormap parula
h=heatmap(xvalues,yvalues,male_HiF_9);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.01 0.01];
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,4)
h=heatmap(xvalues,yvalues,female_HiF_9);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.01 0.01];
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
colormap parula
set(gca,'FontSize',18,'FontName','Times New Roman')
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)'},...
    'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
    'FontSize', 18, 'FontWeight', 'normal')

fg = figure(3);
subplot(2,2,1)
colormap parula
h=heatmap(xvalues,yvalues,male_HiC_24);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.2 0.2];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,2)
colormap parula
h=heatmap(xvalues,yvalues,female_HiC_24);
%h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.2 0.2];
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,3)
colormap parula
h=heatmap(xvalues,yvalues,male_HiF_24);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.2 0.2];
set(gca,'FontSize',18,'FontName','Times New Roman')

subplot(2,2,4)
h=heatmap(xvalues,yvalues,female_HiF_24);
 %h.ColorScaling = 'scaledrows';
h.ColorLimits = [-0.2 0.2];
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
colormap parula
set(gca,'FontSize',18,'FontName','Times New Roman')
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)'},...
    'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
    'FontSize', 18, 'FontWeight', 'normal')
