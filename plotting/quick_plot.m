function [fig] = quick_plot(VARS,BGC,DIAG,GRID);
% Simple function to return a figure based on inputs
%
% Usage:
% - fig = quick_plot(VARS,BGC,DIAG,GRID)
%
% Inputs:
% - VARS = cell array of bgc.* or diag.* variables
% - BGC  = output structure from model run (bgc)
% - DIAG = output structure from model run (diag)
% - GRID = grid structure from output
close all;

% Initiate figure
if length(VARS)==1
	fig = piofigs('sfig',1);
	sb(1) = subplot(1,1,1);
else
	fig = piofigs('lfig',1./length(VARS));
	for i = 1:length(VARS)
		sb(i) = subplot(1,length(VARS),i);
	end
end

% Get cmap
try
	cmap = cmocean('balance',length(BGC.(VARS{1})));
catch
	cmap = cmocean('balance',length(DIAG.(VARS{1})));
end

% Make figure
for i = 1:length(VARS)
	set(fig,'CurrentAxes',sb(i));
	try
		this_var = BGC.(VARS{i});
	catch
		this_var = DIAG.(VARS{i});
	end
	for j = 1:length(cmap)
		plot(this_var(j,:),-GRID.z_r,'color',cmap(j,:),'linewidth',0.1);
		hold on
	end
    ylim([0 max(abs(GRID.z_r))]);
    set(gca,'YDir','Reverse');
	if i == 1
		ylabel('Depth ($m$)','Interpreter','Latex');
	else
		set(gca,'YTickLabel',[]);
	end
	title(VARS{i},'Interpreter','none');
    set(gca,'FontSize',6);
    grid on
end

pltjpg(1);
