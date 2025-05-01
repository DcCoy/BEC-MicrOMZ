% Script to plot tracer and diagnostic averages over final year of simulations
clear all
warning off

% Extract tracers 
fdir  = ['../output/'];
fname = ['MicrOMZ_obligate_baseline_Km.mat'];
fname = ['MicrOMZ_obligate_MM_baseline_new_rst.mat'];
fname = ['MicrOMZ_obligate_baseline_uniform_Km.mat'];
load([fdir,fname]);

cmd = ['mkdir plots/']; system(cmd);
cmd = ['mkdir plots/',fname(1:end-4)]; system(cmd);
cmd = ['mkdir plots/',fname(1:end-4),'/profiles']; system(cmd);

% List tracers, titles, ranges for plot1
tracers = {'o2','nh4','no3','no2','n2o','n2'};
ttits   = {'O$_2$','NH$^{+}_4$','NO$^{-}_3$','NO$^{-}_2$','N$_2$O','N$_2$'};

% List diagnostics, titles, ranges for plot2
diags = {};
dtits = {};

% Find index of end_year - 1
tidx = find(inputs.hist_time >= inputs.hist_time(end)-1);

% Tracer plots
fig = piofigs('lfig',1./length(tracers));
for i = 1:length(tracers)
	% Get subplot
	sb(i) = subplot(1,length(tracers),i);
end
for i = 1:length(tracers)
	set(fig,'CurrentAxes',sb(i));
	% Extract stats	
	this_mean = nanmean(bgc.([tracers{i},'_avg'])(tidx(1):tidx(end),:),1);
	this_std  = nanstd(bgc.([tracers{i},'_avg'])(tidx(1):tidx(end),:),0,1);
	low_dat   = this_mean - this_std;
	hih_dat   = this_mean + this_std;
	% Plot	
	plot(this_mean',-grid.z_r,'k','linewidth',0.5);
	hold on
	[patch1] = shaded_1d(low_dat',hih_dat',-grid.z_r,'none',rgb('DimGray'),0.7);
	sb(i).YDir = 'reverse';
	sb(i).FontSize = 6;
	if i == 1
		ylabel('Depth (m)','Interpreter','Latex');
	else
		set(sb(i),'YTickLabel',[]);
	end
	xlabel('$\mu$M','Interpreter','Latex');
	title(ttits{i},'Interpreter','Latex');
	grid on
	ylim([0 1000]);
end
export_fig('-png',['plots/',fname(1:end-4),'/profiles/tracers'],'-m5');
close all
	

