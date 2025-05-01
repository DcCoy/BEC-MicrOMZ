function plot_tracer(fname,outdir,vars,tits)
% Matlab function to plot the results from the MicrOMZ 1D model 
warning off

% Get output file
load(fname);

% Get nrec
nrec = inputs.nt_hist; 

% Get timesteps to plot
if nrec > 500
	dt = 10;
elseif nrec > 100
    dt = 5;
else
    dt = 1;
end

% Get colormap for old timesteps
old_cmap = cmocean('balance',(nrec/dt)+(floor((nrec/dt).*0.2)));
old_cmap = old_cmap(floor((nrec/dt).*0.1):end-floor((nrec/dt).*0.1),:);

% Make directories, remove existing plots
cmd = ['mkdir plots']; system(cmd);
cmd = ['mkdir plots/',outdir]; system(cmd);
cmd = ['rm plots/',outdir,'/tracer*png']; system(cmd);

% Get max bacterial biomass (for limits)
for i = 1:length(vars)
    tmp   = squeeze(bgc.(vars{i}));
    lims{i} = [0 max(max(tmp(:)),1e-24)];
end

% Initialize figure
fig = piofigs('lfig',(1/length(vars)));
for i = 1:length(vars)
	sb(i)  = subplot(1,length(vars),i);
	sb(i).Position(2) = sb(i).Position(2) + 0.15;
	sb(i).Position(4) = sb(i).Position(4) - 0.15;
	pos{i} = sb(i).Position;
end

% Get data up until tstep 
for i = 1:length(vars)
	idx = find(strcmp(vars{i},opt.tracers)==1);
	set(fig,'CurrentAxes',sb(i));
    plot(NaN,NaN,'linestyle','none','marker','none','color','none');
    l = legend(tits{i},'Location','Southeast');
    l.Interpreter = 'Latex';
    l.AutoUpdate = 'off';
    l.Box = 'off';
	this_var = squeeze(bgc.(vars{i}));
	hold on
	tcnt = 1;
	for t = 1:dt:nrec
		if t == 1
			plot(bgc.(vars{i})(t,:),-grid.z_r,'--','linewidth',0.5,'color',rgb('Black'));
		else
			plot(bgc.(vars{i})(t,:),-grid.z_r,'-','linewidth',0.5,'color',old_cmap(tcnt,:));
			tcnt = tcnt + 1;
		end
	end
	plot(bgc.(vars{i})(end,:),-grid.z_r,'-','linewidth',0.5,'color',rgb('Black'));
	xlim([lims{i}]);
	ylim([0 max(abs(grid.z_r))]);
	set(gca,'YDir','Reverse');
	if i > 1
		set(gca,'YTickLabel',[]);
	else
		ylabel('Depth ($m$)','Interpreter','Latex');
	end
	set(gca,'FontSize',6);
	grid on
	sb(i).Position = pos{i};
end
	
export_fig('-png',['plots/',outdir,'/tracers'],'-m5','-p0.05');
close all
        
