% List fname
clear all
warning off
fdir  = ['../output/'];
fname = ['MicrOMZ_facultative-28-Apr-2024.mat'];
load([fdir,fname]);

% Make directory
mkdir plots
cmd = ['mkdir plots/',fname(1:end-4),'/']; system(cmd);
cmd = ['mkdir plots/',fname(1:end-4),'/individual']; system(cmd);
outdir = ['plots/',fname(1:end-4),'/individual/'];

% Get number of history records
nrec = inputs.nt_hist;

% Get colormap
cmap = cmocean('balance',nrec);

% Initialize figure
fig = piofigs('sfig',0.5);
sb(1) = subplot(1,2,1); 
sb(2) = subplot(1,2,2); 
sb(1).Position(1) = sb(1).Position(1) + 0.01;
sb(2).Position(1) = sb(2).Position(1) - 0.01;
pos{1} = sb(1).Position;
pos{2} = sb(2).Position;
close all

% Get diag fields
diag_fields = fields(diag);
for i = 1:length(diag_fields)
	if size(diag.(diag_fields{i}),2) == 1
		bad(i) = 1;
	else
		bad(i) = 0;
	end
end
diag_fields(bad==1) = [];

% Get variable to plot
for i = 1:length(diag_fields)

    % Get max bacterial biomass (for limits)
    tmp   = squeeze(diag.(diag_fields{i}));
    lims  = [0 max(max(tmp(:)),1e-24)];
    dl    = diag.(diag_fields{i}) - diag.(diag_fields{i})(1,:);
    dlims = [-max(abs(dl(:))) max(abs(dl(:)))];
	if max(dlims)==0
		continue
	end

    % Initialize figure
    fig = piofigs('sfig',0.5);
    sb(1) = subplot(1,2,1); 
    sb(2) = subplot(1,2,2); 

    % Get data up until tstep 
    set(fig,'CurrentAxes',sb(1));
    for t = 1:nrec
        plot(diag.(diag_fields{i})(t,:),-grid.z_r,'-','linewidth',0.5,'color',cmap(t,:));
        hold on
    end
    plot(diag.(diag_fields{i})(1,:),-grid.z_r,'--','linewidth',0.5,'color',rgb('Black'));
    plot(diag.(diag_fields{i})(end,:),-grid.z_r,'-','linewidth',0.5,'color',rgb('Black'));
    xlim([lims]);
    ylim([0 max(abs(grid.z_r))]);
    set(gca,'YDir','Reverse');
    ylabel('Depth ($m$)','Interpreter','Latex');
    set(gca,'FontSize',6);
    grid on
    title(upper(diag_fields{i}));
    sb(1).Position = pos{1};

    % Get difference
    set(fig,'CurrentAxes',sb(2));
    for t = 1:nrec
        plot(diag.(diag_fields{i})(t,:)-diag.(diag_fields{i})(1,:),-grid.z_r,'-','linewidth',0.5,'color',cmap(t,:));
        hold on
    end
    xlim([dlims]);
    ylim([0 max(abs(grid.z_r))]);
    set(gca,'YDir','Reverse');
    set(gca,'YTickLabel',[]);
    set(gca,'FontSize',6);
    grid on
    title('Difference from t=0');

    % Get colorbar
    cbar = colorbar('location','eastoutside');
    colormap(cmap);
    cbar.Ticks  = [cbar.Ticks(1) cbar.Ticks(end)];
    cbar.TickLabels = [0 inputs.hist_time_vec(end)/(86400*365)];
    sb(2).Position = pos{2};
    cbar.Position(3) = cbar.Position(3)*0.5;
    yl = ylabel(cbar,'Simulation years');
    yl.Position(1) = yl.Position(1) - 1;
        
    export_fig('-png',[outdir,diag_fields{i}],'-m5');
    close all
end
