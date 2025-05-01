% List fname
clear all
warning off
fdir  = ['../output/'];
fname = ['MicrOMZ_paper.mat'];
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

% Get variable to plot
for i = 1:length(opt.tracers)
    this_var = opt.tracers{i};

    % Get max bacterial biomass (for limits)
    tmp   = squeeze(bgc.(this_var));
    lims  = [0 max(max(tmp(:)),1e-24)];
    dl    = bgc.(this_var) - bgc.(this_var)(1,:);
    dlims = [-max(abs(dl(:))) max(abs(dl(:)))];

    % Initialize figure
    fig = piofigs('sfig',0.5);
    sb(1) = subplot(1,2,1); 
    sb(2) = subplot(1,2,2); 

    % Get data up until tstep 
    set(fig,'CurrentAxes',sb(1));
    for t = 1:nrec
        plot(bgc.(this_var)(t,:),-grid.z_r,'-','linewidth',0.5,'color',cmap(t,:));
        hold on
    end
    plot(bgc.(this_var)(1,:),-grid.z_r,'--','linewidth',0.5,'color',rgb('Black'));
    plot(bgc.(this_var)(end,:),-grid.z_r,'-','linewidth',0.5,'color',rgb('Black'));
    xlim([lims]);
    ylim([0 max(abs(grid.z_r))]);
    set(gca,'YDir','Reverse');
    ylabel('Depth ($m$)','Interpreter','Latex');
    set(gca,'FontSize',6);
    grid on
    title(upper(this_var));
    sb(1).Position = pos{1};

    % Get difference
    set(fig,'CurrentAxes',sb(2));
    for t = 1:nrec
        plot(bgc.(this_var)(t,:)-bgc.(this_var)(1,:),-grid.z_r,'-','linewidth',0.5,'color',cmap(t,:));
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
        
    export_fig('-png',[outdir,this_var],'-m5');
    close all
end
