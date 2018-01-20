%% FIG4FIELDSPECTRA
% This short code checks the spectra generated by the field model.
% Version: 2013jan15 by Cliff Kerr (cliffk@neurosim.downstate.edu)

% Load data
disp('Loading...')
base='../output/output'; % Base file name -- needs to match runsim
types={'tc','bg','pa'};
fullnames={'TC','BG','PD'};
nfiles=length(types);
traces=cell(nfiles,1);
ntraces=zeros(nfiles,1);
file=cell(nfiles,1);
reorder=[1 2 8 9 3:7]; % The populations are listed in a strange order in the file
for f=1:nfiles
    file{f}=sprintf('%s-%s-fld.txt',base,types{f});
    data=load(file{f});
    ntraces(f)=size(data,2)-1;
    traces{f}=cell(ntraces(f),1);
    if f==nfiles, colors=vectocolor(1:ntraces(f))/1.1; end % Make green and yellow a bit darker with the 1.1
    for t=1:ntraces(f)
        traces{f}{t}=data(:,t+1)*1000; % Convert from s to ms
        traces{f}{t}=traces{f}{t}-mean(traces{f}{t});
        traces{f}{t}=traces{f}{t}/mean(abs(traces{f}{t}));
    end
    if f==1
        traces{f}=traces{f}([1 2 4 3]); % Reorder E, I, S, R
    elseif f==2 || f==3 % Reorder so cortical populations are first
        traces{f}=traces{f}(reorder);
    end
end
time=data(:,1)/1000; % Convert from ms to s
dt=diff(time(1:2)); % Assume time is uniformly spaced
names={{'   Excitatory','   Inhibitory','   Relay','   Reticular'},{'Excitatory','Inhibitory','Striatal D1','Striatal D2','GPi','GPe','Subthalamic','Relay','Reticular'}};
names{2}=names{2}(reorder);
names{3}=names{2};



% Calculate spectra
disp('Calculating spectra...')
spectra=cell(nfiles,1);
for f=1:nfiles
    spectra{f}=cell(ntraces(f),1);
    for t=1:ntraces(f)
        [spectra{f}{t}, F]=spectrum('pl',traces{f}{t},round(1/dt),1);
        spectra{f}{t}(1)=spectra{f}{t}(2); % To avoid problems with smoothing
        spectra{f}{t}=cksmooth(spectra{f}{t},200);
    end
end


%% Plotting
disp('Plotting...')
figure('position',[200 200 800 700])
maxtime=0.5;
for f=1:nfiles
    cksubplot([2 nfiles],[1 f],[90 80],[3 5],[5 -1],[1 1]); hold on
    for t=1:ntraces(f)
        if t<=4, width=2; else width=1; end % This code selects out E, I, S, R neurons to make stand out
        plot(time,traces{f}{t}-(1+(f==1))*4*t,'color',colors(t,:),'linewidth',width)
    end
    legend(names{f},'location','eastoutside')
    if f==nfiles, xlabel('Time (s)'), else set(gca,'xticklabel',[]), end
    set(gca,'xtick',0:0.1:maxtime)
    set(gca,'yticklabel',[])
    ylabel(sprintf('Normalized\nvoltage'))
    ylim([-40 0])
    title(fullnames{f})
    xlim([0 maxtime])
    
    cksubplot([2 nfiles],[2 f],[90 80],[5 5],[8 -1],[1 1]); hold on
    for t=ntraces(f):-1:1
        if t<=4, width=2; else width=1; end
        plot(F,spectra{f}{t},'color',colors(t,:),'linewidth',width)
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-3 1e2])
    xlim([1 100])
    if f==nfiles, xlabel('Frequency (Hz)'), else set(gca,'xticklabel',[]), end
    ylabel('Normalized power')
    title(fullnames{f})
    set(gca,'ticklength',0.02*[1 1])
end


%% Plot labels
ax=axes('position',[0 0 1 1],'units','normalized'); hold on
text(0.01,0.97,'A','fontsize',16,'fontweight','bold')
text(0.49,0.97,'B','fontsize',16,'fontweight','bold')
annotation('arrow',[0.02 0.045],[0.45 0.45],'color',[1 0.3 0])
annotation('arrow',[0.02 0.045],[0.43 0.43],'color',[1 0.3 0])
annotation('arrow',[0.02 0.045],[0.135 0.135],'color',[0 0.3 1])
annotation('arrow',[0.02 0.045],[0.08 0.08],'color',[0 0.3 1])
annotation('arrow',[0.81 0.81],[0.46 0.49],'color',[0 0.6 0])
annotation('arrow',[0.81 0.81],[0.14 0.17],'color',[0 0.6 0])
xlim([0 1])
ylim([0 1])
set(ax,'visible','off')
