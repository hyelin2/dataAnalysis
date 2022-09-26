function [rasterx, rastery] = plotraster(stimulustimes, spiketimes, tmin, tmax)

numberoftrials=length(stimulustimes);

close all

figure(1)

tbins=[tmin:0.0001:tmax];

rasterx=[];
rastery=[];

for trial=1:numberoftrials;
    
    trialtime=stimulustimes(trial); 
    
    spikeinds=find(spiketimes<(trialtime+tmax) & spiketimes>(trialtime+tmin));
    
    relative_spiketimes=spiketimes(spikeinds)-trialtime;

    [n,bins]=histc(relative_spiketimes, tbins);	

    rasterx=[rasterx tbins(bins)];
    rastery=[rastery trial*ones(1,sum(n))];   
         
end

%plot([rasterx; rasterx],[rastery-0.4 rastery+0.4],'Color', 'k')  %raster plot with vertical lines (slow on Octave).

scatter(rasterx, rastery, 2, 'k', 'o')

h = get(gcf, 'currentaxes');
set(h, 'fontsize', 16, 'linewidth', 0.5);

xlabel('time (s)')
ylabel('trial number')

axis([tmin tmax 0 numberoftrials+0.5])

