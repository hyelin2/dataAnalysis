spiketimes=P27.spiketimes;
units=P27.units;
eventtimes=P27.eventtimes;
binduration=0.1;
tmin=-1;
tmax=1;


% if size(spiketimes,1)<size(spiketimes,2)
%   spiketimes = spiketimes';
% end

ntrials=length(eventtimes); nunits=length(units);

timebins=[tmin:binduration:(tmax)];

spikespertrial=zeros(nunits, length(timebins));

spiketrain_trial=[];

for trial=1:ntrials;
    
    trialtime=eventtimes(trial); 

    for unit=1:length(units)

    unitspiketimes=spiketimes{units(unit)};
    
    spikeinds=find(unitspiketimes<(trialtime+max(timebins)) & unitspiketimes>(trialtime+min(timebins)));

    relative_spiketimes=unitspiketimes(spikeinds)-trialtime;

    [n,bins]=histc(relative_spiketimes, timebins);	

    spikespertrial(unit,:)= n;

    end
    
    spiketrain_trial{trial}=spikespertrial;

    spikespertrial=zeros(nunits, length(timebins));

end