%% Data load and setting
load P27.mat;
units = P27.units;
neuron_spike = P27.spiketimes(units);
eventtime = P27.eventtimes;


%% [eventtime(i)-1 < neuraldata < eventtime(i)+1]   
% Make Trial(i)
for n = 1:150
    for i = 1:110
        spiketime = neuron_spike{1,n};
        Idx = find(spiketime > eventtime(i)-1 & spiketime < eventtime(i));
        trialA = spiketime(Idx);
        Idx = find(spiketime > eventtime(i) & spiketime < eventtime(i)+1);
        trialB = spiketime(Idx);
        eval(["trial" + num2str(i) + " = [trialA trialB];"]);
    end

%% Shape Trial vector
    for i = 1:110
eval(["trial"+ num2str(i) + "(1,500) = 0;"]);
end
% Make Sum Vector
    sum1 = vertcat(trial1, trial2);
    for i = 2:109
    eval(["sum" + num2str(i) + " = vertcat(sum" + num2str(i-1) + ",trial" + num2str(i+1) + ");"]);
end

%clear sum
    for i = 1:108
    eval(["clear sum" + num2str(i) + ";"]);
end

%clear trial
    for i = 1:110
    eval(["clear trial" + num2str(i) + ";"]);
end
        eval(["trialTotime" + num2str(n) + "= sum109;"]);
        clear sum109;
end

%% Make Time to Neuron vector
% Make Trial Vector
Tsum1 = horzcat(trialTotime1, trialTotime2);
for i = 2:149
 eval(["Tsum" + num2str(i) + " = horzcat(Tsum" + num2str(i-1) + ",trialTotime" + num2str(i+1) + ");"]);
end

%clear Tsum
for i = 1:148
    eval(["clear Tsum" + num2str(i) + ";"]);
end
% clear trialTotime
for i = 1:150
    eval(["clear trialTotime" + num2str(i) + ";"]);
end

% 150 Neuron trial1 spike train 
for n = 1:110
    eval(["Trial" + num2str(n) + "Neuron = reshape(Tsum149(" + num2str(n) + ",:),[150,500]);" ]);
end
clear Tsum149;

% Sum of Trial(cell vector 생성)
n = 110;
trial_sum = cell(n,1);
for n = 1:n
    trial_sum{n} = eval(["Trial" + num2str(n) + "Neuron ; "]);
    eval(["clear Trial" + num2str(n) + "Neuron;"]);
end

n = 110;
Count_sum = cell(n,1);
for n = 1:n

%% Counting time to neuron(bin = 0.1)

bin = 0.01;
for n = 1:110
    for j = 1:150
        for i = 1:200
            binNum = [eventtime(n)-1:bin:eventtime(n)+1];
            eval(["Spike_Count" + num2str(n) + "(j,i) = numel(find(binNum(i) < trial_sum{n,1}(j,:) & binNum(i+1) > trial_sum{n,1}(j,:)));"]);
        end
    end
end
    Count_sum{n} = eval(["Spike_Count" + num2str(n) + "; "]);
    eval(["clear Spike_Count" + num2str(n) + ";"]);
end


%% Make PCA each trial
% for n = 1:109
% Count_sum{n} = eval(["Spike_Count" + num2str(n) + ";"]);
% end

Count_sum{1,1} = spiketrain_trial{1,1}
% Count_sum{1,1} = Spike_Count9
B = struct('data', spiketrain_trial ,'condition','reach1','epochstars',1,'epochcolors',[0,0,1]);
B = B';
DataHigh(B,'DimReduce')

% % Make PCA all trial
% % Make PCA
% B = struct('data', trial_sum ,'condition','reach1','epochstars',1,'epochcolors',[0,0,1]);
% B = B';
% DataHigh(B,'DimReduce')

