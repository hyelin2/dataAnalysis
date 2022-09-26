%% neuron-time(-1~1) plot 그리는 함수
load P27.mat;
units = P27.units;
neuron_spike = P27.spiketimes(units);
eventtime = P27.eventtimes;
for n = 1:150
%% [eventtime(i)-1 < neuraldata < eventtime(i)+1]   
% Make Trial(i)
for i = 1:110
% for n = 1 : numel(neuron_spike)
% spiketime = neuron_spike{1,n};
% eval(["Idn = find(spiketime > eventtime(" + num2str(i) + ")-1 & spiketime < eventtime(" + num2str(i) + "));" ]);
spiketime = neuron_spike{1,n};
Idx = find(spiketime > eventtime(i)-1 & spiketime < eventtime(i));
trialA = spiketime(Idx);
% eval(["Idx = find(spiketime > eventtime(" + num2str(i) + " & spiketime < eventtime(" + num2str(i) + ")+1);" ]);
Idx = find(spiketime > eventtime(i) & spiketime < eventtime(i)+1);
trialB = spiketime(Idx);
% end
eval(["trial" + num2str(i) + " = [trialA trialB];"]);
% eval(["trial" + num2str(i) + " = rescale(trial" + num2str(i) + ",-1,1);"]);
% trial110 = [trialA trialB];
% trial110 = rescale(trial110,-1,1);
end

%% 

% trial vector 보간
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
    Count_sum{n} = eval(["Spike_Count" + num2str(n) + "; "]);
    eval(["clear Spike_Count" + num2str(n) + ";"]);
end


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


% 
% for i = 1:20
% binNum = [eventtime(n)-1:bin:eventtime(n)+1];    
% Spike_Count1(1,i) = numel(find(binNum(i) < trial_sum{1,1}(1,:) & binNum(i+1) > trial_sum{1,1}(1,:)));
% Spike_Count1(2,i) = numel(find(binNum(i) < trial_sum{1,1}(2,:) & binNum(i+1) > trial_sum{1,1}(2,:)));
% end


% % Make Counting Vector
% Csum1 = horzcat(Spike_Count1, Spike_Count2);
% for i = 2:109
%  eval(["Csum" + num2str(i) + " = horzcat(Csum" + num2str(i-1) + ",Spike_Count" + num2str(i+1) + ");"]);
% end
% 
% 
% 
% for i = 1:108
%     eval(["clear Csum" + num2str(i) + ";"]);
% end
% Csum109 = reshape(Csum109,[110,2000]);


% bin = zeros(21,150);
%     for l = 1:21
%         for m = 1:numel(Trial1Neuron)-1
%         if binNum(l) <= Trial1Neuron(m) & binNum(l+1) >= Trial1Neuron(m)
%           a = numel(Trial1Neuron(m));
%           bin(l,1) = a;
%         else
%           bin(l,1) = 0;
%           m=m+1;
%         end
%         end
%     
%     end

% for i = 1:110
% % eval(["idx = find(Trial" + num2str(i) + "Neuron);"]);
% eval(["Trial" + num2str(i) + "Neuron(idx) = 1;"]);
% end




% 
% Spike_sum1 = horzcat(Spike_Count1, Spike_Count2);
% for i = 2:109
%  eval(["Spike_sum" + num2str(i) + " = horzcat(Spike_sum" + num2str(i-1) + ",Spike_Count" + num2str(i+1) + ");"]);
% end
% Spike_Sum = reshape(Spike_sum109,[110,2000]);

% % Make PCA each trial
for n = 1:96
Count_sum{n} = eval(["Spike_Count" + num2str(n) + ";"]);
end

Count_sum{1,1} = Spike_Count9
B = struct('data', Count_sum ,'condition','reach1','epochstars',1,'epochcolors',[0,0,1]);
B = B';
DataHigh(B,'DimReduce')


% 
% % Make PCA all trial
% % Make PCA
% B = struct('data', trial_sum ,'condition','reach1','epochstars',1,'epochcolors',[0,0,1]);
% B = B';
% DataHigh(B,'DimReduce')

