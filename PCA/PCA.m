
%% load data and setting
load P27.mat
units = P27.units;
neuron_spike = P27.spiketimes(units);
% spiketime1 = P27.spiketimes{1,2};
% spiketime2 = P27.spiketimes{1,3};
% spiketime = [spiketime1(1:2025); spiketime2(1:2025)]
% spiketime1_p = spiketime1{1, 1} - mean(spiketime1{1, 1});
% spiketime2_p = spiketime2{1,1} - mean(spiketime2{1,1});
% plot(spiketime1_p,spiketime1_p,'b+');

% % cell size = 0.1 counting
% 
% for m = 1:150
%     data_cell= neuron_spike{1,m};
%      for n = 1:length(bin)
%           for j = 1:37100
%         spike_count(m,j) = length(find(bin(n)<data_cell & bin(n+1)>data_cell));
%          end
%      end
% end


% %% Make test data set
% data_cell= neuron_spike{1,2};
% spike_count = length(find(0.4< data_cell & 0.5>data_cell))
% Idx = cell(150,1);
% for n = 1:150
% DataReduct = neuron_spike{1,n};
% result = DataReduct(:,1:833);
% Result = num2cell(result);
% end
% Result;


%% PCA
hold on
grid on

% trial_PCAdata11 = [Spike_Count99;Spike_Count98;Spike_Count97;Spike_Count96] 
% plot(Spike_Count99,Spike_Count98)

v = nonzeros(Trial1Neuron)
v = sort(v)

hold off
[coeff11, mu11] = pca(Trial1Neuron);
plot(v)


[coeff17,mu17] = pca(trial_PCAdata17);
plot(coeff17)


[coeff,mu] = pca(trial_PCAdata18);
plot(coeff18)

% plot3(coeff11, coeff18, coeff17)
