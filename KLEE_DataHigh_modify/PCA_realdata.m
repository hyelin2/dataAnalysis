%% Data prep

%% Data load and scailng

%% PCA

load("spike_train.mat");
load("s1_data_raw.mat");
plot(vels)
plot(pos)
plot(acc)
st1_p = spike_train1 - mean(spike_train1);
st2_p = spike_train2 - mean(spike_train2);


%% Make Plot
figure(1);
plot(spike_train1,spike_train2,'b+');
hold on; 
plot([-1 4],[0 0],'k:');
plot([0 0],[-1 4],'k:');
title('Original data');
xlim([-1 4]); ylim([-1 4]);
axis square;


figure(2);
plot(st1_p,st2_p,'b+');
hold on; plot([-1.5 1.5],[0 0],'k:'); 
plot([0 0],[-1.5 1.5],'k:');
title('Data with the means subtracted'); 
axis square;