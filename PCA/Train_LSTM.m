%% Data load and plot
load P27.mat
units = P27.units;
neuron_spike = P27.spiketimes(units);
Neuron_spike = tall(neuron_spike);
data = spiketrain_trial;
data = [data{:}];

figure
plot(data)
xlabel("spike")
ylabel("time")
title("spiketime")

%% Train and testing(LSTM)

numTimeStepsTrain = floor(0.6*numel(data));

dataTrain = data(1:numTimeStepsTrain+1);
dataTest = data(numTimeStepsTrain+1:end);

%% 데이터 표준화 / fitting을 위해 평균 0, 분산 1이 되도록 데이터 표준화
mu = mean(dataTrain);
sig = std(dataTrain);

dataTrainStandardized = (dataTrain - mu) / sig;

%% 예측변수와 응답 변수준비
XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);

%% LSTM 신경망 아키텍쳐 정의

 numFeatures = 1;
 numResponses = 1;
 numHiddenUnits = 200; % 200개의 hidden 유닛 설정
 miniBatchSize = 64;
 layers = [...
     sequenceInputLayer(numFeatures)
     lstmLayer(numHiddenUnits)
     fullyConnectedLayer(numResponses)
     regressionLayer];

 %% 훈련 옵션 지정
 options = trainingOptions('adam', ...
     'MaxEpochs',500,...
     "GradientThreshold",1,...
     'InitialLearnRate',0.005,...
     'LearnRateSchedule','piecewise',...
     'LearnRateDropPeriod',125,...
     'LearnRateDropFactor',0.2,...
     'Verbose',0,...
     'Plots','training-progress');

 %% trainig network 지정, 신경망 훈련
 net = trainNetwork(XTrain, YTrain, layers, options);

 %% Predict update state
 dataTestStandardized = (dataTest - mu) / sig;
 XTest = dataTestStandardized(1:end-1); %test data 표준화

 net = predictAndUpdateState(net,XTrain);
 [net,YPred] = predictAndUpdateState(net,YTrain(end));

 numTimeStepsTest = numel(XTest);

 for i = 2:numTimeStepsTest
     [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','auto');
 end

YPred = sig*YPred + mu;
% 예측 표준화 해제, 훈련 진행상황은 RMSE(제곱평균제곱근)를 보고함.

YTest = dataTest(2:end);
rmse = sqrt(mean((YPred-YTest).^2))


%% 전망값을 사용하여 훈련 시계열 플로팅
figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[data(numTimeStepsTrain) YPred], '.-')
hold off

xlabel("Spike")
ylabel("Time")
title("Forecast")
legend(["Observed" "Forecast"])

%% 전망 값과 테스트 데이터 비교
figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred, '.-')
hold off
legend(["Observed" "Forecast"])
ylabel("Time")
title("Forecast")

subplot(2,1,2)
stem(YPred - YTest)
xlabel("Spike")
ylabel("Error")
title("RMSE = " + rmse)

%% 관측값으로 Neural network update
net = resetState(net);
net = predictAndUpdateState(net, XTrain);

YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net, XTest(:,i),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;
% accuracy, 잔차 계산(평균제곱근오차)
rmse = sqrt(mean((YPred-YTest).^2))
% mae = mae()
%% test data와 비교
figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred, '.-')
hold off
legend(["Observed" "Predicted"])
ylabel("Time")
title("Forecast with Update")

subplot(2,1,2)
stem(YPred - YTest)
xlabel("Spike")
ylabel("Error")
title("RMSE = " + rmse)



