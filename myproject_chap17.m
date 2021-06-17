%Hadi KM, khordad 1400
% project chapter 17

clc;clear all;close all;
load('Chap17_Data.mat')

%making raster plots and PETH for all neurons [from instruction_cue 1 sec]
%and [-500ms from movement intil +500]

%% finding active neurons
total_spikes = zeros(length(unit),1);
for neuron_num=1:length(unit)
    total_spikes(neuron_num)= length(unit(neuron_num).times);
end
activity = sort(total_spikes,'descend');
disp('most active neurons:');disp(find(total_spikes>activity(5)))
%%
neuron_num = 129 %the desired neuron . 66 and next is MI (motor cortex)

num_bins = length(0.001:0.005:1); %200
num_trials = length(instruction)
raster_ins = zeros(num_trials,num_bins); %instruction raster init
raster_go = zeros(num_trials,num_bins); % same for go
PETH_bins = 20;
peth_ins = zeros(PETH_bins,1); %Initialize the PETH with zeros
peth_go = zeros(PETH_bins,1); %Initialize the PETH with zeros

for trial=1:num_trials
    raster_ins(trial,:) = histcounts(unit(neuron_num).times, instruction(trial):0.005:instruction(trial)+1); %raster
    raster_go(trial,:)  = histcounts(unit(neuron_num).times, go(trial)-0.5:0.005:go(trial)+0.5); %raster
    peth_ins = peth_ins + histcounts(unit(neuron_num).times,instruction(trial):1/PETH_bins:instruction(trial)+1);
    peth_go =  peth_go  + histcounts(unit(neuron_num).times,go(trial)-0.5:1/PETH_bins:go(trial)+0.5);
end

%% plotting
figure('Name',['neuron:' num2str(neuron_num)],'Position',[100 100 800 500]) %Create figure for plotting
subplot(2,2,1)
imagesc(~raster_ins ) %'B' inverts 0s and 1s
colormap('gray')
title('Raster instruction')
xlabel('from instruction time to +1 sec')
ylabel('Trial')

subplot(2,2,2)
imagesc(~raster_go ) %'B' inverts 0s and 1s
colormap('gray')
title('Raster go')
xlabel('from go time to +1 sec')
ylabel('Trial')

subplot(2,2,3)
bar(0.001:1/PETH_bins:1,peth_ins(1,:)); %Plot PETH as a bar graph
title('PETH instruction')
xlabel('Time (sec)') %Label x-axis
ylabel('# of spikes') %Label y-axis

subplot(2,2,4)
bar(0.001:1/PETH_bins:1,peth_go(1,:)); %Plot PETH as a bar graph
title('PETH go')
xlabel('Time (sec)') %Label x-axis
ylabel('# of spikes') %Label y-axis

%% plotting directions for the neuron
figure('Name',['neuron:' num2str(neuron_num) 'with directions 1 to 8'],'Position',[100 0 1500 1000]) %Create figure for plotting
pett_bins = 30;
tuning_curve_ins = zeros(8,1);
tuning_curve_go = zeros(8,1);

for direction_iter = 1:8
    directional_trials = find(direction==direction_iter);
    rass_ins=zeros(length(directional_trials),num_bins);
    rass_go=zeros(length(directional_trials),num_bins);
    pett_ins = zeros(pett_bins,1); %Initialize the PETH with zeros
    pett_go = zeros(pett_bins,1); %Initialize the PETH with zeros
    
    for trial_iter = 1:length(directional_trials)
        
        rass_ins(trial_iter,:) = histcounts(unit(neuron_num).times, instruction(directional_trials(trial_iter)):0.005:instruction(directional_trials(trial_iter))+1); %raster
        rass_go(trial_iter,:)  = histcounts(unit(neuron_num).times, go(directional_trials(trial_iter))-0.5:0.005:go(directional_trials(trial_iter))+0.5); %raster
        pett_ins = pett_ins + histcounts(unit(neuron_num).times,instruction(directional_trials(trial_iter)):1/pett_bins:instruction(directional_trials(trial_iter))+1);
        pett_go =  pett_go  + histcounts(unit(neuron_num).times,go(directional_trials(trial_iter))-0.5:1/pett_bins:go(directional_trials(trial_iter))+0.5);
    end
    
    subplot(4,8,direction_iter)
    imagesc(~rass_ins ) %'B' inverts 0s and 1s
    colormap('gray')
    title('Raster instruction')
    xlabel('from instruction time to +1 sec')
    ylabel('Trial')
    clear rass_ins
    
    subplot(4,8,direction_iter+16)
    imagesc(~rass_go ) %'B' inverts 0s and 1s
    colormap('gray')
    title('Raster go')
    xlabel('from go time to +1 sec')
    ylabel('Trial')
    
    subplot(4,8,direction_iter+8)
    bar(0.001:1/pett_bins:1,pett_ins(1,:)); %Plot PETH as a bar graph
    title('PETH instruction')
    xlabel('Time (sec)')
    ylabel('# of spikes')
    ylim([0 50])
    tuning_curve_ins(direction_iter) = sum(pett_ins(1,:),2);
    
    subplot(4,8,direction_iter+24)
    bar(0.001:1/pett_bins:1,pett_go(1,:)); %Plot PETH as a bar graph
    title('PETH go')
    xlabel('Time (sec)')
    ylabel('# of spikes')
    ylim([0 50])
    tuning_curve_go(direction_iter) = sum(pett_go(1,:),2);
    
    
end

%% tuning curve
figure('Name',['neuron:' num2str(neuron_num) '   tuning curve']) %Create figure for plotting
subplot(1,2,1)
a1=gca;
x = (1:8)'*45;
plot((1:8)*45,tuning_curve_ins)
title('tuning curve for instruction')
xlabel('degrees')
ylabel('# of spikes')

% fitting cosine curve
mystring = 'p(1)+p(2)*cos(theta-p(3))'; %Cosine function in string form
myfun = inline(mystring, 'p', 'theta' ); %Converts string to a function
p= nlinfit(x, tuning_curve_ins, myfun, [1 1 0] ); %Least squares curve fit to inline function �myfun�
hold on %Allows 2 plots of the same graph
yFit = myfun(p,x); %Calculates fitted regression line
plot(x,yFit,'r') %Plots regression

predictor=[ones(8,1) sin(x) cos(x)]; %Bundle predictor variables
p=regress(tuning_curve_ins,predictor) %Linear regression
yFit=predictor*p; %Calculate fit values
theta=atan2(p(2),p(3)); %Find preferred direction from fit weights


%plot the next tuning curve

subplot(1,2,2)
plot((1:8)*45,tuning_curve_go)
title('tuning curve for go')
xlabel('degrees')
ylabel('# of spikes')

% fitting cosine curve
mystring = 'p(1)+p(2)*cos(theta-p(3))'; %Cosine function in string form
myfun = inline(mystring, 'p', 'theta' ); %Converts string to a function
p= nlinfit(x, tuning_curve_go, myfun, [1 1 0] ); %Least squares curve fit to inline function �myfun�
hold on %Allows 2 plots of the same graph
yFit = myfun(p,x); %Calculates fitted regression line
plot(x,yFit,'r') %Plots regression

predictor=[ones(8,1) sin(x) cos(x)]; %Bundle predictor variables
p=regress(tuning_curve_go,predictor) %Linear regression
yFit=predictor*p; %Calculate fit values
theta=atan2(p(2),p(3)); %Find preferred direction from fit weights

legend('tuning curve','fit')
legend('boxoff')