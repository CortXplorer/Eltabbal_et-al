%% for Hyase treated data 
clear all ; clc; 
load('DataforGroupanalysisHYASE.mat'); 
for i = 1 : size(UBSSA_HBSSA{1},2)
    BF_list(i) = UBSSA_HBSSA{1}{1,i}.BF_table(2);
end
for anim = 1 :size(UBSSA_HBSSA{1},2)
    for cond = 1 :size(UBSSA_HBSSA,2)
        data = mean(UBSSA_HBSSA{cond}{anim}.SingleTrialCSD{1, BF_list(anim)},3); 
        csdbeforeafter(:,cond) = mean(data(:,200:250),2); 
    end
    [crosscorrvalues(anim,:),~,~] = crosscorr(csdbeforeafter(:,1),csdbeforeafter(:,2)); 
    clear csdbeforeafter
    x = -20:20; 
    plot(x, crosscorrvalues(anim,:)); hold on 
    [~, ind ]= max(crosscorrvalues(anim,:));
    maxcorrchannel(anim)= abs(ind-20); 
end
title('Channels cross-correlation before and after Hyase')
xlabel('Channels')
ylabel('correlation')
maximum_shiftprepost = mean(maxcorrchannel); 
std_shiftprepost = std(maxcorrchannel); 
        
fprintf('The mean of electrode shift is %f channel with standard deviation of %f',maximum_shiftprepost, std_shiftprepost)        
%% for control treated data 
clear all ; clc; 
load('ControlECMDATA.mat')
UB_controlA{1}{1}= FreqTunDataContainer{1}.ECMc001_0007; 
UB_controlA{1}{2}= FreqTunDataContainer{1}.ECMc002_0008;
UB_controlA{1}{3}= FreqTunDataContainer{1}.ECMc003_0010; 

UB_controlA{2}{1}= FreqTunDataContainer{2}.ECMc001_0026; 
UB_controlA{2}{2}= FreqTunDataContainer{2}.ECMc002_0014; 
UB_controlA{2}{3}= FreqTunDataContainer{2}.ECMc003_0023; 
% 

for i = 1 : size(UB_controlA{1},2)
    BF_list(i) = UB_controlA{1}{1,i}.BF_table(2);
end
for anim = 1 :size(UB_controlA{1},2)
    for cond = 1 :size(UB_controlA,2)
        data = mean(UB_controlA{cond}{anim}.SingleTrialCSD{1, BF_list(anim)},3); 
        csdbeforeafter(:,cond) = mean(data(:,200:250),2); 
    end
    [crosscorrvalues(anim,:),~,~] = crosscorr(csdbeforeafter(:,1),csdbeforeafter(:,2)); 
    clear csdbeforeafter
    x = -20:20; 
    plot(x, crosscorrvalues(anim,:)); hold on 
    [~, ind ]= max(crosscorrvalues(anim,:));
    maxcorrchannel(anim)= abs(ind-20); 
end
title('Channels cross-correlation before and after Hyase')
xlabel('Channels')
ylabel('correlation')
maximum_shiftprepost = mean(maxcorrchannel); 
std_shiftprepost = std(maxcorrchannel); 
        
fprintf('The mean of electrode shift is %f channel with standard deviation of %f',maximum_shiftprepost, std_shiftprepost) 




