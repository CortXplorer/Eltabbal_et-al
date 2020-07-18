function Tabbal_et_al_CSD_group_analysis
%this function will group different variables for each condition and make
%comparative group statistics (Author: EL-Tabbal and Happel July 2020)

path = 'Add where you want to save all figures add your "Dataforgroupanalysishyase.mat in the folder containing the figure folder';
cd(path)
cd ..  
load ('DataforGroupanalysisHYASE.mat');
UBSSA= UBSSA_HBSSA{1,1};
HBSSA= UBSSA_HBSSA{1,2};
for i = 1 : size(UBSSA_HBSSA{1},2)
    BF_list(i) = UBSSA_HBSSA{1}{1,i}.BF_table(2); % BF based on the granular Avgsink in the table 10Feb2020)
end
layers = {'late_infragranular' , 'early_infragranular' , 'granular' , 'supragraular'};
clear i 
%% the effect of hyase on CSD parameters
% 1st peak amplitude ( figure not included in the paper but good to see the
% effects and compare ) 
amplitude_peak_before_hyase_matrix = nan(size(UBSSA,2), size(UBSSA{1}.AvgSink_Peak,2), size(UBSSA{1}.AvgSink_Peak,3)) ;
amplitude_peak_after_hyase_matrix = nan(size(HBSSA,2), size(HBSSA{1}.AvgSink_Peak,2), size(HBSSA{1}.AvgSink_Peak,3)) ;
peakBAmatrix = {amplitude_peak_before_hyase_matrix,amplitude_peak_after_hyase_matrix};

for i = 1:size(peakBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            peakBAmatrix{1,i}(anim,:,k) = cell2mat(UBSSA_HBSSA{1,i}{1,anim}.AvgSink_Peak(2,:,k));
        end
    end
end
clear i anim k
BF_nonBF_before_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_after_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_peakBAmatrix = {BF_nonBF_before_hyse_matrix,BF_nonBF_after_hyse_matrix};
for i = 1:size(BF_nonBF_peakBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if BF_list(anim) == 1
                BF_nonBF_peakBAmatrix{1,i}(anim,[3:5],k) = peakBAmatrix{1,i}(anim,[BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) == 2
                BF_nonBF_peakBAmatrix{1,i}(anim,[2:5],k) = peakBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) > 2 && BF_list(anim) < 7
                BF_nonBF_peakBAmatrix{1,i}(anim,:,k) = peakBAmatrix{1,i}(anim,[BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) == 7
                BF_nonBF_peakBAmatrix{1,i}(anim,[2:5],k) = peakBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)-2],k);
            else
                BF_nonBF_peakBAmatrix{1,i}(anim,[2,3,5],k) = peakBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)-2],k);
            end
        end
    end
end
clear i anim k amplitude_peak_after_hyase_matrix amplitude_peak_before_hyase_matrix peakBAmatrix BF_nonBF_after_hyse_matrix BF_nonBF_before_hyse_matrix
% ploting peak amplitude tuning before and after hyase in different layers
colors = {'b','r'} ;
FIGURENAME = 'PEAK AMPLITUDE TUNING CURVES BEFORE AND AFTER INJECTION (July 2020) '
h_fig = figure
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
for k = 1:size(BF_nonBF_peakBAmatrix{1,1},3);
    for i = 1 : size(BF_nonBF_peakBAmatrix,2);
        a= nanmean(BF_nonBF_peakBAmatrix{1,i}(:,:,k));
        b= nanstd(BF_nonBF_peakBAmatrix{1,i}(:,:,k));
        c = b./sqrt(9);
        subplot(2,2,k);
        title(layers{k});
        errorbar(a,c,colors{i});
        xticks([1,2,3,4,5,6,7]);
        set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
        xlim([0,6]);
        ylabel('Peak amplitude');
        legend('Untreated','Hyase' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
cd(path)
saveas(h_fig,FIGURENAME,'fig')
close all
clear a b c colors i k  FIGURENAME

%% 2nd Peak Latency ( figure not included in the paper but good to see the
% effects and compare )
amplitude_peaklatency_before_hyase_matrix = nan(size(UBSSA,2), size(UBSSA{1}.AvgSink_Peak,2), size(UBSSA{1}.AvgSink_Peak,3)) ;
amplitude_peaklatency_after_hyase_matrix = nan(size(HBSSA,2), size(HBSSA{1}.AvgSink_Peak,2), size(HBSSA{1}.AvgSink_Peak,3)) ;
peaklatencyBAmatrix = {amplitude_peaklatency_before_hyase_matrix,amplitude_peaklatency_after_hyase_matrix};
UBSSA_HBSSA = {UBSSA, HBSSA};
for i = 1:size(peaklatencyBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            peaklatencyBAmatrix{1,i}(anim,:,k) = cell2mat(UBSSA_HBSSA{1,i}{1,anim}.AvgSink_Peak(1,:,k));
        end
    end
end
clear i anim k
BF_nonBF_before_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_after_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_peaklatencyBAmatrix = {BF_nonBF_before_hyse_matrix,BF_nonBF_after_hyse_matrix};
for i = 1:size(BF_nonBF_peaklatencyBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if BF_list(anim) == 1
                BF_nonBF_peaklatencyBAmatrix{1,i}(anim,[3:5],k) = peaklatencyBAmatrix{1,i}(anim,[BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) == 2
                BF_nonBF_peaklatencyBAmatrix{1,i}(anim,[2:5],k) = peaklatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) > 2 && BF_list(anim) < 7
                BF_nonBF_peaklatencyBAmatrix{1,i}(anim,:,k) = peaklatencyBAmatrix{1,i}(anim,[BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim)== 7
                BF_nonBF_peaklatencyBAmatrix{1,i}(anim,[2:5],k) = peaklatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)-2],k);
            else
                BF_nonBF_peaklatencyBAmatrix{1,i}(anim,[2,3,5],k) = peaklatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)-2],k);
                
            end
        end
    end
end
clear i anim k amplitude_peaklatency_after_hyase_matrix amplitude_peaklatency_before_hyase_matrix peaklatencyBAmatrix BF_nonBF_after_hyse_matrix BF_nonBF_before_hyse_matrix
% getting rid of all the -1 (values that crossing the cutoffs are replaced by -1)
for i = 1 : size(BF_nonBF_peaklatencyBAmatrix,2)
    a= find (BF_nonBF_peaklatencyBAmatrix{1,i} == -1)
    BF_nonBF_peaklatencyBAmatrix{1,i}(a)= nan
end
clear a i
% ploting peak amplitude tuning before and after hyase in different layers
colors = {'b','r'} ;
FIGURENAME = 'PEAK Latency TUNING CURVES BEFORE AND AFTER INJECTION (july 2020)'
h_fig = figure
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
for k = 1:size(BF_nonBF_peaklatencyBAmatrix{1,1},3);
    for i = 1 : size(BF_nonBF_peaklatencyBAmatrix,2);
        a= nanmean(BF_nonBF_peaklatencyBAmatrix{1,i}(:,:,k));
        b= nanstd(BF_nonBF_peaklatencyBAmatrix{1,i}(:,:,k));
        c = b./sqrt(size(UBSSA,2));
        subplot(2,2,k);
        title(layers{k});
        errorbar(a,c,colors{i});
        xticks([1,2,3,4,5,6,7]);
        set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
        xlim([0,6]);
        ylabel('Peak latency(s)');
        legend('Untreated','Hyase' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
cd(path)
saveas(h_fig,FIGURENAME,'fig')
close all
clear a b c colors i k  FIGURENAME
%% 3rd Onset latency  ( figure not included in the paper but good to see the
% effects and compare )
amplitude_onsetlatency_before_hyase_matrix = nan(size(UBSSA,2), size(UBSSA{1}.AvgSink_Peak,2), size(UBSSA{1}.AvgSink_Peak,3)) ;
amplitude_onsetlatency_after_hyase_matrix = nan(size(HBSSA,2), size(HBSSA{1}.AvgSink_Peak,2), size(HBSSA{1}.AvgSink_Peak,3)) ;
onsetlatencyBAmatrix = {amplitude_onsetlatency_before_hyase_matrix,amplitude_onsetlatency_after_hyase_matrix};
UBSSA_HBSSA = {UBSSA, HBSSA};
for i = 1:size(onsetlatencyBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            onsetlatencyBAmatrix{1,i}(anim,:,k) = cell2mat(UBSSA_HBSSA{1,i}{1,anim}.Latency_Sink(2,:,k));
        end
    end
end
clear i anim k
BF_nonBF_before_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_after_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_onsetlatencyBAmatrix = {BF_nonBF_before_hyse_matrix,BF_nonBF_after_hyse_matrix};
for i = 1:size(BF_nonBF_onsetlatencyBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if BF_list(anim) == 1
                BF_nonBF_onsetlatencyBAmatrix{1,i}(anim,[3:5],k) = onsetlatencyBAmatrix{1,i}(anim,[BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) == 2
                BF_nonBF_onsetlatencyBAmatrix{1,i}(anim,[2:5],k) = onsetlatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) > 2 && BF_list(anim) < 7
                BF_nonBF_onsetlatencyBAmatrix{1,i}(anim,:,k) = onsetlatencyBAmatrix{1,i}(anim,[BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) == 7
                BF_nonBF_onsetlatencyBAmatrix{1,i}(anim,[2:5],k) = onsetlatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)-2],k);
            else
                BF_nonBF_onsetlatencyBAmatrix{1,i}(anim,[2,3,5],k) = onsetlatencyBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)-2],k);
                
            end
        end
    end
end
clear i anim k amplitude_onsetlatency_after_hyase_matrix amplitude_onsetlatency_before_hyase_matrix onsetlatencyBAmatrix BF_nonBF_after_hyse_matrix BF_nonBF_before_hyse_matrix
for i = 1 : size(BF_nonBF_onsetlatencyBAmatrix,2)
    a= find (BF_nonBF_onsetlatencyBAmatrix{1,i} == -1)
    BF_nonBF_onsetlatencyBAmatrix{1,i}(a)= nan
end
clear a i
% ploting peak amplitude tuning before and after hyase in different layers
colors = {'b','r'} ;
FIGURENAME = 'ONSET Latency TUNING CURVES BEFORE AND AFTER INJECTION (july2020)'
h_fig = figure
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
for k = 1:size(BF_nonBF_onsetlatencyBAmatrix{1,1},3);
    for i = 1 : size(BF_nonBF_onsetlatencyBAmatrix,2);
        a= nanmean(BF_nonBF_onsetlatencyBAmatrix{1,i}(:,:,k));
        b= nanstd(BF_nonBF_onsetlatencyBAmatrix{1,i}(:,:,k));
        c = b./sqrt(size(UBSSA,2));
        subplot(2,2,k);
        title(layers{k});
        errorbar(a,c,colors{i});
        xticks([1,2,3,4,5,6,7]);
        set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
        xlim([0,6]);
        ylabel('Onset latency(s)');
        legend('Untreated','Hyase' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
cd(path)
saveas(h_fig,FIGURENAME,'fig')
close all
clear a b c colors i k  FIGURENAME

%% 4th RMS_AvgSink This part of the script generate figure 2.A. in the paper 
RMS_AvgSink_before_hyase_matrix = nan(size(UBSSA,2), size(UBSSA{1}.AvgSink_Peak,2), size(UBSSA{1}.AvgSink_Peak,3)) ;
RMS_AvgSink_after_hyase_matrix = nan(size(HBSSA,2), size(HBSSA{1}.AvgSink_Peak,2), size(HBSSA{1}.AvgSink_Peak,3)) ;
RMS_AvgSinkBAmatrix = {RMS_AvgSink_before_hyase_matrix,RMS_AvgSink_after_hyase_matrix};
UBSSA_HBSSA = {UBSSA, HBSSA};
for i = 1:size(RMS_AvgSinkBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            RMS_AvgSinkBAmatrix{1,i}(anim,:,k) = UBSSA_HBSSA{1,i}{1,anim}.RMS_AvgSink(1,:,k);
        end
    end
end
clear i anim k
BF_nonBF_before_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_after_hyse_matrix= nan(size(UBSSA,2),5,size(UBSSA{1}.AvgSink_Peak,3));
BF_nonBF_RMS_AvgSinkBAmatrix = {BF_nonBF_before_hyse_matrix,BF_nonBF_after_hyse_matrix};
for i = 1:size(BF_nonBF_RMS_AvgSinkBAmatrix,2)
    for anim = 1:size(UBSSA,2) % number of animals
        for k = 1:size(UBSSA{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if BF_list(anim) == 1
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(anim,[3:5],k) = RMS_AvgSinkBAmatrix{1,i}(anim,[BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim)== 2
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(anim,[2:5],k) = RMS_AvgSinkBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim) > 2 && BF_list(anim) < 7
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(anim,:,k) = RMS_AvgSinkBAmatrix{1,i}(anim,[BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2],k);
            elseif BF_list(anim)== 7
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(anim,[2:5],k) = RMS_AvgSinkBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)-2],k);
            else
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(anim,[2,3,5],k) = RMS_AvgSinkBAmatrix{1,i}(anim,[BF_list(anim)-1,BF_list(anim),BF_list(anim)-2],k);
                
            end
        end
    end
end
clear i anim k RMS_AvgSink_after_hyase_matrix RMS_AvgSink_before_hyase_matrix RMS_AvgSinkBAmatrix BF_nonBF_after_hyse_matrix BF_nonBF_before_hyse_matrix

% ploting RMS_AvgSink_ tuning before and after hyase in different layers
colors = {'b','r'} ;
figure
FIGURENAME = 'RMS_AvgSink TUNING CURVES BEFORE AND AFTER INJECTION ( July2020) '
h_fig = figure
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])

for k = 1:size(BF_nonBF_RMS_AvgSinkBAmatrix{1,1},3);
    for i = 1 : size(BF_nonBF_RMS_AvgSinkBAmatrix,2);
        a= nanmean(BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(:,:,k));
        b= nanstd(BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(:,:,k));
        c = b./sqrt(size(UBSSA,2));
        subplot(2,2,k);
        title(layers{k});
        errorbar(a,c,colors{i});
        xticks([1,2,3,4,5,6,7]);
        set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
        ylabel('RMS-Avg-Sink amplitude');
        xlim([0,6])
        legend('Untreated','Hyase' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
cd(path)
saveas(h_fig,FIGURENAME,'fig')
close all
clear a b c colors i k  FIGURENAME
cd(path);
cd ..

save ('BF_nonBF_RMS_AvgSinkBAmatrix' ,  'BF_nonBF_RMS_AvgSinkBAmatrix')
save ('BF_nonBF_onsetlatencyBAmatrix' , 'BF_nonBF_onsetlatencyBAmatrix' )
save ('BF_nonBF_peakBAmatrix' , 'BF_nonBF_peakBAmatrix')
save ('BF_nonBF_peaklatencyBAmatrix' , 'BF_nonBF_peaklatencyBAmatrix')

%%  similarity layer index (  This part of the script generate figure 2.B. in the paper)
%  first get the layer index which is a ration of the value in specific
%  layer relative to granular layer .. then calculate layer symmetry index
%  by taking difference of both layers over the sum
variables = {'BF_nonBF_peakBAmatrix.mat','BF_nonBF_onsetlatencyBAmatrix.mat','BF_nonBF_RMS_AvgSinkBAmatrix.mat'};
for k = 1 : size (variables,2);
    matrix = struct2cell(load(variables{k}));
    for i = 1 : size(matrix{1,1},2);
        early_infra_index{i} = matrix{1,1}{1,i}(:,:,2)./matrix{1,1}{1,i}(:,:,3);
        Supra_infra_index{i}= matrix{1,1}{1,i}(:,:,4)./matrix{1,1}{1,i}(:,:,3);
        Sime_index{i} = (Supra_infra_index{i}- early_infra_index{i})./ (early_infra_index{i} + Supra_infra_index{i}) ;
    end
    similarity_index{k} = Sime_index;
end
clear i k matrix variables late_infra_index Supra_infra_index Sime_index

Peak_sime_index = similarity_index{1};
Onset_sime_index = similarity_index{2};
RMS_sime_index= similarity_index{3};
clear similarity_index
% ploting sime_index
colors = {'b' , 'r'};
figure();
for k = 1 : size(RMS_sime_index,2) ;
    a= nanmean(RMS_sime_index{1,k});
    b= nanstd(RMS_sime_index{1,k});
    c = b./sqrt(size(RMS_sime_index{1,k},1));
    
    title('RMS amplitude similarity index (july2020)');
    errorbar(a,c,colors{k});
    xticks([1,2,3,4,5,6,7]);
    set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
    xlim([0,6]);
    ylabel('Similarity index');
    legend('Untreated','Hyase');
    legend('boxoff');
    hold on
end

clear i k  a b c colors
save ('Figure2B(july2020)','RMS_sime_index' )
cd(path)
saveas(gcf,'Similarity index at BF RMSampl (July202)','pdf')
saveas(gcf,'Similarity index at BF RMSampl (July202)','fig')
%close all
%% RMS AvRecCSD and RElREsCSD (  This part of the script generate figure 3. in the paper Max, kindly add your script here)

RMS_AvRecCSD = {};  RMS_RelResCSD = {};

RMS_AvRecCSD{1} = nan(size(UBSSA_HBSSA{1},2),5);
RMS_AvRecCSD{2} = nan(size(UBSSA_HBSSA{1},2),5);

RMS_RelResCSD{1} = nan(size(UBSSA_HBSSA{1},2),5);
RMS_RelResCSD{2} = nan(size(UBSSA_HBSSA{1},2),5);

for cond = 1:size(UBSSA_HBSSA,2)
    for  anim = 1 : size(UBSSA_HBSSA{1,1},2)
        if BF_list(anim) == 1
            RMS_AvRecCSD{cond}(anim,3:5)= UBSSA_HBSSA{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,3:5)= UBSSA_HBSSA{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
        elseif BF_list(anim)== 2
            RMS_AvRecCSD{cond}(anim,2:5) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            RMS_RelResCSD{cond}(anim,2:5) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            
        elseif BF_list(anim) > 2 && BF_list(anim) < 7
            RMS_AvRecCSD{cond}(anim,:) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            RMS_RelResCSD{cond}(anim,:) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            
        elseif BF_list(anim)== 7
            RMS_AvRecCSD{cond}(anim,2:5) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,2:5) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
        else
            RMS_AvRecCSD{cond}(anim,[2,3,5]) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,[2,3,5]) = UBSSA_HBSSA{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)]);
        end
        
    end
end
cd(path)
save('RMS_AvRecCSD_july2020','RMS_AvRecCSD')
save('RMS_RelResCSD_july2020','RMS_RelResCSD')
clear cond anim 
namesRMSfigures = {'RMS of AVGRECCSD', 'RMS of RelResCSD'}
combineRMS = {RMS_AvRecCSD,RMS_RelResCSD} ; 
for i = 1: length(combineRMS)
    for k = 1 :size(combineRMS{i},2)
        meanRMS= nanmean(combineRMS{i}{k},1);
        stderror = nanstd(combineRMS{i}{k},0,1);
        subplot(1,2,i);
        errorbar(meanRMS,stderror)
        title(namesRMSfigures{i})
        xticks([1,2,3,4,5,6,7]);
        set(gca, 'xticklabel',{'BF-2'  'BF-1'  'BF'  'BF+1'  'BF+2' ''})
        xlim([0,6]);
        ylabel('RMS amplitude');
        legend('Untreated','Hyase');
        legend('boxoff');
        hold on;
    end
end

        
        
        







