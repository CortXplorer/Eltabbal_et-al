%% this script reproduce figure 6 A and B in the paper 
load ('SpontaneousECMDATA.mat');
UntreatedSpont= struct2cell(SpontaneousData{1,1})';
HYaseSpont= struct2cell(SpontaneousData{1,2})';
Untreated_Hyase = {UntreatedSpont, HYaseSpont };%Contains 2 cells w/ untreated and HYase data
eval_singletrial=0;    %if set true, each individual spontaneous columnar event detection will be shown

% Group data for granular ,early / late infragranular , supragranular layer traces
% resulting cell array: first 4 cells for each condition and in each condition four cortical layers w/ n animals
for i = 1 : size (Untreated_Hyase ,2) % for each condition
    for ii = 1: size (Untreated_Hyase{1,1}{1,1}.Binned_spontaneous_RMS_AvgSink , 1) % for each layer
        for iii = 1: size(Untreated_Hyase{1,1},2) % for each animal             
            % Create Spont_AvgSink w/ RMS values for {i}conditions and {iii} animals for {ii} layers 
            Spont_AvgSink{1,i}{iii,ii} = Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_RMS_AvgSink(ii,:);
        end 
    end    
end
clear i ii iii iiii  

% Group data for AvgRecCSD and RelResCSD data
for i = 1 : size (Untreated_Hyase ,2) % for each condition
    for ii = 1: size(Untreated_Hyase{1,1},2) % for each animal
            for iii = 1 :  size (Untreated_Hyase{1,1}{1,1}.Binned_spontaneous_RMS_AvgRecCSD, 2) %for each bin            
              % Create Spont_RelResAVREC w/ RMS values for 1.dim: (1) AVREC and
              % (2) RMS and 2 dim. {:,1} before and {:,2} after HYase for 7 animals
              Spont_RelResAVREC{1,i}{1,ii} = Untreated_Hyase{1,i}{1,ii}.Binned_spontaneous_RMS_AvgRecCSD;
              Spont_RelResAVREC{2,i}{1,ii} = Untreated_Hyase{1,i}{1,ii}.Binned_spontaneous_RMS_RelResCSD;
            end
    end
end
clear i ii iii iiii  



%% generating grand means before plotting 
% 1st step read out single data
for i = 1 : size(Spont_AvgSink,2) %for each condition
     for ii = 1 : size(Spont_AvgSink{1,1},2) % for each layer
         for iii = 1 : size(Spont_AvgSink{1,1},1) % for each animal   
             % Creates Mean_Spont_AvgSink w/ 2x4 cell for each cond and
             % layers for nsubXnrep each
             Mean_Spont_AvgSink{i,ii}(iii,:) = Spont_AvgSink{i}{iii,ii};
         end
     end
end
clear i ii iii                              

for i = 1 : size(Spont_RelResAVREC,1) %for each condition 
     for ii = 1 : size(Spont_RelResAVREC,2) % for each parameter
         for iii = 1 : size(Spont_RelResAVREC{1,1},2) % for each animal   
             % Creates Mean_Spont_RelResAVREC w/ 2x2 cell for each cond and
             % parameter for nsubXnrep each
             Mean_Spont_RelResAVREC{i,ii}(iii,:) = Spont_RelResAVREC{ii,i}{1,iii};
         end
     end
end
clear i ii iii  

% 2nd step paste all single data in one matrix
for i = 1 : size(Mean_Spont_AvgSink,1) %for each condition 
     for ii = 1 : size(Mean_Spont_AvgSink,2) % for each parameter
        mean_RMS_Sink{i,ii} = reshape(Mean_Spont_AvgSink{i,ii},1,[]);
     end
end

for i = 1 : size(Mean_Spont_RelResAVREC,1) %for each condition 
     for ii = 1 : size(Mean_Spont_RelResAVREC,2) % for each parameter
        mean_AVREC_Spont{i,ii} = reshape(Mean_Spont_RelResAVREC{ii,i},1,[]);
     end
end
clear i ii

% Plotting of RMS values
% Plot#1 - RMS of AvgRecCSD & RelResCSD (across rep. of animals)
titles = {'RMS AvgRecCSD', 'RMS RelResCSD'};
h_fig(1) = figure;
FIGURENAME = 'SpontActiv_RMS AvgRecCSD'
set(gcf,'Name',FIGURENAME ,'NumberTitle','off','position',[100,100,300,400]);
subplot(6,2,[1 3 5 7 9])
boxplot([mean(Mean_Spont_RelResAVREC{1,1},1)',mean(Mean_Spont_RelResAVREC{2,1},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2);
% delete outliers
h=findobj(gca,'tag','Outliers');delete(h);
%set ylim
%ylim([0.00525 .00625])
title(titles(1),'FontWeight','bold');
ylabel('RMS Amplitude of AvgRecCSD');
subplot(6,2,[2 4 6 8 10])
boxplot([mean(Mean_Spont_RelResAVREC{1,2},1)',mean(Mean_Spont_RelResAVREC{2,2},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2);
h=findobj(gca,'tag','Outliers');delete(h);
%ylim([0.0375 .05])
ylabel('RMS Amplitude of RelResCSD');
title(titles(2),'FontWeight','bold');

% Include statistics to plot
[n_AvgRecCSD,p_AvgRecCSD]=ttest(mean(Mean_Spont_RelResAVREC{1,1},1)',mean(Mean_Spont_RelResAVREC{2,1},1)');
[n_RelResCSD,p_RelResCSD]=ttest(mean(Mean_Spont_RelResAVREC{1,2},1)',mean(Mean_Spont_RelResAVREC{2,2},1)');
subplot(6,2,11);set(gca,'Visible','off')
text('Position',[0 0.5],'string','T-test; p-value AvgRecCSD:','FontWeight','bold')
text('Position',[0 0],'string',p_AvgRecCSD);
subplot(6,2,12);set(gca,'Visible','off')
text('Position',[0 0.5],'string','T-test; p-value RelResCSD:','FontWeight','bold')
text('Position',[0 0],'string',p_AvgRecCSD);

clear titles

% Plot#2 - RMS of AvgSinks 
titles = {'IG late', 'IG early', 'G', 'SG'};
h_fig(2) = figure;
FIGURENAME = 'SpontActiv_RMS AvgSink'
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,400,600,400]);
subplot(1,4,1)
boxplot([mean(Mean_Spont_AvgSink{1,1},1)',mean(Mean_Spont_AvgSink{2,1},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2)
% delete outliers
h=findobj(gca,'tag','Outliers');delete(h);
ylabel('RMS Amplitude of AvgSink');
ylim([0.0002 0.0007]);
title(titles(1),'FontWeight','bold');
% plot t-Test results
[n_IGlate,p_IGlate]=ttest(mean(Mean_Spont_AvgSink{1,1},1)',mean(Mean_Spont_AvgSink{2,1},1)');
annotation('textbox',[0.05 0 0.1 0.1],'string','T-Test:','EdgeColor','none','FontWeight','bold');
annotation('textbox',[0.12 0 0.1 0.1],'string',p_IGlate,'EdgeColor','none','FontWeight','bold');

subplot(1,4,2)
boxplot([mean(Mean_Spont_AvgSink{1,2},1)',mean(Mean_Spont_AvgSink{2,2},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2)
% delete outliers
h=findobj(gca,'tag','Outliers');delete(h);
ylim([0.0002 0.0007]);
title(titles(2),'FontWeight','bold');
% plot t-Test results
[n_IGearly,p_IGearly]=ttest(mean(Mean_Spont_AvgSink{1,2},1)',mean(Mean_Spont_AvgSink{2,2},1)');
annotation('textbox',[0.33 0 0.1 0.1],'string',p_IGearly,'EdgeColor','none','FontWeight','bold');

subplot(1,4,3)
boxplot([mean(Mean_Spont_AvgSink{1,3},1)',mean(Mean_Spont_AvgSink{2,3},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2)
% delete outliers
h=findobj(gca,'tag','Outliers');delete(h);
ylim([0.0002 0.0007]);
title(titles(3),'FontWeight','bold');
% plot t-Test results
[n_G,p_G]=ttest(mean(Mean_Spont_AvgSink{1,3},1)',mean(Mean_Spont_AvgSink{2,3},1)');
annotation('textbox',[0.55 0 0.1 0.1],'string',p_G,'EdgeColor','none','FontWeight','bold','LineWidth',2);

subplot(1,4,4)
boxplot([mean(Mean_Spont_AvgSink{1,4},1)',mean(Mean_Spont_AvgSink{2,4},1)'],'plotstyle','compact','color',['b';'r'],'widths',1);
set(gca,'XTick',1:2,'XTickLabel',{'Untreated','HYase'},'Color','none','box','off','LineWidth',2)
% delete outliers
h=findobj(gca,'tag','Outliers');delete(h);
ylim([0.0002 0.0007]);
title(titles(4),'FontWeight','bold');
% plot t-Test results
[n_SG,p_SG]=ttest(mean(Mean_Spont_AvgSink{1,4},1)',mean(Mean_Spont_AvgSink{2,4},1)');
annotation('textbox',[0.75 0 0.1 0.1],'string',p_SG,'EdgeColor','none','FontWeight','bold');

%% Qunatification of Spontaneous Columnar Events (SCE's)

SpontColEvtUntreated=[];%Counting variable for SCE's
SpontColEvtHYase=[];%Counting variable for SCE's
for i = 1 : size (Untreated_Hyase ,2) % for each condition
    for ii = 1: size (Untreated_Hyase{1,1}{1,1}.Binned_spontaneous_AvgRecCSD , 2) % for each trial
        for iii = 1: size(Untreated_Hyase{1,1},2) % for each animal             
            % Create Spont_AvgSink w/ RMS values for {i}conditions and {iii} animals for {ii} layers 
            % Spont_AvgSink{1,i}{iii,ii} = Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_RMS_AvgSink(ii,:);
            
                SingleTrial_AvgRecCSD{i}{iii,ii}=Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_AvgRecCSD{ii};
                peaks = Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_AvgRecCSD{ii};
                a = mean(peaks);a = a + 3*std(peaks);
                [pks,locs] = findpeaks(peaks,'minpeakheight',a,'minpeakdistance',150);        % Find peaks and their indices
                if i==1
                    SpontColEvtUntreated=[SpontColEvtUntreated;size(pks,2)];
                else
                    SpontColEvtHYase=[SpontColEvtHYase;size(pks,2)];
                end

                
              if i==2
                
              if eval_singletrial
              % Plot individual example
                        h_fig=figure
                        subplot(2,1,1)
                        imagesc(Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_CSD{ii});hold on
                        subplot(2,1,2)
                        plot(SingleTrial_AvgRecCSD{i}{iii,ii});hold on;
                        plot(locs,peaks(locs),'k^','markerfacecolor',[1 0 0]);
                        keyboard
                        close(h_fig)
              end
              end
        end 
    end 
end


figure
set(gcf,'position',[100,100,900,900]);
set(gcf,'Name','SpontaneousColumnarEvents' ,'NumberTitle','off') ;
colors = {'b','r'};
FIGURENAME='SpontaneousColumnarEvents';

% % % choose good single trial example to plot
% Untreated
      i=1;ii=14;iii=2;
subplot(6,5,[1:4 6:9])
    imagesc(Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_CSD{ii});hold on
    box off
    clear clim
    clim=get(gca,'CLim')*.75;
    caxis([-max(clim) max(clim)]);
    title('Spont. Column. Events (SCE)','fontsize',12) 
    ylabel('Cortical depth [mm]','fontsize',12)
    set(gca, 'YTick', [0 8 16 24])
    set(gca, 'yticklabel',[0 .5 1 1.5])
    set(gca,'fontsize',10,'color','none','LineWidth',2);
subplot(6,5,[11:14])   
    peaks = Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_AvgRecCSD{ii};
    a = median(peaks);a = a + 3*std(peaks);
    [pks,locs] = findpeaks(peaks,'minpeakheight',a,'minpeakdistance',100);  
    plot(SingleTrial_AvgRecCSD{i}{iii,ii});hold on;
    ylim([0 max(SingleTrial_AvgRecCSD{i}{iii,ii})]);
    set(gca,'fontsize',10,'color','none','LineWidth',2);
    box off
    plot(locs,peaks(locs),'k^','markerfacecolor',[1 0 0]);      
    %xlabel('Time [ms]','fontsize',12)
    ylabel('AvgRecCSD [mV/mm^2]','fontsize',12)


     i=2;ii=12;iii=2; % this is good i=2;ii=2;iii=3

subplot(6,5,[16:19 21:24])
    imagesc(Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_CSD{ii});hold on
    box off
    %clim=get(gca,'CLim')*.75;
    caxis([-max(clim) max(clim)]);
    %title('Spont. Column. Events (SCE) after HYase') 
    ylabel('Cortical depth [mm]','fontsize',12)
    set(gca, 'YTick', [0 8 16 24])
    set(gca, 'yticklabel',[0 .5 1 1.5])
    set(gca,'fontsize',10,'color','none','LineWidth',2);
subplot(6,5,[26:29])   
    peaks = Untreated_Hyase{1,i}{1,iii}.Binned_spontaneous_AvgRecCSD{ii};
    a = median(peaks);a = a + 3*std(peaks);
    [pks,locs] = findpeaks(peaks,'minpeakheight',a,'minpeakdistance',100);  
    plot(SingleTrial_AvgRecCSD{i}{iii,ii});hold on;
    ylim([0 max(SingleTrial_AvgRecCSD{i}{iii,ii})]);
    set(gca,'fontsize',10,'color','none','LineWidth',2);
    box off
    plot(locs,peaks(locs),'k^','markerfacecolor',[1 0 0]);      
    xlabel('Time [ms]','fontsize',12)
    ylabel('AvgRecCSD [mV/mm^2]','fontsize',12)
    
subplot(6,5,[5 10])
% divide sum per 3sec by 3 to get SCE's / sec
boxplot([SpontColEvtUntreated/6,SpontColEvtHYase/6],'plotstyle','compact','color',['b';'r'],'widths',1,'whisker',0.75);
box off
title('SCE / sec')
set(gca,'fontsize',12,'color','none','LineWidth',2);%,'XTickLabel',{'Untreated','Hyase'})
set(gca,'XTick',[1 2],'XTickLabels',{'un.' 'HY'});
% TTest
[n_SCE,p_SCE]=ttest(SpontColEvtUntreated',SpontColEvtHYase');
annotation('textbox',[0.8 0.55 0.1 0.1],'string','T-Test:','EdgeColor','none','FontWeight','bold');
annotation('textbox',[0.8 0.52 0.1 0.1],'string',p_SCE,'EdgeColor','none','color','r');

% Insert labels
annotation('textbox',[0.128 0.845 0.1 0.1],'string','Untreated','EdgeColor','none','fontsize',8,'FontWeight','bold');%,'BackgroundColor',[.8 .8 .8])
annotation('textbox',[0.128 0.42 0.1 0.1],'string','HYase','EdgeColor','none','fontsize',8,'FontWeight','bold');%,'BackgroundColor',[.8 .8 .8])

