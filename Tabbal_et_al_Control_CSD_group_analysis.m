function Tabbal_Control_CSD_group_analysis
%this function reproduce supplementary figure 2 
p_cur = CleanSlate_00('Tabbal_Control_CSD_group_analysis','open')
%%
load ('ControlECMDATA.mat');
BeforeNaCl= struct2cell(FreqTunDataContainer {1,1})';
AfterNaCl= struct2cell(FreqTunDataContainer {1,2})';


for i = 1 : size(BeforeNaCl,2)
    BF_list(i) = BeforeNaCl{1,i}.BF_table(3); 
end
layers = {'late_infragranular' , 'early_infragranular' , 'granular' , 'supragraular'};
clear i
%% the effect of nacl on CSD  RMS_AvgSink
RMS_AvgSink_before_NaCl_matrix = nan(size(BeforeNaCl,2), size(BeforeNaCl{1}.AvgSink_Peak,2), size(BeforeNaCl{1}.AvgSink_Peak,3)) ;
RMS_AvgSink_after_NaCl_matrix = nan(size(AfterNaCl,2), size(AfterNaCl{1}.AvgSink_Peak,2), size(AfterNaCl{1}.AvgSink_Peak,3)) ;
RMS_AvgSinkBAmatrix = {RMS_AvgSink_before_NaCl_matrix,RMS_AvgSink_after_NaCl_matrix};
BeforeNaCl_AfterNaCl = {BeforeNaCl, AfterNaCl};
for i = 1:size(RMS_AvgSinkBAmatrix,2)
    for ii = 1:size(BeforeNaCl,2) % number of animals
        for k = 1:size(BeforeNaCl{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if ii==2
                RMS_AvgSinkBAmatrix{1,i}(ii,2:end,k) = BeforeNaCl_AfterNaCl{1,i}{1,ii}.RMS_AvgSink(1,:,k);
            else
                RMS_AvgSinkBAmatrix{1,i}(ii,:,k) = BeforeNaCl_AfterNaCl{1,i}{1,ii}.RMS_AvgSink(1,:,k);
            end
        end
    end
end
clear i ii k
BF_nonBF_before_NaCl_matrix= nan(size(BeforeNaCl,2),5,size(BeforeNaCl{1}.AvgSink_Peak,3));
BF_nonBF_after_NaCl_matrix= nan(size(BeforeNaCl,2),5,size(BeforeNaCl{1}.AvgSink_Peak,3));
BF_nonBF_RMS_AvgSinkBAmatrix = {BF_nonBF_before_NaCl_matrix,BF_nonBF_after_NaCl_matrix};
for i = 1:size(BF_nonBF_RMS_AvgSinkBAmatrix,2)
    for ii = 1:size(BeforeNaCl,2) % number of animals
        for k = 1:size(BeforeNaCl{1}.AvgSink_Peak,3); % k number of sinks 1-late infra 2-early infra 3-granular 4-supra
            if BF_list(ii) == 1
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(ii,[3:5],k) = RMS_AvgSinkBAmatrix{1,i}(ii,[BF_list(ii),BF_list(ii)+1,BF_list(ii)+2],k);
            elseif BF_list(ii)== 2
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(ii,[2:5],k) = RMS_AvgSinkBAmatrix{1,i}(ii,[BF_list(ii)-1,BF_list(ii),BF_list(ii)+1,BF_list(ii)+2],k);
            elseif BF_list(ii) > 2 && BF_list(ii) < 7
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(ii,:,k) = RMS_AvgSinkBAmatrix{1,i}(ii,[BF_list(ii)-2,BF_list(ii)-1,BF_list(ii),BF_list(ii)+1,BF_list(ii)+2],k);
            elseif BF_list(ii)== 7
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(ii,[2:5],k) = RMS_AvgSinkBAmatrix{1,i}(ii,[BF_list(ii)-1,BF_list(ii),BF_list(ii)+1,BF_list(ii)-2],k);
            else
                BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(ii,[2,3,5],k) = RMS_AvgSinkBAmatrix{1,i}(ii,[BF_list(ii)-1,BF_list(ii),BF_list(ii)-2],k);
                
            end
        end
    end
end
clear i ii k RMS_AvgSink_after_hyase_matrix RMS_AvgSink_before_hyase_matrix RMS_AvgSinkBAmatrix BF_nonBF_after_hyse_matrix BF_nonBF_before_hyse_matrix

% ploting RMS_AvgSink_ tuning before and after hyase in different layers
colors = {'b','r'} ;
figure
FIGURENAME = 'RMS_AvgSink TUNING CURVES BEFORE AND AFTER INJECTION '
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
for k = 1:size(BF_nonBF_RMS_AvgSinkBAmatrix{1,1},3);
    for i = 1 : size(BF_nonBF_RMS_AvgSinkBAmatrix,2);
        a= nanmean(BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(:,:,k));
        b= nanstd(BF_nonBF_RMS_AvgSinkBAmatrix{1,i}(:,:,k));
        c = b./sqrt(size(BeforeNaCl,2));
        subplot(2,2,k);
        title(layers{k});
        errorbar(a,c,colors{i});
        set(gca, 'xticklabel',{'' '-2BF' '-1Bf' 'BF' '+1BF' 'SSA' ''})
        ylabel('RMS-Avg-Sink amplitude');
        xlim([0,6])
        legend('Before NaCl','After NaCl' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
cd(p_cur)
set(gcf,'Units','inches')
fig_pos = get(gcf,'position');
fig_margin = 0.01;
set(gcf,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
saveas(gcf,FIGURENAME,'pdf')
saveas(gcf,FIGURENAME,'fig')
save('RMS_AvgSinkBAmatrix.mat','BF_nonBF_RMS_AvgSinkBAmatrix');
%close all
clear a b c colors i k  FIGURENAME
clear BF_nonBF_after_NaCl_matrix BF_nonBF_before_NaCl_matrix RMS_AvgSinkBAmatrix fig_margin fig_pos RMS_AvgSink_after_NaCl_matrix RMS_AvgSink_before_NaCl_matrix BF_nonBF_RMS_AvgSinkBAmatrix

%% RMS AvRecCSD and RElREsCSD
RMS_AvRecCSD = {};  RMS_RelResCSD = {};

RMS_AvRecCSD{1} = nan(size(BeforeNaCl_AfterNaCl{1},2),5);
RMS_AvRecCSD{2} = nan(size(BeforeNaCl_AfterNaCl{1},2),5);

RMS_RelResCSD{1} = nan(size(BeforeNaCl_AfterNaCl{1},2),5);
RMS_RelResCSD{2} = nan(size(BeforeNaCl_AfterNaCl{1},2),5);

for cond = 1:size(BeforeNaCl_AfterNaCl,2)
    for  anim = 1 : size(BeforeNaCl_AfterNaCl{1,1},2)
        if BF_list(anim) == 1
            RMS_AvRecCSD{cond}(anim,3:5)= BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,3:5)= BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
        elseif BF_list(anim)== 2
            RMS_AvRecCSD{cond}(anim,2:5) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            RMS_RelResCSD{cond}(anim,2:5) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            
        elseif BF_list(anim) > 2 && BF_list(anim) < 7
            RMS_AvRecCSD{cond}(anim,:) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            RMS_RelResCSD{cond}(anim,:) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-2,BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)+2]);
            
        elseif BF_list(anim)== 7
            RMS_AvRecCSD{cond}(anim,2:5) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,2:5) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)+1,BF_list(anim)]);
        else
            RMS_AvRecCSD{cond}(anim,[2,3,5]) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_AvgRecCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)]);
            RMS_RelResCSD{cond}(anim,[2,3,5]) = BeforeNaCl_AfterNaCl{1,cond}{1,anim}.RMS_RelResCSD([BF_list(anim)-1,BF_list(anim),BF_list(anim)]);
        end
        
    end
end
save('RMS_AvRecCSD_Control','RMS_AvRecCSD');
save('RMS_RelResCSD_Control','RMS_RelResCSD');
Avgrec_RelResCSD = {RMS_AvRecCSD,RMS_RelResCSD};
titleavrec = {'RMS-AvRecCSD','RMS-RelResCSD'} ;
figure
colors = {'b' , 'r'};
for k = 1:size(Avgrec_RelResCSD,2);
    for i = 1 : size(Avgrec_RelResCSD{k},2);
        a= nanmean(Avgrec_RelResCSD{k}{i});
        b= nanstd(Avgrec_RelResCSD{k}{i});
        c = b./sqrt(size(BeforeNaCl,2));
        subplot(2,1,k);
        title(titleavrec{k});
        errorbar(a,c,colors{i});
        set(gca, 'xticklabel',{'' '-2BF' '-1Bf' 'BF' '+1BF' 'SSA' ''})
        ylabel('RMS-AVrec or RelResCSD');
        xlim([0,6])
        legend('Before NaCl','After NaCl' ,'Location','southwest','FontSize',6);
        legend('boxoff');
        hold on
    end
end
end

