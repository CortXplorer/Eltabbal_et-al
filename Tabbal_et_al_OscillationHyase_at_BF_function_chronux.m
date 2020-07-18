function Tabbal_et_al_OscillationHyase_at_BF_function00
%this functiongroup different variables and reproduce figure 4 in the paper
%

path = '/Users/tabbal/Dropbox (OIST)/CorteXexplorer LAB/Data Organized/Hyase_ECM data/Figures';
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
%% LFP analysis in this part i will generate power spectrograms for granular ,early / late infragranular , supragranular
% the resulting cell array , the first 4 cells for each condition and in
% each condition there are four cortical layers , in each u find n animals
% traces
for k = 1 : size (UBSSA_HBSSA ,2) % for each condition
    for ii = 1: size (UBSSA_HBSSA{1,1}{1,1}.AvgSink , 3) % for each layer
        for i = 1: size(UBSSA_HBSSA{k},2) % for each animal
            all_animals_BF_avgsink{1,k}{1,ii}(:,i) = UBSSA_HBSSA{1,k}{1,i}.AvgSink(:,BF_list(i),ii) ;
            if BF_list(i) > 6
                all_animals_SSA_BF_avgsink{1,k}{1,ii}(:,i)= UBSSA_HBSSA{1,k}{1,i}.AvgSink(:,BF_list(i)-2,ii) ;
            else
                all_animals_SSA_BF_avgsink{1,k}{1,ii}(:,i)= UBSSA_HBSSA{1,k}{1,i}.AvgSink(:,BF_list(i)+2,ii) ;
            end
        end
    end
end
clear i ii k HASSA HBSSA UASSA UBSSA FreqTunDataContainer 
all_animals_BF_avgsink_hyase = {all_animals_BF_avgsink{1,1} , all_animals_BF_avgsink{1,2}} ;
clear all_animals_SSA_BF_avgsink all_animals_BF_avgsink
cd(path)
save('all_animals_BF_avgsink_hyase_July2020','all_animals_BF_avgsink_hyase')
%% generating a spectrogram and power spectrum using Chronux functions ( Multitaper analysis based on chronux toolbox
% you will need to use chronux toolbox and make sure that you remove MVG
% Grangercausality toolbox cuz there are some conflict between some
% functions this part of the script will generate figure 4A 
params.tapers=[3 5]; params.pad= 0; params.Fs=1000; params.fpass=[1 100]; params.err = [1 0.05]; params.trialave=0;
for layer = 1 : size (all_animals_BF_avgsink_hyase{1,1},2)
    for cond = 1 : size (all_animals_BF_avgsink_hyase ,2)
        [power_spec{1,cond}{1,layer},f2]= mtspectrumc(all_animals_BF_avgsink_hyase{1,cond}{1,layer},params);
    end
end
clear layer  cond
% normalized power spectrum
for cond = 1 : size (all_animals_BF_avgsink_hyase,2)
    for layer = 1: size (all_animals_BF_avgsink_hyase{1,1} ,2)
        for k = 1 : size (power_spec{1,cond}{1,layer},2)
        norm_power_spec{1,cond}{1,layer}(:,k)= power_spec{1,cond}{1,layer}(:,k)./sum(power_spec{1,cond}{1,layer}(:,k));
        end
    end
end
clear cond layer 
%plotting
h_fig = figure;
FIGURENAME = 'powerspectrum Chronux version';
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
colors = {'b','r'};
layers  = {'Late infragranular layer', 'Early infragranular layer' , 'Granular layer', 'Supragranular layer'};
for layer = 1 : size (all_animals_BF_avgsink_hyase{1,1} ,2);
    for cond = 1 : size (all_animals_BF_avgsink_hyase,2)
        meanS = mean (norm_power_spec{1,cond}{1,layer},2);
        Serror = std(norm_power_spec{1,cond}{1,layer},0,2)./sqrt(size (power_spec{1,cond}{1,layer},2));
        subplot(2,2,layer);
        shadedErrorBar(f2,meanS,Serror,colors{cond},[0]);
        hold on
    end 
    title(layers{layer});
    set(gca,'YScale','log')
    ylabel('Normalized power')
    xlabel('Freq (Hz)')
    xlim([0 101])
    ylim([10e-6 10e-1])
end 

legend('Untreated','','','','Hyase')
set(legend,'Location','best');
cd(path)
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
saveas(h_fig,FIGURENAME,'pdf')
saveas(h_fig,FIGURENAME,'fig')
save('EVOKEDPSD_10thFeb2020','norm_power_spec')
clear layers cond k  layer meanS params k colors Serror Serr1 screenposition
%% binning power spectrum ( this is more artificial way of determining bands we did this analysis based on conventional definition of bands, however
% we thought it is more accurate to depart this artifical bands selections
% into a more statistically unbiased detection of the bands that showed
% siginificant differences. 
freq_names = {'Delta band', 'Theta band', 'alpha band','beta band' , 'low Gamma band ' ,'High-Gamma band'};
dlg_title = 'Freq. bands'
num_lines = 1 ;
default = {'[1,4] ', '[4,10]', '[10,20]','[20,30]','[30,60]' ,'[60,90]'};
Freq_bands = inputdlg(freq_names,dlg_title, num_lines , default)
clear num_lines dlg_title default

Delta = eval(Freq_bands{1}); Delta = [min(find(f2>=Delta(1)&f2<Delta(1)+1)),min(find(f2>=Delta(2)&f2<Delta(2)+1))];
Theta = eval(Freq_bands{2}); Theta = [min(find(f2>=Theta(1)&f2<Theta(1)+1)),min(find(f2>=Theta(2)&f2<Theta(2)+1))];
alpha = eval(Freq_bands{3}); alpha = [min(find(f2>=alpha(1)&f2<alpha(1)+1)),min(find(f2>=alpha(2)&f2<alpha(2)+1))];
beta = eval(Freq_bands{4});  beta = [min(find(f2>=beta(1)&f2<beta(1)+1)),min(find(f2>=beta(2)&f2<beta(2)+1))];
low_gamma = eval(Freq_bands{5});low_gamma = [min(find(f2>=low_gamma(1)&f2<low_gamma(1)+1)),min(find(f2>=low_gamma(2)&f2<low_gamma(2)+1))];
high_gamma =eval(Freq_bands{6}); high_gamma = [min(find(f2>=high_gamma(1)&f2<high_gamma(1)+1)),min(find(f2>=high_gamma(2)&f2<high_gamma(2)+1))];
% 
for cond = 1: size (all_animals_BF_avgsink_hyase,2) %for each condition
    for layer = 1 : size (all_animals_BF_avgsink_hyase {1,1},2)
        data_binned.delta{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(Delta(1):Delta(2),:);
        data_binned.theta{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(Theta(1):Theta(2),:);
        data_binned.alpha{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(alpha(1):alpha(2),:);
        data_binned.beta{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(beta(1):beta(2),:);
        data_binned.lowgamma{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(low_gamma(1):low_gamma(2),:);
        data_binned.highgamma{1,cond}{1,layer}= norm_power_spec{1,cond}{1,layer}(high_gamma(1):high_gamma(2)-1,:);
    end 
end
clear layer cond 
data_binned = {data_binned.delta, data_binned.theta, data_binned.alpha, data_binned.beta, data_binned.lowgamma, data_binned.highgamma};
for freq = 1 : size(data_binned,2)
    for cond = 1 :size(data_binned{1,1},2)
        for lay = 1 : size(data_binned{1,1}{1,1},2)
            data_binned_plot{1,freq}{1,cond}(:,lay) = sum(data_binned{1,freq}{1,cond}{1,lay},1)';
        end
    end
end 
clear freq cond lay 
% plotting each frequency band (2 conditions 4 layers )

h_fig = figure
FIGURENAME = 'binned powerspectrum Chronux version' 
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
colors = {'b','r'}
layers  = {'L-IG', 'E-IG' , 'G', 'SG'}
for freq = 1 : size (data_binned_plot,2)
    subplot(2,3,freq)
    for cond = 1 : size(data_binned_plot{1,1},2)
    a = mean(data_binned_plot{1,freq}{1,cond},1);
    serr = std (data_binned_plot{1,freq}{1,cond},1)./sqrt(size(data_binned_plot{1,freq}{1,cond},1));
    errorbar(a,serr,colors{cond});
    hold on 
    end
    ax = gca;
    set(gca, 'XTick', [1 2 3 4])
    set(gca, 'xticklabel',layers,'fontsize',12);
    ax.XTickMode = 'manual';
    xlabel('cortical layer','fontsize',12)
    ylabel('Normlaized Power','fontsize',12);
    title(strcat(freq_names{freq}, ' - power spectrum '),'fontsize',12);
end
cd(path)
legend('Untreated','Hyase')
set(gcf,'Units','inches');
saveas(gcf,FIGURENAME,'fig');
%close all 
%% Ttest for each frequency bin frequencies reproduce figure4B 
path = '/Users/tabbal/Dropbox (OIST)/CorteXexplorer LAB/Data Organized/Hyase_ECM data/Figures';
EvokedPSD = load('EVOKEDPSD_10thFeb2020.mat'); EvokedPSD = EvokedPSD.norm_power_spec;
SpontPSD = load('SpontaneousPSD_10thFEB202.mat'); SpontPSD = SpontPSD.mean_trials_norm_power_spec;

ALLPSD = {EvokedPSD,SpontPSD}; 
% adding one step to take the logarithm of the normalized spectrum 
for psd = 1 :size(ALLPSD,2)
    for cond = 1 :size(ALLPSD{psd},2)
        for layer = 1 :size(ALLPSD{psd}{cond},2)
            AllPSDlog{psd}{cond}{layer} = log(ALLPSD{psd}{cond}{layer});
        end
    end
end
EvokedPSD = AllPSDlog {1}; SpontPSD = AllPSDlog {2}; 
clear layer psd cond 
layers  = {'Late infragranular layer', 'Early infragranular layer' , 'Granular layer', 'Supragranular layer'}
h_fig = figure
FIGURENAME = 'Ttest Evoked Vs Spontaneous  for each frequecy bin'
set(h_fig,'Units','inches')
fig_pos = get(h_fig,'position');
fig_margin = 0.01;
set(h_fig,'papersize',fig_pos(3:4)*[1+fig_margin*2],'paperposition',[fig_pos(3:4)*fig_margin fig_pos(3:4)]);
set(gcf,'Name',FIGURENAME ,'NumberTitle','off')
set(gcf,'position',[100,100,1300,700])
colors={'g','r'}
f2 = 1:size(EvokedPSD{1}{1},1);
for layer = 1 : size (EvokedPSD{1,1} ,2);
    for freqbin = 1 :size (EvokedPSD{1}{1},1)
        [Eh,Ep,ci,Estats] = ttest(EvokedPSD{1,1}{1,layer}(freqbin,:),EvokedPSD{1,2}{1,layer}(freqbin,:));
        [Sh,Sp,ci,Sstats] = ttest(SpontPSD{1,1}{1,layer}(freqbin,:),SpontPSD{1,2}{1,layer}(freqbin,:));
        

        Etvalue{layer}(freqbin) = Estats.tstat;
        Ehvec{layer}(freqbin) = Eh ;Epvec{layer}(freqbin)= Ep;
        
        Stvalue{layer}(freqbin) = Sstats.tstat;
        Shvec{layer}(freqbin) = Sh ; Spvec{layer}(freqbin)= Sp; 
        
       
        clear Istats Ih Ip Estats Eh Ep Sstats Sh Sp
    end
    [VectSh, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Spvec{layer},0.05,'pdep','yes');
    [VectEh, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Epvec{layer},0.05,'pdep','yes');
    
    subplot(2,2,layer);
    plot(f2,Etvalue{layer},'g')
    hold on
    plot(f2,Stvalue{layer},'r')
    hold on 
    Eindex= find(VectEh==1);  Sindex= find(VectSh==1); 
    if ~isempty(Eindex)
        y = repmat([8.5],size(f2(Eindex)));
        scatter(f2(Eindex),y,'g','.')
        clear y
    else
    end
    
    if ~isempty(Sindex)
        y = repmat([6.5],size(f2(Sindex)));
        scatter(f2(Sindex),y,'r','.')
        clear y
    else
    end    
    
    title(layers{layer});
    ylabel('ttest value')
    xlabel('Freq (Hz)')
    xlim([0 101])
    ylim([-10 10])
    legend('eVokedPSD','SpontaneousPSD','location','southeast')
end
cd(path)
saveas(gcf,FIGURENAME,'fig');
saveas(gcf,FIGURENAME,'pdf');
%close all
clear FIGURENAME ci colors Eh Ehvec Eindex  h_fig Ipvec Epvec layer layers f2 Etvalue  meantrials_norm_power_spec norm_power_spec path fig_pos fig_margin freqbin f2


end