%% This is a script to analyse bivariate/multivariate granger predictive causality
% this is the singletrial  for data on the BF stimulation  this script will
% reproduce figure 5 with all panels
%% load LFP data, here the format of the resulting matrix X will be two cell for each coniditon 1st is the untreated
%seond is the Hyase, in each there will be matrix ( channels x time x
%trials (animals) )
%cd('/Users/tabbal/Dropbox (OIST)/CorteXexplorer LAB/Data Organized/Granger causality')
%fig_path = 'C:\Users\TABBAL\Dropbox (OIST)\CorteXexplorer LAB\Data Organized\Granger causality'
load ('all_animals_BF_singlesink_hyase_10thFeb2020.mat');
mkdir GrangerCausality;
fig_path = [cd, '\GrangerCausality'];

%% detrending data using the difference method
nData = all_animals_BF_singlesink_hyase;
for cond = 1 : size(nData,2) %condition
    for lay = 1: size(nData{1},2) %layer
        for anim = 1 : size(nData{1}{1},2) % animal
            for trial = 1 :size(nData{cond}{lay}{anim},1) %trial
                detrendeddata{cond}{anim}(lay,:,trial)= diff(nData{cond}{lay}{anim}(trial,:));
            end
        end
    end
end
clear lay anim  cond i k trial
%%
%stationary_matrices = cell(size(detrendeddata,2),size(detrendeddata{1},2))
for cond = 1 :size(detrendeddata,2)
    for anim = 1:size(detrendeddata{1},2)
        for trial = 1:size(detrendeddata{cond}{anim},3)
            sd(cond,trial) = check_stationarity(detrendeddata{cond}{anim}(:,:,trial));
        end
    end
end
clear  i k
clear sd s ans
save('SingleTrialBFdetrendeddata','detrendeddata'); save('SingleTrialBFData','nData');
%% deternded data use (100 ms to 500 ms around stimulus) 
%% deternded data use (100 ms to 500 ms around stimulus)
%avg_trials = 1;
for cond = 1 :size(detrendeddata,2)
    for anim = 1:size(detrendeddata{1},2)
        for trial = 1:size(detrendeddata{cond}{anim},3)  %/avg_trials)-1
%             k =size(detrendeddata{cond}{anim},3);
%             r = randi(k,[1,avg_trials]);
             %firstpoint = trial*avg_trials+1; secondpoint = firstpoint+(avg_trials-1) ; 
            detrendeddata_monte{cond}{anim}(:,:,trial) = detrendeddata{cond}{anim}(:,100:499,trial); %replace trial with r and add mean or (trial with commenting r and k) 100:499
        end
    end
end
clear cond anim k r avg_trials r trial detrendeddata
detrendeddata = detrendeddata_monte;
%% declaring some variables to be used in the mvgc toolbox
nvars = 4 % number of channels
icregmode = 'LWR'; regmode   = 'LWR'; morder    = 'AIC'; momax     = 20; fs = 1000;
nobs      = 400;  % number of datapoints
ntrials   = 0;    % in this case will change
tstat     = '';
alpha     = 0.05; mhtc ='FDR' ; acmaxlags=400; fres=100; specm = 'MT'; etests= true; stlags=[]; acorr=true;
nperms    = 2; % for permutations you dont need it as you will generate model for each animal in different conditions and then compare
bsize = []; cm =flipud(bone);
%%
F = cell(size(detrendeddata,2),1)
F_spec = cell(size(detrendeddata,2),1)
for anim = 1:size(detrendeddata{1},2)
    for cond = 1 :size(detrendeddata,2)
        X = detrendeddata{cond}{anim};
        ntrials = size(detrendeddata{cond}{anim},3)
        % defining the model order to be used in the AR model using AIC test used funciton from mvgc toolbox
        ptic('\n*** tsdata_to_infocrit\n');
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        ptoc('*** tsdata_to_infocrit took ');
        
        % Plot information criteria.
        
        [~,bmo_AIC] = min(AIC);
        [~,bmo_BIC] = min(BIC);
        
        figure(1); clf;
        plot((1:momax)',[AIC BIC]);
        title('Model order estimation');
        legend('AIC','BIC');
        
        fprintf('\nbest model order (AIC) = %d\n',moAIC);
        fprintf('best model order (BIC) = %d\n',moBIC);
        
        if morder == 'AIC'
            morder=bmo_AIC;
        else
            morder = bmo_BIC
        end
        
        %% Estimating the Covaiance matrix and coefficients of the vector autoregression model(VAR model) using the mvgc toolbox
        ptic('\n*** tsdata_to_var... ');
        [A,SIG,E] = tsdata_to_var(X,morder,regmode);
        ptoc;
        %assert(~isbad(A),'VAR estimation failed'); %qualtiy check if regression was done or not
        
        ptic('*** var_to_autocov... ');
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        ptoc;
        
        var_info(info,true); % Report results and check for errors.
        

        %% getting the pairwise granger causality of the model
        
        ptic('*** autocov_to_pwcgc... ');
        F{cond}{anim} = autocov_to_pwcgc(G);
        ptoc;
        % Check for failed GC calculation
        assert(~isbad(F{cond}{anim},false),'GC calculation failed');
        %% spectral
        % Calculate spectral pairwise-conditional causalities at given frequency
        % resolution - again, this only requires the autocovariance sequence.
        
        ptic('\n*** autocov_to_spwcgc... ');
        F_spec{cond}{anim} = autocov_to_spwcgc(G,fres);
        ptoc;
        
        % Check for failed spectral GC calculation
        
        assert(~isbad(F_spec{cond}{anim},false),'spectral GC calculation failed');
        
        % Plot spectral causal graph.
        
        figure(3); clf;
        plot_spw(F_spec{cond}{anim},fs);
        
        %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)
        
        % Check that spectral causalities average (integrate) to time-domain
        % causalities, as they should according to theory.
        
        fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
        Fint = smvgc_to_mvgc(F_spec{cond}{anim}); % integrate spectral MVGCs
        mad = maxabs(F{cond}{anim}-Fint);
        madthreshold = 1e-5;
        if mad < madthreshold
            fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
        else
            fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
        end
        %
    end
    
    
end
%%
x = 1:101; colors = {'b','r'}; figure; 
for cond = 1:size(F_spec,1)
    for anim = 1:size(F_spec{cond},2)
        lateinfra_earlyinfra(cond,:,anim) = (F_spec{cond}{anim}(1,2,:));
        lateinfra_gran(cond,:,anim) = (F_spec{cond}{anim}(3,2,:));
        lateinfra_supra(cond,:,anim) = (F_spec{cond}{anim}(4,2,:));
        
        earlyinfra_supra(cond,:,anim) = (F_spec{cond}{anim}(4,1,:));
        earlyinfra_gran(cond,:,anim) = (F_spec{cond}{anim}(3,1,:));
        earlyinfra_lateinfra(cond,:,anim) = (F_spec{cond}{anim}(2,1,:));
        
        supra_earlyinfra(cond,:,anim) = (F_spec{cond}{anim}(1,4,:));
        supra_lateinfra(cond,:,anim) = (F_spec{cond}{anim}(2,4,:));
        supra_granular(cond,:,anim) = (F_spec{cond}{anim}(3,4,:));
        
        gran_earlyinfra(cond,:,anim) = (F_spec{cond}{anim}(1,3,:));
        gran_lateinfra(cond,:,anim) = (F_spec{cond}{anim}(2,3,:));
        gran_supra(cond,:,anim) = (F_spec{cond}{anim}(4,3,:));
        
    end
    subplot(3,4,1)
    y = mean(lateinfra_earlyinfra(cond,:,:),3);
    sem = std(lateinfra_earlyinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
    title('lateinfra-earlyinfra')
    hold on 
    
    subplot(3,4,2)
    y = mean(lateinfra_gran(cond,:,:),3);
    sem = std(lateinfra_gran(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
    title('lateinfra-gran')
    hold on 
    
    subplot(3,4,3)
    y = mean(lateinfra_supra(cond,:,:),3);
    sem = std(lateinfra_supra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
    title('lateinfra-supra')
    hold on 
    
        subplot(3,4,4)
    y = mean(earlyinfra_supra(cond,:,:),3);
    sem = std(earlyinfra_supra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
     title('earlyinfra-supra')
    hold on 
    
        subplot(3,4,5)
    y = mean(earlyinfra_gran(cond,:,:),3);
    sem = std(earlyinfra_gran(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
     title('earlyinfra-gran')
    hold on 
    
        subplot(3,4,6)
    y = mean(earlyinfra_lateinfra(cond,:,:),3);
    sem = std(earlyinfra_lateinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
         title('earlyinfra-lateinfra')

    hold on 
        subplot(3,4,7)
    y = mean(supra_earlyinfra(cond,:,:),3);
    sem = std(supra_earlyinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
             title('supra-earlyinfra')
    hold on 
        subplot(3,4,8)
    y = mean(supra_lateinfra(cond,:,:),3);
    sem = std(supra_lateinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
                 title('supra-lateinfra')

    hold on 
        subplot(3,4,9)
    y = mean(supra_granular(cond,:,:),3);
    sem = std(supra_granular(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
                     title('supra-granular')

    hold on 
        subplot(3,4,10)
    y = mean(gran_earlyinfra(cond,:,:),3);
    sem = std(gran_earlyinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
     title('gran-earlyinfra')

    hold on 
        subplot(3,4,11)
    y = mean(gran_lateinfra(cond,:,:),3);
    sem = std(gran_lateinfra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
     title('gran-lateinfra')
    hold on 
        subplot(3,4,12)
    y = mean(gran_supra(cond,:,:),3);
    sem = std(gran_supra(cond,:,:),0,3)/sqrt(9);
    shadedErrorBar(x,y,sem,colors{cond},1);
    title('gran-supra')
    hold on 
    
end


%%
G_matrix = {}
for cond = 1:size(F,1)
    for anim = 1:size(F{1},2)
        G_matrix{cond}(:,:,anim)= F{cond}{anim} ;
    end
end
%save('partial_granger_singletrial','G_matrix');
logG{1}= log(G_matrix{1});
logG{2}= log(G_matrix{2});
[h,p,ci,stats]= ttest(logG{1},logG{2},'dim',3);
[hh, crit_p, adj_ci_cvrg, adj_p]= fdr_bh(p,0.05,'pdep','yes');
%figure()
%%
G_matrix_prism = cell(1,2);
for cond = 1 :size(G_matrix,2)
    for anim =1 :size(G_matrix{cond},3)
        aa = G_matrix{cond}(:,:,anim);
        aa = aa(~isnan(aa));
        %G_matrix_prism{cond}(:,anim)= log(aa);
        G_matrix_prism{cond}(:,anim)= aa;
    end
end

save('G_matrix_prism_23Nov2020','G_matrix_prism')
%%
condname = {'Conditional pairwise GC before Hyase','Conditional pairwise GC after Hyase'};
netowrkname = {'G-causal network graph before Hyase','G-causal network graph after Hyase'};
for cond = 1:size(F,1)
    count = 0
    median_F{cond} = median(G_matrix{cond},3);
    a = reshape(median_F{cond},1,16); a = a(~isnan(a));
    b = reshape(p,1,16); b = b(~isnan(b));
    indb = find(b<=crit_p); %a(indb) = 0.0000001 ; % replace the nonsignficant values of  granger values with zero to plot in the network graph
    s = [1 1 1 2 2 2 3 3 3 4 4 4];
    t = [2 3 4 1 3 4 1 2 4 1 2 3];
    x = [2 1 1.5 2]; y = [1 2 3 4]; z = [0 0 0 0];
    pair{1} = [2,3,4]; pair{2}=[1,3,4];pair{3}= [1,2,4]; pair{4}= [1,2,3];
    weight{1}= a(1:3);weight{2}= a(4:6);weight{3}= a(7:9);weight{4}= a(10:12);
    names = {'Layer Vb','Layer VI','Layer III/IV','Layer I/II'};
    %k(cond)= subplot(3,2,cond)
    fig1 = figure()
    
    G = digraph(s,t,a,names);
    h = plot(G,'XData',x,'YData',y,'ZData',z,'NodeFontSize',16,'MarkerSize',10,'NodeColor',[0 0 0]);view(2)
    set(gca,'visible','off')
    for i = 1 :4
        maxweight(i) = max(weight{i});
    end
    maxweigtall = max(maxweight);
    for i = 1 :4
        weight{i} = weight{i}./maxweigtall;
    end
    for i = 1:4
        for k = 1:3
            s1 = i; s2 = pair{i}(k); w = weight{i}(k)
            count = count+1
            if ismember(count,indb)
                highlight(h,[s1 s2],'EdgeColor','r','LineWidth',12*w);
            else
                highlight(h,[s1 s2],'EdgeColor','k','LineWidth',12*w);
            end
        end
    end
    annotation(fig1,'textbox',...
        [0.196969696969697 0.834791059280854 0.633608815426998 0.0505344995140912],...
        'String',netowrkname{cond},...
        'FontWeight','bold',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    %k(cond+2)=subplot(3,2,cond+2)
    figure()
    x = repmat(1:size(p,2),size(p,1),1); % generate x-coordinates
    y = x'; % generate y-coordinates
    updatedF = round(median_F{cond},4)
    t = num2cell(updatedF); % extact values into cells
    t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
    colormap(cm);
    imagesc(updatedF);
    caxis([0,0.2]);
    title(condname{cond});
    colorbar()
    xlabel('From','FontSize',10);ylabel('To','FontSize',10);
    set(gca,'Ytick',1:4,'YTickLabel',names,'FontSize',10)
    set(gca,'Xtick',1:4,'XTickLabel',names,'FontSize',10)
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center','Color','r')
end
%k(5) = subplot(3,2,5); % the last (odd) axes
%set(k(5),'Position',[0.3, 0.1, 0.3, 0.2])
% Generate Labels
figure()
t = num2cell(p); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
colormap(cm);
imagesc(p); colorbar();
title(sprintf('The critical p-value is %f',crit_p));
xlabel('From','FontSize',10);ylabel('To','FontSize',10);
set(gca,'Ytick',1:4,'YTickLabel',names,'FontSize',10)
set(gca,'Xtick',1:4,'XTickLabel',names,'FontSize',10)
text(x(:), y(:), t, 'HorizontalAlignment', 'Center','Color','r')