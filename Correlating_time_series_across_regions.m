clc;
clear all;
% close all;
%% Settings
conditions=1:5;
categories=1:3;
n_trials_per_cond=3;

Similar_conditions=1; % 1=yes; 2=no;
Which_conditions=[1 2];

applying_ERPs=0;
ERP_type=2;

SNR=0.99;

Regions=[1 3]; % regions are in order {'LPS','LPI','LAS','LAI','RPS','RPI','RAS','RAI'}


load_simulated_dataset=1;

directional=0; % should the data have a direction or not?

%% Finding the right seeds that result in directed flow of signals
if load_simulated_dataset==0
    counter=zeros(1,max(conditions));
    for cond=conditions
        while 1
            t=rng;
            [~,Sig,error_flag]=gen_ar_biv_for_finding_seeds(300,5,t,directional);
            if directional
                direction(conditions)=[Sig(1,1)~=1 && Sig(2,1)==1 && Sig(1,2)~=1 && Sig(2,2)~=1]; % preferred direction settings
            else
                direction(conditions)=[Sig(2,1)~=1 && Sig(1,2)~=1]; % no interaction bw sources
            end
            if direction(cond)
                s(cond)=t.Seed;
                rng('shuffle');
                break;
            end
            rng('shuffle');
        end
        counter(cond)=counter(cond)+1;
    end
    %% Generating datasets and applying ERPs
    for cond=conditions
        for cate=categories
            dataset_string = ['RDM_benchmark_cond_',num2str(cond),'_cate_',num2str(cate)];
            if Similar_conditions==1 && sum(cond==Which_conditions)
                generate_datasets_ar(n_trials_per_cond,dataset_string,s(min(Which_conditions)),SNR,Regions);
            else
                generate_datasets_ar(n_trials_per_cond,dataset_string,s(cond),SNR,Regions);
            end
            [cond cate]
        end
    end
    ccc
    x=[1:200];
    
    % A=[0.1:0.1:1];
    A=ones(1,10);
    % d=[10:20];
    d=[10:10:100];
    
    % Combinations_params=randperm([1:max(conditions)],[2 max(conditions)]);
    Combinations_params=[1:length(conditions);1:length(conditions)];
    
    
    % A=0.1;
    % d=10;
    s=4;
    
    
    LB=-1;
    UB=5;
    N=40;
    
    counter_cond_cate=0;
    for cond=conditions
        for cate=categories
            counter_cond_cate=counter_cond_cate+1;
            
            dataset_string = ['RDM_benchmark_cond_',num2str(cond),'_cate_',num2str(cate)];
            for trial = 1:n_trials_per_cond
                load(['data/' dataset_string '/MEG/dataset_' num2str(trial),'/data.mat'],'MEG_data');
                
                if ERP_type==1
                    % Gaussian
                    ERP=exp(-(((x-d(Combinations_params(2,cond))).^2)./(2*(s.^2))));
                    ERP=A(Combinations_params(1,cond)).*[zeros(1,d(Combinations_params(2,cond))) ERP];
                else
                    %%mexiacan hat
                    [ERP,~] = mexihat(LB,UB,N);
                    ERP=A(Combinations_params(1,cond)).*[zeros(1,d(Combinations_params(2,cond))) ERP];
                end
                
                for ch=1:size(MEG_data,1)
                    
                    if applying_ERPs
                        tmp=conv(ERP,MEG_data(ch,:));
                    else
                        tmp=MEG_data(ch,:);
                    end
                    
                    MEG_datas(ch,1:length(tmp),counter_cond_cate,trial)=tmp;
                end
            end
        end
    end
    save('simulated_dataset.mat','MEG_datas');
else
    load('simulated_dataset.mat');
end

%% Channel selection
load('C:\Users\mq20185770\Documents\MATLAB\Postdoc\Info_Connectivity\BBCB_code\data\sa.mat')
frontal_char='LF'; % left frontal
occipital_char='LO'; % left occipital

counter_fr=0;
counter_oc=0;
for ch=1:length(sa.MEG_clab_electrodes)
    if strncmp(sa.MEG_clab_electrodes{1,ch}(2:3),frontal_char,2)
        counter_fr=counter_fr+1;
        channel_indx_frontal(counter_fr)=ch;
    end
    if strncmp(sa.MEG_clab_electrodes{1,ch}(2:3),occipital_char,2)
        counter_oc=counter_oc+1;
        channel_indx_occipital(counter_oc)=ch;
    end
    
end

%% RSA analysis
COR_total=nan*ones(size(MEG_datas,3),size(MEG_datas,3),size(MEG_datas,2),3);
for area=1:3 % 1=frontal, 2=occipital, 3=whole head
    if area==1
        channels=channel_indx_frontal;
    elseif area==2
        channels=channel_indx_occipital;
    elseif area==3
        channels=[1:length(sa.MEG_clab_electrodes)];
    end
    
    cors_t=nan*ones(size(MEG_datas,3),size(MEG_datas,3),size(MEG_datas,2),size(MEG_datas,3));
    for time=1:size(MEG_datas,2)
        for cond1=1:size(MEG_datas,3)
            for cond2=cond1:size(MEG_datas,3)
                trials=0;
                for trial1=1:size(MEG_datas,4)
                    if cond1==cond2
                        span=1:size(MEG_datas,4);
                    else
                        span=trial1:size(MEG_datas,4);
                    end
                    
                    for trial2=span
                        trials=trials+1;
                        cors_t(cond1,cond2,time,trials)=corr(squeeze(nanmean(MEG_datas(channels,time,cond1,trial1),2)),squeeze(nanmean(MEG_datas(channels,time,cond2,trial2),2)));
                    end
                end
            end
        end
        COR_total(:,:,time,area)=nanmean(cors_t(:,:,time,:),4);
    end
        figure;
        imagesc(nanmean(COR_total(:,:,:,area),3),[-1 1]);
        colorbar;
end

for time=1:size(MEG_datas,2)
    for area=1:3
        for cond1=1:size(MEG_datas,3)
            COR_total(cond1,cond1,:,area)=nan;
        end
    end
end
%% Generating model RDMs
Connected_conditions=[1 1;2 2;3 3; 4 4;5 5];
Model_RDM_tmp=nan*ones(max(conditions)*max(categories),size(Connected_conditions,1),size(Connected_conditions,1));

for row=1:size(Connected_conditions,1)
    for cond1=1:max(conditions)*max(categories)
        for cond2=cond1+1:max(conditions)*max(categories)
            if sum(cond1==[((Connected_conditions(row,1)-1)*max(categories)+1):(Connected_conditions(row,1)*max(categories))])>0 && sum(cond2==[((Connected_conditions(row,2)-1)*max(categories)+1):(Connected_conditions(row,2)*max(categories))])>0
                Model_RDM_tmp(cond1,cond2,row)=1;
            else
                Model_RDM_tmp(cond1,cond2,row)=0;
            end
        end
    end
    %     figure;
    %     imagesc(squeeze(Model_RDM_tmp(:,:,row)))
end
Model_RDM_conditions=nansum(Model_RDM_tmp,3);
figure;
imagesc(Model_RDM_conditions,[-1 1])



Connected_conditions=[1 1;1 2;2 2;3 3;3 4;3 5;4 5;4 4;5 5];
Model_RDM_tmp=nan*ones(max(conditions)*max(categories),size(Connected_conditions,1),size(Connected_conditions,1));

for row=1:size(Connected_conditions,1)
    for cond1=1:max(conditions)*max(categories)
        for cond2=cond1+1:max(conditions)*max(categories)
            if sum(cond1==[((Connected_conditions(row,1)-1)*max(categories)+1):(Connected_conditions(row,1)*max(categories))])>0 && sum(cond2==[((Connected_conditions(row,2)-1)*max(categories)+1):(Connected_conditions(row,2)*max(categories))])>0
                Model_RDM_tmp(cond1,cond2,row)=1;
            else
                Model_RDM_tmp(cond1,cond2,row)=0;
            end
        end
    end
    %     figure;
    %     imagesc(squeeze(Model_RDM_tmp(:,:,row)))
end
Model_RDM_task=nansum(Model_RDM_tmp,3);
figure;
imagesc(Model_RDM_task,[-1 1])


%% Using Multi-variate Granger causality analysis using MVGC toolbox (using signal time series)

% ([84 88 92 97 102 51 60 65 70 83])  % occipital and frontal electrodes

% mvgc_Hamid(squeeze(nanmean(MEG_datas([channel_indx_occipital channel_indx_frontal],:,:,1),4)))

%% Using Multi-variate Granger causality analysis using Clark's method

% past_smoothing=5; % *5ms
sizem=size(MEG_datas,3);

for t=1:size(COR_total,3)
    Front_task (t)= corr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),reshape(Model_RDM_task,[sizem*sizem 1]),'row','complete');
    Front_conditions (t)= corr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),'row','complete');
    
    Back_task (t)= corr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),reshape(Model_RDM_task,[sizem*sizem 1]),'row','complete');
    Back_conditions (t)= corr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),'row','complete');
end
% plot([Front_task;Front_conditions;Back_task;Back_conditions]')

% mvgc_Hamid([Front_task;Front_conditions;Back_task;Back_conditions])

%% Erin's
delay=20; % *5ms
past_smoothing=5; % *5ms
sizem=size(MEG_datas,3);

for t=delay+1:size(COR_total,3)
    FF (t)= partialcorr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),'row','complete');
    FB (t)= partialcorr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),'row','complete');
end


% significance testing
for t=delay+1:size(COR_total,3)
    for iteration=1:1000
        FF_iter (t,iteration)= partialcorr(randsample(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),sizem*sizem),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),'row','complete');
        FB_iter (t,iteration)= partialcorr(randsample(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),sizem*sizem),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),'row','complete');
        FF_minus_FB_iter (t,iteration)=FF_iter (t,iteration)-FB_iter (t,iteration);
    end
    [t]
end


significance_threshold=0.05;


figure;
subplot(3,1,1)
plot(FF)
hold on;
for t=delay+1:size(COR_total,3)
    if FF(t)>0
        FF_sign(t)=1-sum(FF(t)>FF_iter(t,:))./size(FF_iter,2);
    else
        FF_sign(t)=1-sum(FF(t)<FF_iter(t,:))./size(FF_iter,2);
    end
end
FF_sign(1:delay)=1;
try
FF_sign(delay+1:end)=mafdr(FF_sign(delay+1:end));
end
plot([1:size(FF_iter,1)],(FF_sign<significance_threshold).*(sign(FF)),'*')

subplot(3,1,2)
plot(FB)
hold on;
for t=delay+1:size(COR_total,3)
    if FB(t)>0
        FB_sign(t)=1-sum(FB(t)>FB_iter(t,:))./size(FB_iter,2);
    else
        FB_sign(t)=1-sum(FB(t)<FB_iter(t,:))./size(FB_iter,2);
    end
end
FB_sign(1:delay)=1;
try
FB_sign(delay+1:end)=mafdr(FB_sign(delay+1:end));
end
plot([1:size(FB_iter,1)],(FB_sign<significance_threshold).*(sign(FB)),'*')



subplot(3,1,3)
plot(FF-FB)
hold on;
for t=delay+1:size(COR_total,3)
    if (FF(t)-FB(t))>0
        FF_FB_sign(t)=1-sum((FF(t)-FB(t))>(FF_iter(t,:)-FB_iter(t,:)))./size(FF_iter,2);
    else
        FF_FB_sign(t)=1-sum((FF(t)-FB(t))<(FF_iter(t,:)-FB_iter(t,:)))./size(FF_iter,2);
    end
end
FF_FB_sign(1:delay)=1;
try
FF_FB_sign(delay+1:end)=mafdr(FF_FB_sign(delay+1:end));
end

plot([1:size(FF_iter,1)],(FF_FB_sign<significance_threshold).*(sign(FF-FB)),'*')

%% Hamid's

delay=20; % *5ms
sizem=size(MEG_datas,3);

for t=delay+1:size(COR_total,3)
    
    Front(t)= corr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),'row','complete');
    Front_minus_back(t)= partialcorr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),'row','complete');
    FF(t)=Front(t)-Front_minus_back(t);
    
    Back(t)= corr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),'row','complete');
    Back_minus_front(t)= partialcorr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),reshape(Model_RDM_conditions,[sizem*sizem 1]),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),'row','complete');
    FB(t)=Back(t)-Back_minus_front(t);
    
end

% Significance testing
for t=delay+1:size(COR_total,3)
    for iteration=1:1000
        
        Front_iter= corr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),randsample(reshape(Model_RDM_conditions,[sizem*sizem 1]),sizem*sizem),'row','complete');
        Front_minus_back_iter= partialcorr(reshape(COR_total(:,:,t,1),[sizem*sizem 1]),randsample(reshape(Model_RDM_conditions,[sizem*sizem 1]),sizem*sizem),reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]),'row','complete');
        FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
        
        Back_iter= corr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),randsample(reshape(Model_RDM_conditions,[sizem*sizem 1]),sizem*sizem),'row','complete');
        Back_minus_front_iter= partialcorr(reshape(COR_total(:,:,t,2),[sizem*sizem 1]),randsample(reshape(Model_RDM_conditions,[sizem*sizem 1]),sizem*sizem),reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]),'row','complete');
        FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
        
        FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
    end
    [t]
end


significance_threshold=0.95;


figure;
subplot(3,1,1)
plot(FF)
hold on;
for t=delay+1:size(COR_total,3)
    if FF(t)>0
        FF_sign(t)=1-sum(FF(t)>FF_iter(t,:))./size(FF_iter,2);
    else
        FF_sign(t)=1-sum(FF(t)<FF_iter(t,:))./size(FF_iter,2);
    end
end
FF_sign(1:delay)=1;
try
FF_sign(delay+1:end)=mafdr(FF_sign(delay+1:end));
end
plot([1:size(FF_iter,1)],(FF_sign<significance_threshold).*(sign(FF)),'*')

subplot(3,1,2)
plot(FB)
hold on;
for t=delay+1:size(COR_total,3)
    if FB(t)>0
        FB_sign(t)=1-sum(FB(t)>FB_iter(t,:))./size(FB_iter,2);
    else
        FB_sign(t)=1-sum(FB(t)<FB_iter(t,:))./size(FB_iter,2);
    end
end
FB_sign(1:delay)=1;
try
FB_sign(delay+1:end)=mafdr(FB_sign(delay+1:end));
end
plot([1:size(FB_iter,1)],(FB_sign<significance_threshold).*(sign(FB)),'*')



subplot(3,1,3)
plot(FF-FB)
hold on;
for t=delay+1:size(COR_total,3)
    if (FF(t)-FB(t))>0
        FF_FB_sign(t)=1-sum((FF(t)-FB(t))>(FF_iter(t,:)-FB_iter(t,:)))./size(FF_iter,2);
    else
        FF_FB_sign(t)=1-sum((FF(t)-FB(t))<(FF_iter(t,:)-FB_iter(t,:)))./size(FF_iter,2);
    end
end
FF_FB_sign(1:delay)=1;
FF_FB_sign(delay+1:end)=mafdr(FF_FB_sign(delay+1:end));

plot([1:size(FF_iter,1)],(FF_FB_sign<significance_threshold).*(sign(FF-FB)),'*')



%% Tim's unexplained variance
delay=20;

for t=delay+1:size(COR_total,3)
    
    % FF   
    
    RDM_targ_past=reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]);
    RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
    
    RDM_targ_prsnt=reshape(COR_total(:,:,t,1),[sizem*sizem 1]);
    RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
    
    RDM_srce_past=reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]);
    RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
    
    
    % Linear regression
    [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
    [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
    FF_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
    
    % General Linear Model
    
    T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
    T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
    glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
    glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
    
    FF_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
    
    % FB
    
    RDM_targ_past=reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]);
    RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
    
    RDM_targ_prsnt=reshape(COR_total(:,:,t,2),[sizem*sizem 1]);
    RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
    
    RDM_srce_past=reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]);
    RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
    
    
    % Linear regression
    [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
    [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
    FB_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
    
    % General Linear Model
    
    T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
    T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
    glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
    glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
    
    FB_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
end


for t=delay+1:size(COR_total,3)
    for iteration=1:1000
           
        % FF
        
        RDM_targ_past=reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]);
        RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
        
        RDM_targ_prsnt=reshape(COR_total(:,:,t,1),[sizem*sizem 1]);
        RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
        RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));
        
        
        RDM_srce_past=reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]);
        RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
        
        
        % Linear regression
        [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
        [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
        FF_LR_iter (t,iteration)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
        
%         % General Linear Model
%         
%         T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%         T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%         glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         
%         FF_GLM_iter (t,iteration)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
        
        % FB
        
        RDM_targ_past=reshape(COR_total(:,:,t-delay,2),[sizem*sizem 1]);
        RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
        
        RDM_targ_prsnt=reshape(COR_total(:,:,t,2),[sizem*sizem 1]);
        RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
        RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));

        RDM_srce_past=reshape(COR_total(:,:,t-delay,1),[sizem*sizem 1]);
        RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
        
        
        % Linear regression
        [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
        [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
        FB_LR_iter (t,iteration)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
        
%         % General Linear Model
        
%         T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%         T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%         glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         
%         FB_GLM_iter (t,iteration)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
    end
    [t]
end


figure;
subplot(3,1,1)
plot(FF_LR)
hold on;
for t=delay+1:size(COR_total,3)
    if FF_LR(t)>0
        FF_LR_sign(t)=1-sum(FF_LR(t)>FF_LR_iter(t,:))./size(FF_LR_iter,2);
    else
        FF_LR_sign(t)=1-sum(FF_LR(t)<FF_LR_iter(t,:))./size(FF_LR_iter,2);
    end
end
FF_LR_sign(1:delay)=1;
FF_LR_sign(delay+1:end)=mafdr(FF_LR_sign(delay+1:end));
plot([1:size(FF_iter,1)],(FF_LR_sign<significance_threshold).*(sign(FF_LR))*5,'*')

subplot(3,1,2)
plot(FB_LR)
hold on;
for t=delay+1:size(COR_total,3)
    if FB_LR(t)>0
        FB_LR_sign(t)=1-sum(FB_LR(t)>FB_LR_iter(t,:))./size(FB_LR_iter,2);
    else
        FB_LR_sign(t)=1-sum(FB_LR(t)<FB_LR_iter(t,:))./size(FB_LR_iter,2);
    end
end
FB_LR_sign(1:delay)=1;
FB_LR_sign(delay+1:end)=mafdr(FB_LR_sign(delay+1:end));
plot([1:size(FB_LR_iter,1)],(FB_LR_sign<significance_threshold).*(sign(FB_LR))*5,'*')



subplot(3,1,3)
plot(FF_LR-FB_LR)
hold on;
for t=delay+1:size(COR_total,3)
    if (FF_LR(t)-FB_LR(t))>0
        FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))>(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
    else
        FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))<(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
    end
end
FF_FB_sign_LR(1:delay)=1;
FF_FB_sign_LR(delay+1:end)=mafdr(FF_FB_sign_LR(delay+1:end));

plot([1:size(FF_iter,1)],(FF_FB_sign_LR<significance_threshold).*(sign(FF_LR-FB_LR))*5,'*')

