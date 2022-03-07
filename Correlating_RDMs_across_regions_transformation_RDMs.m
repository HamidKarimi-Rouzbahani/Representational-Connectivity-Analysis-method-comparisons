clc;
% clear all;
close all;

repititions=1;
for repitition=1:repititions
    
    A=randn(10,16,100); % channel_in_region * condition_in_region * time sample
    RDM_source=nan*ones((size(A,2)*size(A,2)-size(A,2))./2,size(A,3),2);
    
    for numb=1:2
        
        cors_source=nan*ones(size(A,2),size(A,2));
        
        for time=1:size(A,3)
            for i=1:size(A,2)
                for j=i+1:size(A,2)
                    cors_source(i,j)=corr(A(:,i,time),A(:,j,time));
                end
            end
            RDM_source(:,time,numb)=cors_source(~isnan(cors_source));
            cors_save(:,:,time,numb)=cors_source;
        end
        %     imagesc(cors_source);
        %     imagesc(squeeze(nanmean(cors_save(:,:,:,1),3)))
        
        X=eye(size(A,2));
        
        channels=[1:size(A,2)];
        
        
        if numb==1 % occipital
            
            Identity={'A','B','C','D','E','F','G','H','A','B','C','D','E','F','G','H'};
            %         Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
            %         Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
            
            for i=1:size(A,2)
                for j=1:size(A,2)
                    if strcmp(Identity{i},Identity{j})==1
                        X(i,j)=0.3+randn*0.1;
                    else
                        X(i,j)=randn*0.1;
                    end
                end
            end
            
        else % frontal
            
            %         Identity={'A','B','C','D','E','F','G','H','A','B','C','D','E','F','G','H'};
            Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
            Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
            
            for i=1:size(A,2)
                for j=1:size(A,2)
                    if sum([strcmp(Semantic{i},Semantic{j})+strcmp(Task{i},Task{j})])==2
                        X(i,j)=0.3+randn*0.1;
                    elseif sum([strcmp(Semantic{i},Semantic{j})+strcmp(Task{i},Task{j})])==1
                        X(i,j)=0.2+randn*0.1;
                    else
                        X(i,j)=randn*0.1;
                    end
                end
            end
            
        end
        
        cors_source=nan*ones(size(A,2),size(A,2));
        for time=1:size(A,3)
            Y=squeeze(A(:,:,time))*X;
            for i=1:size(A,2)
                for j=i+1:size(A,2)
                    cors_source(i,j)=corr(Y(:,i),Y(:,j));
                end
            end
            RDM_source(:,time,numb)=cors_source(~isnan(cors_source));
            cors_save(:,:,time,numb)=cors_source;
        end
        %     for numb=1:2
        
        subplot(1,2,numb)
        %     imagesc(cors_source,[-1 1]);
        % figure
%         imagesc(squeeze(nanmean(cors_save(:,:,:,numb),3)))
    end
    
    %% Correlating RDMs in time
    clc;
    direction=2; % 1= (1=>2) FF; 2= (2=>1) FB; 3= (1<=>2) Bidirectional
    
    
    N=size(A,3); % samples size
    P=5;  % model order
    M = size(RDM_source,1); %number of cells in RDM;
    
    sigma = 0.1; %scale of random AR-parameters
    
    N0=ceil(N./20); %length of ignored start
    lambdamax=10;
    
    c=0;
    while (lambdamax>1 || lambdamax< 0.9)
        c=c+1;
        Arsig=[];
        for k=1:P
            aloc = zeros(M*2);
            for i=1:M*2
                
                aloc(i,i)=abs(randn)*sigma;  %everything related to its past
                if direction==1
                    if i>M
                        aloc(i,i-M)=abs(randn)*sigma;  % feed-forward
                    end
                elseif direction==2
                    if i<=M
                        aloc(i,i+M)=abs(randn)*sigma;  % feedback
                    end
                end
            end
            Arsig=[Arsig,aloc];
        end
        E=eye(M*2*P);
        AA=[Arsig;E(1:end-M*2,:)];
        lambda=eig(AA);
        lambdamax=max(abs(lambda));
        %         [c lambdamax]
    end
    
    x=[[squeeze(RDM_source(:,:,1));squeeze(RDM_source(:,:,2))] zeros(M*2,N0)];
    y=x;
    for i=P+1:N+N0
        yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
        y(:,i)=Arsig*yloc+x(:,i);
        sums(i)=sum(y(1:M,i)-y(M+1:end,i));
    end
    data=y(:,N0+1:end);
    
    %% Evaluating if the correlated RDMs are correlated using simple time-shifted correlation
    
    delays=[5];
    
    cors_FF=nan*ones(size(data,2),length(delays));
    cors_FB=nan*ones(size(data,2),length(delays));
    cors3=nan*ones(size(data,2),length(delays));
    cors4=nan*ones(size(data,2),length(delays));
    
    d=0;
    for delay=delays
        d=d+1;
        for time=1+delay:size(data,2)-delay
            %         cors_FF(time,d)=corr(data(M+1:end,time),nanmean(data(1:M,time-delay:time),2),'Type','Spearman'); % occipital-frontal : feedforward
            %         cors_FB(time,d)=corr(data(1:M,time),nanmean(data(M+1:end,time-delay:time),2),'Type','Spearman'); % frontal-occipital : feedback
            cors_FF(time,d)=corr(data(M+1:end,time),nanmean(data(1:M,time-delay:time),2)); % occipital-frontal : feedforward
            cors_FB(time,d)=corr(data(1:M,time),nanmean(data(M+1:end,time-delay:time),2)); % frontal-occipital : feedback
        end
    end
    figure;
    plot(cors_FF)
    hold on;
    plot(cors_FB)
    
    %% Generating model RDMs
    % conditions=1:4;
    % categories=1:4;
    % Connected_conditions=[1 1;2 2;3 3; 4 4];
    %
    % Model_RDM_tmp=nan*ones(max(conditions)*max(categories),size(Connected_conditions,1),size(Connected_conditions,1));
    %
    % for row=1:size(Connected_conditions,1)
    %     for cond1=1:max(conditions)*max(categories)
    %         for cond2=cond1+1:max(conditions)*max(categories)
    %             if sum(cond1==[((Connected_conditions(row,1)-1)*max(categories)+1):(Connected_conditions(row,1)*max(categories))])>0 && sum(cond2==[((Connected_conditions(row,2)-1)*max(categories)+1):(Connected_conditions(row,2)*max(categories))])>0
    %                 Model_RDM_tmp(cond1,cond2,row)=1;
    %             else
    %                 Model_RDM_tmp(cond1,cond2,row)=0;
    %             end
    %         end
    %     end
    % end
    % Model_RDM_identities=nansum(Model_RDM_tmp,3);
    % % figure;
    % % imagesc(Model_RDM_identities,[-1 1])
    %
    %
    %
    % Connected_conditions=[1 1;1 2;2 2;3 3;3 4;4 4];
    % Model_RDM_tmp=nan*ones(max(conditions)*max(categories),size(Connected_conditions,1),size(Connected_conditions,1));
    %
    % for row=1:size(Connected_conditions,1)
    %     for cond1=1:max(conditions)*max(categories)
    %         for cond2=cond1+1:max(conditions)*max(categories)
    %             if sum(cond1==[((Connected_conditions(row,1)-1)*max(categories)+1):(Connected_conditions(row,1)*max(categories))])>0 && sum(cond2==[((Connected_conditions(row,2)-1)*max(categories)+1):(Connected_conditions(row,2)*max(categories))])>0
    %                 Model_RDM_tmp(cond1,cond2,row)=1;
    %             else
    %                 Model_RDM_tmp(cond1,cond2,row)=0;
    %             end
    %         end
    %     end
    % end
    % Model_RDM_task_semantics=nansum(Model_RDM_tmp,3);
    % % figure;
    % % imagesc(Model_RDM_task_semantics,[-1 1])
    % c=0;
    % for cond1=1:max(conditions)*max(categories)
    %     for cond2=cond1+1:max(conditions)*max(categories)
    %         c=c+1;
    %         Model_RDM_task_tmp(c,1)=Model_RDM_task_semantics(cond1,cond2);
    %         Model_RDM_conditions_tmp(c,1)=Model_RDM_identities(cond1,cond2);
    %     end
    % end
    % Model_RDM_task_semantics=Model_RDM_task_tmp;
    % Model_RDM_identities=Model_RDM_conditions_tmp;
    %
    % ccc
    
    Model_RDM_identities=nan*ones(size(A,2));
    Model_RDM_task_semantics=nan*ones(size(A,2));
    
    Identity={'A','B','C','D','E','F','G','H','A','B','C','D','E','F','G','H'};
    %         Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
    %         Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
    
    for i=1:size(A,2)
        for j=i+1:size(A,2)
            if strcmp(Identity{i},Identity{j})==1
                Model_RDM_identities(i,j)=1;
            else
                Model_RDM_identities(i,j)=0;
            end
        end
    end
%     figure;
%     imagesc(Model_RDM_identities)
    Model_RDM_identities=Model_RDM_identities(~isnan(Model_RDM_identities));
    
    %         Identity={'A','B','C','D','E','F','G','H','A','B','C','D','E','F','G','H'};
    Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
    Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
    
    for i=1:size(A,2)
        for j=i+1:size(A,2)
            if sum([strcmp(Semantic{i},Semantic{j})+strcmp(Task{i},Task{j})])==2
                Model_RDM_task_semantics(i,j)=1;
            elseif sum([strcmp(Semantic{i},Semantic{j})+strcmp(Task{i},Task{j})])==1
                Model_RDM_task_semantics(i,j)=0.5;
            else
                Model_RDM_task_semantics(i,j)=0;
            end
        end
    end
%     figure;
%     imagesc(Model_RDM_task_semantics)
    Model_RDM_task_semantics=Model_RDM_task_semantics(~isnan(Model_RDM_task_semantics));
    
    
    %% Using Multi-variate Granger causality analysis using Clark's method
    
    % past_smoothing=5; % *5ms
    
    for t=1:size(data,2)
        Front_task_semantics (t)= corr(data(M+1:end,t),Model_RDM_task_semantics,'row','complete');
        Front_identities (t)= corr(data(M+1:end,t),Model_RDM_identities,'row','complete');
        
        Back_task_semantics (t)= corr(data(1:M,t),Model_RDM_task_semantics,'row','complete');
        Back_identities (t)= corr(data(1:M,t),Model_RDM_identities,'row','complete');
    end
    % figure;
    % plot([Front_task;Front_conditions;Back_task;Back_conditions]')
    significance_level=0.01;
    MVGC(:,:,repitition)=mvgc_Hamid([Front_task_semantics;Front_identities;Back_task_semantics;Back_identities],significance_level);
%     close all
%     ccc
    %% Erin's
    delay=5; % *5ms
    
    for t=delay+1:size(data,2)
        FF (t)= partialcorr(data(M+1:end,t),nanmean(data(1:M,t-delay:t),2),nanmean(data(M+1:end,t-delay:t),2),'row','complete');
        FB (t)= partialcorr(data(1:M,t),nanmean(data(M+1:end,t-delay:t),2),nanmean(data(1:M,t-delay:t),2),'row','complete');
    end
    
    
    % significance testing
    for t=delay+1:size(data,2)
        for iteration=1:1000
            FF_iter (t,iteration)= partialcorr(randsample(data(M+1:end,t),M),nanmean(data(1:M,t-delay:t),2),nanmean(data(M+1:end,t-delay:t),2),'row','complete');
            FB_iter (t,iteration)= partialcorr(randsample(data(1:M,t),M),nanmean(data(M+1:end,t-delay:t),2),nanmean(data(1:M,t-delay:t),2),'row','complete');
            FF_minus_FB_iter (t,iteration)=FF_iter (t,iteration)-FB_iter (t,iteration);
        end
%         [t]
    end
    
    
    significance_threshold=0.05;
    
    
    figure;
    subplot(3,1,1)
    plot(FF)
    hold on;
    for t=delay+1:size(data,2)
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
    for t=delay+1:size(data,2)
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
    for t=delay+1:size(data,2)
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
    
    FFs(:,repitition)=FF;
    FBs(:,repitition)=FB;
    FFs_iter(:,:,repitition)=FF_iter;
    FBs_iter(:,:,repitition)=FB_iter;
    
    close all;
    %% Hamid's
    
    delay=5; % *5ms
    
    chosen_model=Model_RDM_identities;
    % chosen_model=Model_RDM_task_semantics;
    
    for t=delay+1:size(data,2)
        
        Front(t)= corr(data(M+1:end,t),chosen_model,'row','complete');
        Front_minus_back(t)= partialcorr(data(M+1:end,t),chosen_model,nanmean(data(1:M,t-delay:t),2),'row','complete');
        FF(t)=Front(t)-Front_minus_back(t);
        
        Back(t)= corr(data(1:M,t),chosen_model,'row','complete');
        Back_minus_front(t)= partialcorr(data(1:M,t),chosen_model,nanmean(data(M+1:end,t-delay:t),2),'row','complete');
        FB(t)=Back(t)-Back_minus_front(t);
        
    end
    
    % Significance testing
    for t=delay+1:size(data,2)
        for iteration=1:1000
            
            Front_iter= corr(data(M+1:end,t),randsample(chosen_model,M),'row','complete');
            Front_minus_back_iter= partialcorr(data(M+1:end,t),randsample(chosen_model,M),nanmean(data(1:M,t-delay:t),2),'row','complete');
            FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
            
            Back_iter= corr(data(1:M,t),randsample(chosen_model,M),'row','complete');
            Back_minus_front_iter= partialcorr(data(1:M,t),randsample(chosen_model,M),nanmean(data(M+1:end,t-delay:t),2),'row','complete');
            FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
            
            FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
        end
%         [t]
    end
    
    
    significance_threshold=0.05;
    
    
    figure;
    subplot(3,1,1)
    plot(FF)
    hold on;
    for t=delay+1:size(data,2)
        if FF(t)>0
            FF_sign(t)=1-sum(FF(t)>FF_iter(t,:))./size(FF_iter,2);
        else
            FF_sign(t)=1-sum(FF(t)<FF_iter(t,:))./size(FF_iter,2);
        end
    end
    % FF_sign(1:delay)=1;
    try
        FF_sign(delay+1:end)=mafdr(FF_sign(delay+1:end));
    end
    plot([1:size(FF_iter,1)],(FF_sign<significance_threshold).*(sign(FF)),'*')
    ylim([-1 1])
    
    subplot(3,1,2)
    plot(FB)
    hold on;
    for t=delay+1:size(data,2)
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
    for t=delay+1:size(data,2)
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
    
    
        
    FFs_H(:,repitition)=FF;
    FBs_H(:,repitition)=FB;
    FFs_iter_H(:,:,repitition)=FF_iter;
    FBs_iter_H(:,:,repitition)=FB_iter;
    
%     close all;
    ccc
    %% Tim's unexplained variance
    delay=5;
    
    for t=delay+1:size(data,2)
        
        % FF
        
        RDM_targ_past=nanmean(data(M+1:end,t-delay:t),2);
        RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
        
        RDM_targ_prsnt=data(M+1:end,t);
        RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
        
        RDM_srce_past=nanmean(data(1:M,t-delay:t),2);
        RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
        
        
        % Linear regression
        [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
        [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
        FF_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
        
        %     General Linear Model
        
        %     T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
        %     T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
        %     glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
        %     glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
        %
        %     FF_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
        
        % FB
        
        RDM_targ_past=nanmean(data(1:M,t-delay:t),2);
        RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
        
        RDM_targ_prsnt=data(1:M,t);
        RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
        
        RDM_srce_past=nanmean(data(M+1:end,t-delay:t),2);
        RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
        
        
        % Linear regression
        [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
        [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
        FB_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
        
        % General Linear Model
        
        %     T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
        %     T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
        %     glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
        %     glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
        %
        %     FB_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
    end
    
    
    for t=delay+1:size(data,2)
        for iteration=1:1000
            
            % FF
            
            RDM_targ_past=nanmean(data(M+1:end,t-delay:t),2);
            RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
            
            RDM_targ_prsnt=data(M+1:end,t);
            RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
            RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));
            
            
            RDM_srce_past=nanmean(data(1:M,t-delay:t),2);
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
            
            RDM_targ_past=nanmean(data(1:M,t-delay:t),2);
            RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
            
            RDM_targ_prsnt=data(1:M,t);
            RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
            RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));
            
            RDM_srce_past=nanmean(data(M+1:end,t-delay:t),2);
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
%         [t]
    end
    
    
%     figure;
%     subplot(3,1,1)
%     plot(FF_LR)
%     hold on;
    for t=delay+1:size(data,2)
        if FF_LR(t)>0
            FF_LR_sign(t)=1-sum(FF_LR(t)>FF_LR_iter(t,:))./size(FF_LR_iter,2);
        else
            FF_LR_sign(t)=1-sum(FF_LR(t)<FF_LR_iter(t,:))./size(FF_LR_iter,2);
        end
    end
    FF_LR_sign(1:delay)=1;
    try
        FF_LR_sign(delay+1:end)=mafdr(FF_LR_sign(delay+1:end));
    end
%     plot([1:size(FF_iter,1)],(FF_LR_sign<significance_threshold).*(sign(FF_LR)),'*')
%     
%     subplot(3,1,2)
%     plot(FB_LR)
%     hold on;
    for t=delay+1:size(data,2)
        if FB_LR(t)>0
            FB_LR_sign(t)=1-sum(FB_LR(t)>FB_LR_iter(t,:))./size(FB_LR_iter,2);
        else
            FB_LR_sign(t)=1-sum(FB_LR(t)<FB_LR_iter(t,:))./size(FB_LR_iter,2);
        end
    end
    FB_LR_sign(1:delay)=1;
    try
        FB_LR_sign(delay+1:end)=mafdr(FB_LR_sign(delay+1:end));
    end
%     plot([1:size(FB_LR_iter,1)],(FB_LR_sign<significance_threshold).*(sign(FB_LR)),'*')
    
    
    
%     subplot(3,1,3)
%     plot(FF_LR-FB_LR)
%     hold on;
    for t=delay+1:size(data,2)
        if (FF_LR(t)-FB_LR(t))>0
            FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))>(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
        else
            FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))<(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
        end
    end
    FF_FB_sign_LR(1:delay)=1;
    
    try
        FF_FB_sign_LR(delay+1:end)=mafdr(FF_FB_sign_LR(delay+1:end));
    end
    
%     plot([1:size(FF_iter,1)],(FF_FB_sign_LR<significance_threshold).*(sign(FF_LR-FB_LR)),'*')
    
    
    FFs_T(:,repitition)=FF_LR;
    FBs_T(:,repitition)=FB_LR;
    FFs_iter_T(:,:,repitition)=FF_LR_iter;
    FBs_iter_T(:,:,repitition)=FB_LR_iter;
    
    close all;
    [repitition]
%     save('Connectivity_repititions.mat','FFs','FBs','FFs_iter','FBs_iter','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H','FFs_T','FBs_T','FFs_iter_T','FBs_iter_T','MVGC')
end

%% Plotting the average
clear all;
% close all;
clc;
load('Connectivity_repititions.mat','FFs','FBs','FFs_iter','FBs_iter','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H','FFs_T','FBs_T','FFs_iter_T','FBs_iter_T','MVGC');
method=2; % 1= Erin; 2=Hamid, 3=TIm

significance_threshold=0.05;

repititions=size(FFs_T,2)-1;

if method==1
    FFs=FFs;
    FBs=FBs;
    FFs_iter=FFs_iter;
    FBs_iter=FBs_iter;
elseif method==2
    FFs=FFs_H;
    FBs=FBs_H;
    FFs_iter=FFs_iter_H;
    FBs_iter=FBs_iter_H;
elseif method==3
    FFs=FFs_T;
    FBs=FBs_T;
    FFs_iter=FFs_iter_T;
    FBs_iter=FBs_iter_T;
end

figure;
shadedErrorBar([1:100],nanmean(FFs'),nanstd(FFs')./sqrt(repititions),{'color',[0.1 0.1 0.8],'LineWidth',2},1)
hold on;
shadedErrorBar([1:100],nanmean(FBs'),nanstd(FBs')./sqrt(repititions),{'color',[0.8 0.1 0.1],'LineWidth',2},1)
shadedErrorBar([1:100],nanmean(FFs')-nanmean(FBs'),(nanmean(FFs')-nanmean(FBs'))./sqrt(repititions),{'color',[0.1 0.8 0.1],'LineWidth',2},1)

for time=1:100
    if nanmean(FFs(time,:))>0
        FFs_sign(time)=(1-sum(nanmean(FFs(time,:))>nanmean(FFs_iter(time,:,:),3))./1000);
    else
        FFs_sign(time)=(1-sum(nanmean(FFs(time,:))<nanmean(FFs_iter(time,:,:),3))./1000);
    end
    
    if nanmean(FBs(time,:))>0
        FBs_sign(time)=(1-sum(nanmean(FBs(time,:))>nanmean(FBs_iter(time,:,:),3))./1000);
    else
        FBs_sign(time)=(1-sum(nanmean(FBs(time,:))<nanmean(FBs_iter(time,:,:),3))./1000);
    end
    
    if nanmean(FFs(time,:)-FBs(time,:))>0
        FFs_minus_FBs_sign(time)=(1-sum(nanmean(FFs(time,:)-FBs(time,:))>nanmean(FFs_iter(time,:,:)-FBs_iter(time,:,:),3))./1000);
    else
        FFs_minus_FBs_sign(time)=(1-sum(nanmean(FFs(time,:)-FBs(time,:))<nanmean(FFs_iter(time,:,:)-FBs_iter(time,:,:),3))./1000);
    end    
end
plot([1:size(FFs_sign,2)],(FFs_sign<significance_threshold).*(sign(nanmean((FFs),2))')*0.3,'*b')
plot([1:size(FBs_sign,2)],(FBs_sign<significance_threshold).*(sign(nanmean((FBs),2))')*0.3,'*r')
plot([1:size(FFs_minus_FBs_sign,2)],(FFs_minus_FBs_sign<significance_threshold).*(sign(nanmean((FFs-FBs),2))')*0.5,'*g')

    
    