clc;
clear all;
% close all;



samples_sizes=[100,200,300,200];
num_conditions=16;
Num_unq_RDM_cells=(num_conditions*num_conditions-num_conditions)./2;
Num_subjects=100;

for subject=1:Num_subjects
    
    RDM_source_t=nan*ones(Num_unq_RDM_cells,max(samples_sizes),2);
    for numb=1:2
        for time=1:max(samples_sizes)
            
            if numb==1 % random initial RDM
                X=nan*ones(num_conditions);
                for i=1:num_conditions
                    for j=i+1:num_conditions
                        X(i,j)=randn*2;
                    end
                end
                RDM_source_t(:,time,numb)=X(~isnan(X));
                
            elseif numb==2 % desired RDM
                
                desired={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
                X=nan*ones(num_conditions);
                
                for i=1:num_conditions
                    for j=i+1:num_conditions
                        if strcmp(desired{i},desired{j})==1
                            X(i,j)=1;
                        else
                            X(i,j)=0;
                        end
                    end
                end
                
                Model_RDM_desired=X(~isnan(X));
                RDM_source_t(:,time,numb)=X(~isnan(X))+randn(Num_unq_RDM_cells,1);
            end
        end
        
        %         subplot(1,size(RDM_source,3),numb)
        %         imagesc(X,[-1 1]);
    end
    
    data_t=nan*ones(Num_unq_RDM_cells*2,1);
    for current_phase=1:length(samples_sizes)
        N=samples_sizes(current_phase); % samples size
        subsample=randsample([1:max(samples_sizes)],N);
        RDM_source=RDM_source_t(:,subsample,:);
        
        %% Correlating RDMs
        
        P=10;  % model order
        M = Num_unq_RDM_cells; %number of cells in RDM;
        
        %         sigmas=[0.07 0.07 0.07 0.07];
        
        sigma = 0.07; %scale of random AR-parameters
        
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
                    
                    if i<=M
                        aloc(i,i+M)=abs(randn)*sigma;  % feedback
                    end
                    
                end
                Arsig=[Arsig,aloc];
            end
            E=eye(M*2*P);
            AA=[Arsig;E(1:end-M*2,:)];
            lambda=eig(AA);
            lambdamax=max(abs(lambda));
            %                     [c lambdamax]
        end
        
        wind_sze_avod_spikes=5;
        if current_phase==1 %baseline
            x=[[squeeze(RDM_source(:,:,1)); squeeze(RDM_source(:,:,1))] zeros(M*2,N0)];
        elseif current_phase==2 % rise
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(RDM_source(:,:,2))+randn(size(RDM_source(:,:,2)))*5] zeros(M*2,N0)];
        elseif current_phase>2  && current_phase<length(samples_sizes) % fall
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(1-RDM_source(:,:,2))+randn(size(RDM_source(:,:,2)))*5] zeros(M*2,N0)];
        elseif current_phase==length(samples_sizes) % baseline
            %             x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-window_size:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(RDM_source(:,:,1))+randn(size(RDM_source(:,:,1)))*5] zeros(M*2,N0)];
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(1-RDM_source(:,:,2))+randn(size(RDM_source(:,:,2)))*5] zeros(M*2,N0)];
        end
        
        y=x;
        
        for i=P+1:N+N0
            yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
            y(:,i)=Arsig*yloc+x(:,i);
        end
        data_t=horzcat(data_t,y(:,N0+1:end-N0));
        %
        %         plot(y(5,N0+1:end-N0))
        %         hold on;
        %
    end
    data_t=data_t(:,2:end);
    %% generating original and delayed RDMs
    RDM_R1=data_t(1:M,:);
    delay=P;
    RDM_R2=[RDM_R1(:,1:delay) RDM_R1(:,1:end-delay)];
%     [subject]
    
    %     save('Corr_RDM_with_models.mat','Model_RDM_desired','RDM_R1','RDM_R2')
    
    %% Plotting RDM-to-model correlations for verification
    % clear all;
    % close all;
    % clc;
    % load('Corr_RDM_with_models.mat')
    % for subject=1:size(RDM_R1,3)
    %     for time=1:size(RDM_R1,2)
    %         cors_to_desired(time,subject,1)=corr(squeeze(RDM_R1(:,time,subject)),Model_RDM_desired);
    %         cors_to_desired(time,subject,2)=corr(squeeze(RDM_R2(:,time,subject)),Model_RDM_desired);
    %     end
    % end
    % figure;
    % plt=shadedErrorBar([1:size(cors_to_desired,1)],nanmean(cors_to_desired(:,:,1)'),nanstd(cors_to_desired(:,:,1)')./sqrt(size(cors_to_desired(:,:,1),2)),{'color',[0.6 0.1 0.1],'LineWidth',2},1);
    % hold on;
    % plt2=shadedErrorBar([1:size(cors_to_desired,1)],nanmean(cors_to_desired(:,:,2)'),nanstd(cors_to_desired(:,:,2)')./sqrt(size(cors_to_desired(:,:,2),2)),{'color',[0.1 0.1 0.6],'LineWidth',2},1);
    % legend([plt.mainLine plt2.mainLine],'Region-1','Region-2')
    %
    % % finding relevant matrix indices
    % RDM_t=nan*ones(num_conditions);
    % for i=1:num_conditions
    %     for j=i+1:num_conditions
    %         RDM_t(i,j)=1;
    %     end
    % end
    % indices_for_RDMs=find(RDM_t);
    
    % % collapsing across subjects
    % RDM_R1=squeeze(nanmean(RDM_R1,3));
    % RDM_R2=squeeze(nanmean(RDM_R2,3));
    %% Evaluating connectivity
    plotting=0;
    
    %% Erin's
    delay=P; % *5ms
    
    for t=delay+1:size(RDM_R1,2)
        FF (t)= partialcorr(RDM_R2(:,t),nanmean(RDM_R1(:,t-delay:t),2),nanmean(RDM_R2(:,t-delay:t),2),'row','complete');
        FB (t)= partialcorr(RDM_R1(:,t),nanmean(RDM_R2(:,t-delay:t),2),nanmean(RDM_R1(:,t-delay:t),2),'row','complete');
    end
    
    
    % significance testing
    for t=delay+1:size(RDM_R1,2)
        for iteration=1:100
            FF_iter (t,iteration)= partialcorr(randsample(RDM_R2(:,t),M),nanmean(RDM_R1(:,t-delay:t),2),nanmean(RDM_R2(:,t-delay:t),2),'row','complete');
            FB_iter (t,iteration)= partialcorr(randsample(RDM_R1(:,t),M),nanmean(RDM_R2(:,t-delay:t),2),nanmean(RDM_R1(:,t-delay:t),2),'row','complete');
            FF_minus_FB_iter (t,iteration)=FF_iter (t,iteration)-FB_iter (t,iteration);
        end
        %         [t]
    end
    
    
    
    if plotting
        significance_threshold=0.05;
        figure;
        subplot(3,1,1)
        plot(FF)
        hold on;
        for t=delay+1:size(RDM_R1,2)
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
        for t=delay+1:size(RDM_R1,2)
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
        for t=delay+1:size(RDM_R1,2)
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
        
    end
    FFs_E(:,subject)=FF;
    FBs_E(:,subject)=FB;
    FFs_iter_E(:,:,subject)=FF_iter;
    FBs_iter_E(:,:,subject)=FB_iter;
    
%     %% Tim's unexplained variance
%     delay=P;
%     
%     for t=delay+1:size(RDM_R1,2)
%         
%         % FF
%         
%         RDM_targ_past=nanmean(RDM_R2(:,t-delay:t),2);
%         RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
%         
%         RDM_targ_prsnt=RDM_R2(:,t);
%         RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
%         
%         RDM_srce_past=nanmean(RDM_R1(:,t-delay:t),2);
%         RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
%         
%         
%         % Linear regression
%         [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
%         [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
%         FF_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
%         
%         %     General Linear Model
%         
%         %     T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%         %     T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%         %     glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         %     glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         %
%         %     FF_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
%         
%         % FB
%         
%         RDM_targ_past=nanmean(RDM_R1(:,t-delay:t),2);
%         RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
%         
%         RDM_targ_prsnt=RDM_R1(:,t);
%         RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
%         
%         RDM_srce_past=nanmean(RDM_R2(:,t-delay:t),2);
%         RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
%         
%         
%         % Linear regression
%         [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
%         [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
%         FB_LR (t)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
%         
%         % General Linear Model
%         
%         %     T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%         %     T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%         %     glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         %     glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%         %
%         %     FB_GLM (t)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
%     end
%     
%     
%     for t=delay+1:size(RDM_R1,2)
%         for iteration=1:100
%             
%             % FF
%             
%             RDM_targ_past=nanmean(RDM_R2(:,t-delay:t),2);
%             RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
%             
%             RDM_targ_prsnt=RDM_R2(:,t);
%             RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
%             RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));
%             
%             
%             RDM_srce_past=nanmean(RDM_R1(:,t-delay:t),2);
%             RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
%             
%             
%             % Linear regression
%             [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
%             [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
%             FF_LR_iter (t,iteration)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
%             
%             %         % General Linear Model
%             %
%             %         T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%             %         T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%             %         glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%             %         glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%             %
%             %         FF_GLM_iter (t,iteration)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
%             
%             % FB
%             
%             RDM_targ_past=nanmean(RDM_R1(:,t-delay:t),2);
%             RDM_targ_past=RDM_targ_past(~isnan(RDM_targ_past));
%             
%             RDM_targ_prsnt=RDM_R1(:,t);
%             RDM_targ_prsnt=RDM_targ_prsnt(~isnan(RDM_targ_prsnt));
%             RDM_targ_prsnt=randsample(RDM_targ_prsnt,length(RDM_targ_prsnt));
%             
%             RDM_srce_past=nanmean(RDM_R2(:,t-delay:t),2);
%             RDM_srce_past=RDM_srce_past(~isnan(RDM_srce_past));
%             
%             
%             % Linear regression
%             [~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
%             [~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
%             FB_LR_iter (t,iteration)=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare));
%             
%             %         % General Linear Model
%             
%             %         T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
%             %         T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
%             %         glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%             %         glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects');
%             %
%             %         FB_GLM_iter (t,iteration)=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted));
%         end
%         %         [t]
%     end
%     
%     if plotting
%         figure;
%         subplot(3,1,1)
%         plot(FF_LR)
%         hold on;
%         for t=delay+1:size(RDM_R1,2)
%             if FF_LR(t)>0
%                 FF_LR_sign(t)=1-sum(FF_LR(t)>FF_LR_iter(t,:))./size(FF_LR_iter,2);
%             else
%                 FF_LR_sign(t)=1-sum(FF_LR(t)<FF_LR_iter(t,:))./size(FF_LR_iter,2);
%             end
%         end
%         FF_LR_sign(1:delay)=1;
%         try
%             FF_LR_sign(delay+1:end)=mafdr(FF_LR_sign(delay+1:end));
%         end
%         plot([1:size(FF_iter,1)],(FF_LR_sign<significance_threshold).*(sign(FF_LR)),'*')
%         
%         subplot(3,1,2)
%         plot(FB_LR)
%         hold on;
%         for t=delay+1:size(RDM_R1,2)
%             if FB_LR(t)>0
%                 FB_LR_sign(t)=1-sum(FB_LR(t)>FB_LR_iter(t,:))./size(FB_LR_iter,2);
%             else
%                 FB_LR_sign(t)=1-sum(FB_LR(t)<FB_LR_iter(t,:))./size(FB_LR_iter,2);
%             end
%         end
%         FB_LR_sign(1:delay)=1;
%         try
%             FB_LR_sign(delay+1:end)=mafdr(FB_LR_sign(delay+1:end));
%         end
%         plot([1:size(FB_LR_iter,1)],(FB_LR_sign<significance_threshold).*(sign(FB_LR)),'*')
%         
%         
%         
%         subplot(3,1,3)
%         plot(FF_LR-FB_LR)
%         hold on;
%         for t=delay+1:size(RDM_R1,2)
%             if (FF_LR(t)-FB_LR(t))>0
%                 FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))>(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
%             else
%                 FF_FB_sign_LR(t)=1-sum((FF_LR(t)-FB_LR(t))<(FF_LR_iter(t,:)-FB_LR_iter(t,:)))./size(FF_LR_iter,2);
%             end
%         end
%         FF_FB_sign_LR(1:delay)=1;
%         
%         try
%             FF_FB_sign_LR(delay+1:end)=mafdr(FF_FB_sign_LR(delay+1:end));
%         end
%         
%         plot([1:size(FF_iter,1)],(FF_FB_sign_LR<significance_threshold).*(sign(FF_LR-FB_LR)),'*')
%     end
%     
%     FFs_T(:,subject)=FF_LR;
%     FBs_T(:,subject)=FB_LR;
%     FFs_iter_T(:,:,subject)=FF_LR_iter;
%     FBs_iter_T(:,:,subject)=FB_LR_iter;
    
    %% Hamid's
    
    delay=P;
    
    
    for t=delay+1:size(RDM_R1,2)
        
        Front(t)= corr(RDM_R2(:,t),Model_RDM_desired,'row','complete');
        Front_minus_back(t)= partialcorr(RDM_R2(:,t),Model_RDM_desired,nanmean(RDM_R1(:,t-delay:t),2),'row','complete');
        FF(t)=Front(t)-Front_minus_back(t);
        
        Back(t)= corr(RDM_R1(:,t),Model_RDM_desired,'row','complete');
        Back_minus_front(t)= partialcorr(RDM_R1(:,t),Model_RDM_desired,nanmean(RDM_R2(:,t-delay:t),2),'row','complete');
        FB(t)=Back(t)-Back_minus_front(t);
        
    end
    
    for t=delay+1:size(RDM_R1,2)
        for iteration=1:100
            
            Front_iter= corr(RDM_R2(:,t),randsample(Model_RDM_desired,M),'row','complete');
            Front_minus_back_iter= partialcorr(RDM_R2(:,t),randsample(Model_RDM_desired,M),nanmean(RDM_R1(:,t-delay:t),2),'row','complete');
            FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
            
            Back_iter= corr(RDM_R1(:,t),randsample(Model_RDM_desired,M),'row','complete');
            Back_minus_front_iter= partialcorr(RDM_R1(:,t),randsample(Model_RDM_desired,M),nanmean(RDM_R2(:,t-delay:t),2),'row','complete');
            FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
            
            FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
        end
    end
    
    if plotting==1
        significance_threshold=0.05;
        
        figure;
        subplot(3,1,1)
        plot(FF)
        hold on;
        for t=delay+1:size(RDM_R1,2)
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
        for t=delay+1:size(RDM_R1,2)
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
        for t=delay+1:size(RDM_R1,2)
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
        
    end
    
    FFs_H(:,subject)=FF;
    FBs_H(:,subject)=FB;
    FFs_iter_H(:,:,subject)=FF_iter;
    FBs_iter_H(:,:,subject)=FB_iter;
    
    %% Alex Clark's
    
    for t=1:size(RDM_R1,2)
        Front_task_semantics (t)= corr(RDM_R2(:,t),Model_RDM_desired,'row','complete');
        Back_task_semantics (t)= corr(RDM_R1(:,t),Model_RDM_desired,'row','complete');
    end
    % figure;
    % plot([Front_task;Front_conditions;Back_task;Back_conditions]')
    significance_level=0.05;
    plotting=0;
    MVGC(:,:,subject)=mvgc_Hamid([Front_task_semantics;Back_task_semantics],significance_level,plotting);
    
    %     close all;
    [subject]
%     save('Evaluating_connectivity_shifted_RDMs.mat','MVGC','FFs_E','FBs_E','FFs_iter_E','FBs_iter_E','FFs_T','FBs_T','FFs_iter_T','FBs_iter_T','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H')
    save('Evaluating_connectivity_shifted_RDMs.mat','MVGC','FFs_E','FBs_E','FFs_iter_E','FBs_iter_E','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H')
end
%% Plotting the average
% if plotting==1
%     clear all;
%     % close all;
%     clc;
%     method=2; % 1= Erin; 2=Hamid, 3=TIm
%     load('Connectivity_subjects_corrected2_FB_only_Hamid.mat')
%
%     FFs2_H=FFs_H(:,1:50);
%     FBs2_H=FBs_H(:,1:50);
%     FFs2_iter_H=FFs_iter_H(:,:,1:50);
%     FBs2_iter_H=FBs_iter_H(:,:,1:50);
%
%
%     load('Connectivity_subjects_corrected_FB_only_Hamid.mat')
%
%
%     FFs_H(:,1:50)=FFs2_H;
%     FBs_H(:,1:50)=FBs2_H;
%     FFs_iter_H(:,:,1:50)=FFs2_iter_H;
%     FBs_iter_H(:,:,1:50)=FBs2_iter_H;
%
%     clearvars 'FFs2' 'FBs2' 'FFs2_iter' 'FBs2_iter' 'FFs2_H' 'FBs2_H' 'FFs2_iter_H' 'FBs2_iter_H' 'FFs2_T' 'FBs2_T' 'FFs2_iter_T' 'FBs2_iter_T' 'MVGC2'
%
%     significance_threshold=0.05;
%
%     subjects=size(FFs_H,2)-1;
%
%     if method==1
%         FFs=FFs;
%         FBs=FBs;
%         FFs_iter=FFs_iter;
%         FBs_iter=FBs_iter;
%     elseif method==2
%         FFs=FFs_H;
%         FBs=FBs_H;
%         FFs_iter=FFs_iter_H;
%         FBs_iter=FBs_iter_H;
%     elseif method==3
%         FFs=FFs_T;
%         FBs=FBs_T;
%         FFs_iter=FFs_iter_T;
%         FBs_iter=FBs_iter_T;
%     end
%
%     figure;
%     shadedErrorBar([1:size(FFs_H,1)],nanmean(FFs'),nanstd(FFs')./sqrt(subjects),{'color',[0.1 0.1 0.8],'LineWidth',2},1)
%     hold on;
%     shadedErrorBar([1:size(FFs_H,1)],nanmean(FBs'),nanstd(FBs')./sqrt(subjects),{'color',[0.8 0.1 0.1],'LineWidth',2},1)
%     shadedErrorBar([1:size(FFs_H,1)],nanmean(FFs')-nanmean(FBs'),(nanmean(FFs')-nanmean(FBs'))./sqrt(subjects),{'color',[0.1 0.8 0.1],'LineWidth',2},1)
%
%     for time=1:size(FFs_H,1)
%         if nanmean(FFs(time,:))>0
%             FFs_sign(time)=(1-sum(nanmean(FFs(time,:))>nanmean(FFs_iter(time,:,:),3))./size(FFs_iter,2));
%         else
%             FFs_sign(time)=(1-sum(nanmean(FFs(time,:))<nanmean(FFs_iter(time,:,:),3))./size(FFs_iter,2));
%         end
%
%         if nanmean(FBs(time,:))>0
%             FBs_sign(time)=(1-sum(nanmean(FBs(time,:))>nanmean(FBs_iter(time,:,:),3))./size(FBs_iter,2));
%         else
%             FBs_sign(time)=(1-sum(nanmean(FBs(time,:))<nanmean(FBs_iter(time,:,:),3))./size(FBs_iter,2));
%         end
%
%         if nanmean(FFs(time,:)-FBs(time,:))>0
%             FFs_minus_FBs_sign(time)=(1-sum(nanmean(FFs(time,:)-FBs(time,:))>nanmean(FFs_iter(time,:,:)-FBs_iter(time,:,:),3))./size(FFs_iter,2));
%         else
%             FFs_minus_FBs_sign(time)=(1-sum(nanmean(FFs(time,:)-FBs(time,:))<nanmean(FFs_iter(time,:,:)-FBs_iter(time,:,:),3))./size(FFs_iter,2));
%         end
%     end
%     plot([1:size(FFs_sign,2)],(FFs_sign<significance_threshold).*(sign(nanmean((FFs),2))')*0.3,'*b')
%     plot([1:size(FBs_sign,2)],(FBs_sign<significance_threshold).*(sign(nanmean((FBs),2))')*0.25,'*r')
%     plot([1:size(FFs_minus_FBs_sign,2)],(FFs_minus_FBs_sign<significance_threshold).*(sign(nanmean((FFs-FBs),2))')*0.5,'*g')
%
% end
%
