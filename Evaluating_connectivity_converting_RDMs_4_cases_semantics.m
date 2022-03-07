clc;
clear all;
% close all;



samples_sizes=[100,200,300,200];
num_conditions=16;
Num_unq_RDM_cells=(num_conditions*num_conditions-num_conditions)./2;
Num_subjects=100;
iterations=100;

MVGC=nan*ones(2,2,Num_subjects,4);
FFs_E=nan*ones(760,Num_subjects);
FBs_E=nan*ones(760,Num_subjects);
FFs_H=nan*ones(760,Num_subjects,4);
FBs_H=nan*ones(760,Num_subjects,4);
FFs_iter_E=nan*ones(760,iterations,Num_subjects);
FBs_iter_E=nan*ones(760,iterations,Num_subjects);
FFs_iter_H=nan*ones(760,iterations,Num_subjects,4);
FBs_iter_H=nan*ones(760,iterations,Num_subjects,4);

for subject=1:Num_subjects
    RDM_source_t=nan*ones(Num_unq_RDM_cells,max(samples_sizes),2);
    for numb=1:3
        for time=1:max(samples_sizes)
            if numb==1
                
                X=nan*ones(num_conditions);
                for i=1:num_conditions
                    for j=i+1:num_conditions
                        X(i,j)=randn*2;
                    end
                end
                RDM_source_t(:,time,numb)=X(~isnan(X));
                
            elseif numb==2
                
                desired_semantics={'A','A','I','I','A','A','I','I','A','A','I','I','A','A','I','I'};
                X_semantics=nan*ones(num_conditions);
                
                for i=1:num_conditions
                    for j=i+1:num_conditions
                        if strcmp(desired_semantics{i},desired_semantics{j})==1
                            X_semantics(i,j)=1;
                        else
                            X_semantics(i,j)=0;
                        end
                    end
                end
                Model_RDM_semantics=X_semantics(~isnan(X_semantics));
                RDM_source_t(:,time,numb)=X_semantics(~isnan(X_semantics))+randn(Num_unq_RDM_cells,1);
                
            elseif numb==3
                
                desired_task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
                X_task=nan*ones(num_conditions);
                
                for i=1:num_conditions
                    for j=i+1:num_conditions
                        if strcmp(desired_task{i},desired_task{j})==1
                            X_task(i,j)=1;
                        else
                            X_task(i,j)=0;
                        end
                    end
                end
                
                Model_RDM_task=X_task(~isnan(X_task));
                RDM_source_t(:,time,numb)=X_task(~isnan(X_task))+randn(Num_unq_RDM_cells,1);
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
                        aloc(i,i+M)=abs(randn)*sigma;  % FF
                    end
                    
                end
                Arsig=[Arsig,aloc];
            end
            E=eye(M*2*P);
            AA=[Arsig;E(1:end-M*2,:)];
            lambda=eig(AA);
            lambdamax=max(abs(lambda));
            %                                 [c lambdamax]
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
    end
    data_t=data_t(:,2:end);
    RDM_R_source=data_t(1:M,:);
    
    
    data_t=nan*ones(Num_unq_RDM_cells*2,1);
    for current_phase=1:length(samples_sizes)
        N=samples_sizes(current_phase); % samples size
        subsample=randsample([1:max(samples_sizes)],N);
        RDM_source=RDM_source_t(:,subsample,:);
        %% Correlating RDMs
        P=10;  % model order
        M = Num_unq_RDM_cells; %number of cells in RDM;
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
                        aloc(i,i+M)=abs(randn)*sigma;  % FF
                    end
                    
                end
                Arsig=[Arsig,aloc];
            end
            E=eye(M*2*P);
            AA=[Arsig;E(1:end-M*2,:)];
            lambda=eig(AA);
            lambdamax=max(abs(lambda));
            %                                 [c lambdamax]
        end
        wind_sze_avod_spikes=5;
        if current_phase==1 %baseline
            x=[[squeeze(RDM_source(:,:,1)); squeeze(RDM_source(:,:,1))] zeros(M*2,N0)];
        elseif current_phase==2 % rise
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(RDM_source(:,:,3))+randn(size(RDM_source(:,:,3)))*5] zeros(M*2,N0)];
        elseif current_phase>2  && current_phase<length(samples_sizes) % fall
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(1-RDM_source(:,:,3))+randn(size(RDM_source(:,:,3)))*5] zeros(M*2,N0)];
        elseif current_phase==length(samples_sizes) % baseline
            %             x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-window_size:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(RDM_source(:,:,1))+randn(size(RDM_source(:,:,1)))*5] zeros(M*2,N0)];
            x=[[nanmean(data_t(1:Num_unq_RDM_cells,end-wind_sze_avod_spikes:end),2)+zeros(Num_unq_RDM_cells,N); squeeze(1-RDM_source(:,:,3))+randn(size(RDM_source(:,:,3)))*5] zeros(M*2,N0)];
        end
        y=x;
        for i=P+1:N+N0
            yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
            y(:,i)=Arsig*yloc+x(:,i);
        end
        data_t=horzcat(data_t,y(:,N0+1:end-N0));
    end
    data_t=data_t(:,2:end);
    RDM_R_destin=data_t(1:M,:);
    
    %% Generating a second RDM based on the current one out of noise
    for t=1:722
        Corr_to_Model(subject,t,1)=corr(Model_RDM_task,RDM_R_source(:,t));
        Corr_to_Model(subject,t,2)=corr(Model_RDM_semantics,RDM_R_source(:,t));
        Corr_to_Model(subject,t,3)=corr(Model_RDM_task,RDM_R_destin(:,t));
        Corr_to_Model(subject,t,4)=corr(Model_RDM_semantics,RDM_R_destin(:,t));
    end
    
    lambdamax=10;
    N=size(RDM_R_source,2);
    N0=ceil(N./20); %length of ignored start
    c=0;
    while (lambdamax>1 || lambdamax< 0.9)
        c=c+1;
        Arsig=[];
        for k=1:P
            aloc = zeros(M*2);
            for i=1:M*2
                aloc(i,i)=abs(randn)*sigma;  %everything related to its past
                if i<=M
                    aloc(i,i+M)=abs(randn)*sigma; % + FF
                end
            end
            Arsig=[Arsig,aloc];
        end
        E=eye(M*2*P);
        AA=[Arsig;E(1:end-M*2,:)];
        lambda=eig(AA);
        lambdamax=max(abs(lambda));
        %                 [case_number c lambdamax]
    end
    
    x=[[RDM_R_destin;RDM_R_source] zeros(M*2,N0)];
    y=x;
    for i=P+1:N+N0
        yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
        y(:,i)=Arsig*yloc+x(:,i);
    end
    data=y(:,N0+1:end-N0);
    
    RDM_R_source=squeeze(data(M+1:end,:));
    RDM_R_destin=squeeze(data(1:M,:));
    
    clearvars x data
    

    for t=1:722 % after transformation
        Corr_to_Model(subject,t,5)=corr(Model_RDM_task,RDM_R_source(:,t));
        Corr_to_Model(subject,t,6)=corr(Model_RDM_semantics,RDM_R_source(:,t));
        Corr_to_Model(subject,t,7)=corr(Model_RDM_task,RDM_R_destin(:,t));
        Corr_to_Model(subject,t,8)=corr(Model_RDM_semantics,RDM_R_destin(:,t));
    end
    
    %% Evaluating connectivity
    plotting=0;
    
    %% Erin's
    delay=P; % *5ms
    
    for t=delay+1:size(RDM_R_source,2)
        FF (t)= partialcorr(RDM_R_destin(:,t),nanmean(RDM_R_source(:,t-delay:t),2),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
        FB (t)= partialcorr(RDM_R_source(:,t),nanmean(RDM_R_destin(:,t-delay:t),2),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
    end
    
    
    % significance testing
    for t=delay+1:size(RDM_R_source,2)
        for iteration=1:iterations
            FF_iter (t,iteration)= partialcorr(randsample(RDM_R_destin(:,t),M),nanmean(RDM_R_source(:,t-delay:t),2),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
            FB_iter (t,iteration)= partialcorr(randsample(RDM_R_source(:,t),M),nanmean(RDM_R_destin(:,t-delay:t),2),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
            FF_minus_FB_iter (t,iteration)=FF_iter (t,iteration)-FB_iter (t,iteration);
        end
    end
    
    
    FFs_E(1:size(FF,2),subject)=FF;
    FBs_E(1:size(FB,2),subject)=FB;
    FFs_iter_E(1:size(FF_iter,1),:,subject)=FF_iter;
    FBs_iter_E(1:size(FF_iter,1),:,subject)=FB_iter;
    
    %% Tim's unexplained variance
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
    %         for iteration=1:iterations
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
    %
    %% Hamid's
    
    delay=P;
    for case_number=1:4
        if case_number==1
            
            for t=delay+1:size(RDM_R_source,2)
                
                Front(t)= corr(RDM_R_destin(:,t),Model_RDM_semantics,'row','complete');
                Front_minus_back(t)= partialcorr(RDM_R_destin(:,t),Model_RDM_semantics,nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                FF(t)=Front(t)-Front_minus_back(t);
                
                Back(t)= corr(RDM_R_source(:,t),Model_RDM_semantics,'row','complete');
                Back_minus_front(t)= partialcorr(RDM_R_source(:,t),Model_RDM_semantics,nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                FB(t)=Back(t)-Back_minus_front(t);
                
            end
            
            for t=delay+1:size(RDM_R_source,2)
                for iteration=1:iterations
                    
                    Front_iter= corr(RDM_R_destin(:,t),randsample(Model_RDM_semantics,M),'row','complete');
                    Front_minus_back_iter= partialcorr(RDM_R_destin(:,t),randsample(Model_RDM_semantics,M),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                    FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
                    
                    Back_iter= corr(RDM_R_source(:,t),randsample(Model_RDM_semantics,M),'row','complete');
                    Back_minus_front_iter= partialcorr(RDM_R_source(:,t),randsample(Model_RDM_semantics,M),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                    FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
                    
                    FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
                end
            end
            
        elseif case_number==2
            for t=delay+1:size(RDM_R_source,2)
                
                Front(t)= corr(RDM_R_destin(:,t),Model_RDM_task,'row','complete');
                Front_minus_back(t)= partialcorr(RDM_R_destin(:,t),Model_RDM_task,nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                FF(t)=Front(t)-Front_minus_back(t);
                
                Back(t)= corr(RDM_R_source(:,t),Model_RDM_task,'row','complete');
                Back_minus_front(t)= partialcorr(RDM_R_source(:,t),Model_RDM_task,nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                FB(t)=Back(t)-Back_minus_front(t);
                
            end
            
            for t=delay+1:size(RDM_R_source,2)
                for iteration=1:iterations
                    
                    Front_iter= corr(RDM_R_destin(:,t),randsample(Model_RDM_task,M),'row','complete');
                    Front_minus_back_iter= partialcorr(RDM_R_destin(:,t),randsample(Model_RDM_task,M),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                    FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
                    
                    Back_iter= corr(RDM_R_source(:,t),randsample(Model_RDM_task,M),'row','complete');
                    Back_minus_front_iter= partialcorr(RDM_R_source(:,t),randsample(Model_RDM_task,M),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                    FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
                    
                    FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
                end
            end
            
        elseif case_number==3
            
            for t=delay+1:size(RDM_R_source,2)
                
                Front(t)= corr(RDM_R_destin(:,t),Model_RDM_semantics,'row','complete');
                Front_minus_back(t)= partialcorr(RDM_R_destin(:,t),Model_RDM_semantics,nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                FF(t)=Front(t)-Front_minus_back(t);
                
                Back(t)= corr(RDM_R_source(:,t),Model_RDM_task,'row','complete');
                Back_minus_front(t)= partialcorr(RDM_R_source(:,t),Model_RDM_task,nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                FB(t)=Back(t)-Back_minus_front(t);
                
            end
            
            for t=delay+1:size(RDM_R_source,2)
                for iteration=1:iterations
                    
                    Front_iter= corr(RDM_R_destin(:,t),randsample(Model_RDM_semantics,M),'row','complete');
                    Front_minus_back_iter= partialcorr(RDM_R_destin(:,t),randsample(Model_RDM_semantics,M),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                    FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
                    
                    Back_iter= corr(RDM_R_source(:,t),randsample(Model_RDM_task,M),'row','complete');
                    Back_minus_front_iter= partialcorr(RDM_R_source(:,t),randsample(Model_RDM_task,M),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                    FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
                    
                    FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
                end
            end
            
         elseif case_number==4

             for t=delay+1:size(RDM_R_source,2)
                
                Front(t)= corr(RDM_R_destin(:,t),Model_RDM_task,'row','complete');
                Front_minus_back(t)= partialcorr(RDM_R_destin(:,t),Model_RDM_task,nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                FF(t)=Front(t)-Front_minus_back(t);
                
                Back(t)= corr(RDM_R_source(:,t),Model_RDM_semantics,'row','complete');
                Back_minus_front(t)= partialcorr(RDM_R_source(:,t),Model_RDM_semantics,nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                FB(t)=Back(t)-Back_minus_front(t);
                
            end
            
            for t=delay+1:size(RDM_R_source,2)
                for iteration=1:iterations
                    
                    Front_iter= corr(RDM_R_destin(:,t),randsample(Model_RDM_task,M),'row','complete');
                    Front_minus_back_iter= partialcorr(RDM_R_destin(:,t),randsample(Model_RDM_task,M),nanmean(RDM_R_source(:,t-delay:t),2),'row','complete');
                    FF_iter(t,iteration)=Front_iter-Front_minus_back_iter;
                    
                    Back_iter= corr(RDM_R_source(:,t),randsample(Model_RDM_semantics,M),'row','complete');
                    Back_minus_front_iter= partialcorr(RDM_R_source(:,t),randsample(Model_RDM_semantics,M),nanmean(RDM_R_destin(:,t-delay:t),2),'row','complete');
                    FB_iter(t,iteration)=Back_iter-Back_minus_front_iter;
                    
                    FF_minus_FB_iter(t,iteration)=FF_iter(t,iteration)-FB_iter(t,iteration);
                end
            end
            
            
        end
        
        FFs_H(1:size(FF,2),subject,case_number)=FF;
        FBs_H(1:size(FB,2),subject,case_number)=FB;
        FFs_iter_H(1:size(FF_iter,1),:,subject,case_number)=FF_iter;
        FBs_iter_H(1:size(FF_iter,1),:,subject,case_number)=FB_iter;
    end
    %% Alex Clark's
    for case_number=1:4
        for t=1:size(RDM_R_source,2)
            if case_number==1
                Front (t)= corr(RDM_R_destin(:,t),Model_RDM_semantics,'row','complete');
                Back (t)= corr(RDM_R_source(:,t),Model_RDM_semantics,'row','complete');
            elseif case_number==2
                Front (t)= corr(RDM_R_destin(:,t),Model_RDM_task,'row','complete');
                Back (t)= corr(RDM_R_source(:,t),Model_RDM_task,'row','complete');
            elseif case_number==3
                Front (t)= corr(RDM_R_destin(:,t),Model_RDM_semantics,'row','complete');
                Back (t)= corr(RDM_R_source(:,t),Model_RDM_task,'row','complete');
            elseif case_number==4
                Front (t)= corr(RDM_R_destin(:,t),Model_RDM_task,'row','complete');
                Back (t)= corr(RDM_R_source(:,t),Model_RDM_semantics,'row','complete');
            end
        end
        significance_level=0.05;
        plotting=0;
        MVGC(:,:,subject,case_number)=mvgc_Hamid([Front;Back],significance_level,plotting);
    end
    [subject case_number]
    save(['Evaluating_connectivity_converting_RDMs_4_cases_semantics.mat'],'MVGC','FFs_E','FBs_E','FFs_iter_E','FBs_iter_E','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H','Corr_to_Model')
end
%% Plotting the average
clc;
clear all;
% close all;

load(['Evaluating_connectivity_converting_RDMs_4_cases_semantics.mat'],'MVGC','FFs_E','FBs_E','FFs_iter_E','FBs_iter_E','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H','Corr_to_Model');


method=1; % 1= Erin; 2=Hamid; 3=Alex
case_number=1; %% 1= identity both regions, 2=task both regions, 3= frontal identity and occipital task, 4= opposite to 3


if method<3
    if method==1
        ylims(1,:)=[-0.1 0.3];
        ylims(2,:)=[-2 1];
        ylims(3,:)=[-1.5 2];
    elseif method==2
        ylims(1,:)=[-0.1 0.3];
        ylims(2,:)=[-1 1];
        ylims(3,:)=[-1.3 0.5];
    end
    Baseline(1,:)=[-1.25 0];
    Baseline(2,:)=[-0.35 0];
    
    subjects=size(FFs_H,2);
    
    if method==1
        FFs=FFs_E(:,:,case_number);
        FBs=FBs_E(:,:,case_number);
        FFs_iter=FFs_iter_E(:,:,:,case_number);
        FBs_iter=FBs_iter_E(:,:,:,case_number);
    elseif method==2
        FFs=FFs_H(:,:,case_number);
        FBs=FBs_H(:,:,case_number);
        FFs_iter=FFs_iter_H(:,:,:,case_number);
        FBs_iter=FBs_iter_H(:,:,:,case_number);
    end
    
    
    
    for time=1:size(FFs_H,1)
        Sigs(time,1)=bf.ttest(FFs(time,:),nanmean(FFs_iter(time,:,:),3));
        Sigs(time,2)=bf.ttest(FBs(time,:),nanmean(FBs_iter(time,:,:),3));
        Sigs(time,3)=bf.ttest(FFs(time,:)-FBs(time,:),nanmean(FFs_iter(time,:,:),3)-nanmean(FBs_iter(time,:,:),3));
    end
    for time=1:size(FFs_H,1)
        Effects=Sigs;
        for e=1:size(Effects,2)
            if Effects(time,e)>10
                Bayes(time,e)=2.5;
            elseif Effects(time,e)>3 && Effects(time,e)<=10
                Bayes(time,e)=1.5;
            elseif Effects(time,e)>1 && Effects(time,e)<=3
                Bayes(time,e)=0.5;
            elseif Effects(time,e)<1 && Effects(time,e)>=1/3
                Bayes(time,e)=-0.5;
            elseif Effects(time,e)<1/3 && Effects(time,e)>=1/10
                Bayes(time,e)=-1.5;
            elseif Effects(time,e)<1/10
                Bayes(time,e)=-2.5;
            end
        end
    end
    
    figure;
    models={'Identities','Tasks','Identities','Tasks'};
    titles={'Before','Before','After','After'};
    for model=1:4
        subplot(2,2,model)
        if model==1
            source_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,2))),squeeze(nanstd(Corr_to_Model(:,:,2)))./sqrt(subjects),{'color',[0.1 0.1 0.8],'LineWidth',2},1);
            hold on;
            destin_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,4))),squeeze(nanstd(Corr_to_Model(:,:,4)))./sqrt(subjects),{'color',[0.8 0.1 0.1],'LineWidth',2},1);
        elseif model==2
            source_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,1))),squeeze(nanstd(Corr_to_Model(:,:,1)))./sqrt(subjects),{'color',[0.1 0.1 0.8],'LineWidth',2},1);
            hold on;
            destin_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,3))),squeeze(nanstd(Corr_to_Model(:,:,3)))./sqrt(subjects),{'color',[0.8 0.1 0.1],'LineWidth',2},1);
        elseif model==3
            source_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,6))),squeeze(nanstd(Corr_to_Model(:,:,6)))./sqrt(subjects),{'color',[0.1 0.1 0.8],'LineWidth',2},1);
            hold on;
            destin_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,8))),squeeze(nanstd(Corr_to_Model(:,:,8)))./sqrt(subjects),{'color',[0.8 0.1 0.1],'LineWidth',2},1);
        elseif model==4
            source_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,5))),squeeze(nanstd(Corr_to_Model(:,:,5)))./sqrt(subjects),{'color',[0.1 0.1 0.8],'LineWidth',2},1);
            hold on;
            destin_cor=shadedErrorBar([1:size(Corr_to_Model,2)],squeeze(nanmean(Corr_to_Model(:,:,7))),squeeze(nanstd(Corr_to_Model(:,:,7)))./sqrt(subjects),{'color',[0.8 0.1 0.1],'LineWidth',2},1);            
        end
        line([1 700],[0 0],'linestyle','--','color','k')
        line([10 10],[ylims(1,:)],'linestyle','--','color','k')
        legend([source_cor.mainLine destin_cor.mainLine],{'Source','Destination'})
        ylabel(['Corr to ',models{model},' RDM'])
        title(titles{model})
        ylim(ylims(1,:))
        xlim([1 700])
        set(gca,'fontsize', 18);
    end
    
    figure;
    subplot(2,1,1)
    Feedf=shadedErrorBar([1:size(FFs_H,1)],nanmean(FFs'),nanstd(FFs')./sqrt(subjects),{'color',[0.1 0.8 0.8],'LineWidth',2},1);
    hold on;
    Feedb=shadedErrorBar([1:size(FFs_H,1)],nanmean(FBs'),nanstd(FBs')./sqrt(subjects),{'color',[0.8 0.8 0.1],'LineWidth',2},1);
    line([1 size(Corr_to_Model,2)],[0 0],'linestyle','--','color','k')
    line([10 10],ylims(2,:),'linestyle','--','color','k')
    ylabel('Information')
    hold on;
    
    steps=0.05;
    distans=1.5; % times step
    colors={[0.1 0.8 0.8],[0.8 0.8 0.1],[0.8 0.1 0.8]};
    for p=1:size(Bayes,2)-1
        for dots=1:size(Bayes,1)
            if Bayes(dots,p)==-0.5 || Bayes(dots,p)==0.5
                plots(p)=plot(dots,Bayes(dots,p).*steps+Baseline(method,1)-(p-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{p},'linewidth',2,'markersize',5);
            else
                plots(p)=plot(dots,Bayes(dots,p).*steps+Baseline(method,1)-(p-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{p},'Color',colors{p},'linewidth',2,'markersize',5);
            end
            hold on;
        end
        baseline_temp=Baseline(method,1)-(p-1)*(3*2+distans)*steps;
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp],'linestyle','-.','Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-2*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-3*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+2*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+3*steps,'Color',colors{p},'linewidth',2);
    end
    ylim(ylims(2,:))
    xlim([1 700])
    legend([Feedf.mainLine Feedb.mainLine],{'Feedforward','Feedback'})
    set(gca,'fontsize', 18);
    
    
    
    subplot(2,1,2)
    Diffr=shadedErrorBar([1:size(FFs_H,1)],nanmean(FFs')-nanmean(FBs'),(nanmean(FFs')-nanmean(FBs'))./sqrt(subjects),{'color',[0.8 0.1 0.8],'LineWidth',2},1);
    line([1 size(Corr_to_Model,2)],[0 0],'linestyle','--','color','k')
    line([10 10],[ylims(3,:)],'linestyle','--','color','k')
    ylabel('Information flow (+:FF -:FB)')
    xlabel('Time sample')
    hold on;
    
    steps=0.05;
    distans=1.5; % times step
    colors={[0.1 0.8 0.8],[0.8 0.8 0.1],[0.8 0.1 0.8]};
    for p=size(Bayes,2)
        for dots=1:size(Bayes,1)
            if Bayes(dots,p)==-0.5 || Bayes(dots,p)==0.5
                plots(p)=plot(dots,Bayes(dots,p).*steps+Baseline(method,2)-(p-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{p},'linewidth',2,'markersize',5);
            else
                plots(p)=plot(dots,Bayes(dots,p).*steps+Baseline(method,2)-(p-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{p},'Color',colors{p},'linewidth',2,'markersize',5);
            end
            hold on;
        end
        baseline_temp=Baseline(method,2)-(p-1)*(3*2+distans)*steps;
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp],'linestyle','-.','Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-2*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-3*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+2*steps,'Color',colors{p},'linewidth',2);
        line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+3*steps,'Color',colors{p},'linewidth',2);
    end
    ylim(ylims(3,:))
    xlim([1 700])
    legend([Diffr.mainLine],{'Difference (FF-FB)'})
    set(gca,'fontsize', 18);
    
else
    matrixx=squeeze(nanmean(MVGC(:,:,:,case_number),3));
    imagesc(matrixx,[0 1])
    xticks([1:4])
    xtickangle(45)
    xticklabels({'Front','Back'})
    xlabel('From')
    yticks([1:4])
    ytickangle(45)
    yticklabels({'Front','Back'})
    ylabel('To')
    CB=colorbar;
    for i=1:size(matrixx,1)
        for j=1:size(matrixx,2)
            text(j-0.1,i,num2str(matrixx(i,j)),'FontSize',18)
        end
    end
    title ('Information flow')
    set(gca,'fontsize', 18);
end

%%
clc;
clear all;
close all;
load(['Evaluating_connectivity_converting_RDMs_4_cases_semantics.mat'],'MVGC','FFs_E','FBs_E','FFs_iter_E','FBs_iter_E','FFs_H','FBs_H','FFs_iter_H','FBs_iter_H','Corr_to_Model');


method=2; % 1= Erin; 2=Hamid; 3=Alex
case_number=1; %% 1= identity both regions, 2=task both regions, 3= frontal identity and occipital task, 4= opposite to 3

colors={'k','r','g','b','m'};
ylims=[-0.5 0.25];
for case_number=1:4
    
    subjects=size(FFs_H,2);
    
    FFs=FFs_H(:,:,case_number);
    FBs=FBs_H(:,:,case_number);
    FFs_iter=FFs_iter_H(:,:,:,case_number);
    FBs_iter=FBs_iter_H(:,:,:,case_number);
    
    Diffr{case_number}=shadedErrorBar([1:size(FFs_H,1)],nanmean(FFs')-nanmean(FBs'),(nanmean(FFs')-nanmean(FBs'))./sqrt(subjects),{'color',colors{case_number},'LineWidth',2},1);
    hold on;
    
    
    for time=1:722
        Effects(time,case_number)=bf.ttest(FFs(time,:)-FBs(time,:),nanmean(FFs_iter(time,:,:),3)-nanmean(FBs_iter(time,:,:),3));
    end
end

for time=1:size(Effects,1)
    for e=1:size(Effects,2)
        if Effects(time,e)>10
            Bayes(time,e)=2.5;
        elseif Effects(time,e)>3 && Effects(time,e)<=10
            Bayes(time,e)=1.5;
        elseif Effects(time,e)>1 && Effects(time,e)<=3
            Bayes(time,e)=0.5;
        elseif Effects(time,e)<1 && Effects(time,e)>=1/3
            Bayes(time,e)=-0.5;
        elseif Effects(time,e)<1/3 && Effects(time,e)>=1/10
            Bayes(time,e)=-1.5;
        elseif Effects(time,e)<1/10
            Bayes(time,e)=-2.5;
        end
    end
end


Baseline=-0.3;
steps=0.005;
distans=1.5; % times step
for case_number=1:size(Bayes,2)
    for time=1:size(Bayes,1)
        if Bayes(time,case_number)==-0.5 || Bayes(time,case_number)==0.5
            plots(case_number)=plot(time,Bayes(time,case_number).*steps+Baseline-(case_number-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{case_number},'linewidth',2,'markersize',5);
        else
            plots(case_number)=plot(time,Bayes(time,case_number).*steps+Baseline-(case_number-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{case_number},'Color',colors{case_number},'linewidth',2,'markersize',5);
        end
        hold on;
    end
    baseline_temp=Baseline-(case_number-1)*(3*2+distans)*steps;
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp],'linestyle','-.','Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-steps,'Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-2*steps,'Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]-3*steps,'Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+steps,'Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+2*steps,'Color',colors{case_number},'linewidth',2);
    line([1 size(FFs_H,1)],[baseline_temp baseline_temp]+3*steps,'Color',colors{case_number},'linewidth',2);
end

ylim(ylims)
xlim([1 700])
line([1 size(Corr_to_Model,2)],[0 0],'linestyle','--','color','k')
line([10 10],[ylims],'linestyle','--','color','k')
ylabel('Information flow (+:FF -:FB)')
xlabel('Time sample')
% legend([Diffr{1}.mainLine Diffr{2}.mainLine Diffr{3}.mainLine Diffr{4}.mainLine Diffr{5}.mainLine],{'Both','Both Task out','Both Semantics out','Semantics Task out','Task Semantics out'})
legend([Diffr{1}.mainLine Diffr{2}.mainLine Diffr{3}.mainLine Diffr{4}.mainLine],{'Semantics across regions','Task across regions','Semantics FF and Task FB','Task FF and Semantics FB'})
set(gca,'fontsize', 18);




