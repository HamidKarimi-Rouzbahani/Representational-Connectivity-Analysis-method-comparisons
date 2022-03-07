clc;
clear all;
close all;

direction=1; % 1= (1=>2) FF; 2= (2=>1) FB; 3= (1<=>2) Bidirectional


repititions=3;
for repitition=1:repititions
    
    A=randn(10,16,200); % channel_in_region * condition_in_region * time sample
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
                        X(i,j)=0.9+randn*0.1;
                    else
                        X(i,j)=randn*0.1;
                    end
                end
            end
            
        else % frontal
            
            %         Identity={'A','B','C','D','E','F','G','H','A','B','C','D','E','F','G','H'};
            %             Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
            Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
            
            for i=1:size(A,2)
                for j=1:size(A,2)
                    if strcmp(Task{i},Task{j})==1
                        X(i,j)=0.1+randn*0.1;
                    else
                        X(i,j)=randn*0.1;
                    end
                end
            end
%                         X=randn(size(A,2))*0.1;

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
            imagesc(cors_source,[-1 1]);
        % figure
        %         imagesc(squeeze(nanmean(cors_save(:,:,:,numb),3)))
    end
    %% Correlating RDMs in time
    clc;
    
    
    N=size(A,3); % samples size
    P=10;  % model order
    M = size(RDM_source,1); %number of cells in RDM;
    
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
    ccc
    x=[[squeeze(RDM_source(:,:,1));squeeze(RDM_source(:,:,2))] zeros(M*2,N0)];
    y=x;
    for i=P+1:N+N0
        yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
        y(:,i)=Arsig*yloc+x(:,i);
        sums(i)=sum(y(1:M,i)-y(M+1:end,i));
    end
    data=y(:,N0+1:end);
    
    %% Evaluating if the correlated RDMs are correlated using simple time-shifted correlation
    
    %     delays=P;
    %
    %     cors_FF=nan*ones(size(data,2),length(delays));
    %     cors_FB=nan*ones(size(data,2),length(delays));
    %     cors3=nan*ones(size(data,2),length(delays));
    %     cors4=nan*ones(size(data,2),length(delays));
    %
    %     d=0;
    %     for delay=delays
    %         d=d+1;
    %         for time=1+delay:size(data,2)-delay
    %             %         cors_FF(time,d)=corr(data(M+1:end,time),nanmean(data(1:M,time-delay:time),2),'Type','Spearman'); % occipital-frontal : feedforward
    %             %         cors_FB(time,d)=corr(data(1:M,time),nanmean(data(M+1:end,time-delay:time),2),'Type','Spearman'); % frontal-occipital : feedback
    %             cors_FF(time,d)=corr(data(M+1:end,time),nanmean(data(1:M,time-delay:time),2)); % occipital-frontal : feedforward
    %             cors_FB(time,d)=corr(data(1:M,time),nanmean(data(M+1:end,time-delay:time),2)); % frontal-occipital : feedback
    %         end
    %     end
    %     figure;
    %     plot(cors_FF)
    %     hold on;
    %     plot(cors_FB)
    
    
    %% Generating model RDMs
    
    Model_RDM_identities=nan*ones(size(A,2));
    Model_RDM_task=nan*ones(size(A,2));
    
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
    %     Semantic={'M','M','N','N','M','M','N','N','M','M','N','N','M','M','N','N'};
    Task={'T','T','T','T','T','T','T','T','D','D','D','D','D','D','D','D'};
    
    for i=1:size(A,2)
        for j=i+1:size(A,2)
            if strcmp(Task{i},Task{j})==1
                Model_RDM_task(i,j)=1;
            else
                Model_RDM_task(i,j)=0;
            end
        end
    end
    %     figure;
    %     imagesc(Model_RDM_task_semantics)
    Model_RDM_task=Model_RDM_task(~isnan(Model_RDM_task));
    
    %%
    for time=1:size(data,2)
        cors_Front_identities(time,repitition)=corr(data(M+1:end,time),Model_RDM_identities); % occipital-frontal : feedforward
        cors_Front_task(time,repitition)=corr(data(M+1:end,time),Model_RDM_task); % occipital-frontal : feedforward
        cors_Back_identities(time,repitition)=corr(data(1:M,time),Model_RDM_identities); % frontal-occipital : feedback
        cors_Back_task(time,repitition)=corr(data(1:M,time),Model_RDM_task); % frontal-occipital : feedback
    end
     
    [repitition]
    save('Corr_RDM_with_models.mat','cors_Front_identities','cors_Front_task','cors_Back_identities','cors_Back_task')   
end

%% Plotting the average
% clear all;
close all;
clc;
load('Corr_RDM_with_models.mat')
A=shadedErrorBar([1:size(cors_Front_identities,1)],nanmean(cors_Front_identities'),nanstd(cors_Front_identities')./sqrt(size(cors_Front_identities,2)),{'color',[0.1 0.1 0.8],'LineWidth',2},1);
hold on;
B=shadedErrorBar([1:size(cors_Front_task,1)],nanmean(cors_Front_task'),nanstd(cors_Front_task')./sqrt(size(cors_Front_task,2)),{'color',[0.1 0.1 0.6],'LineWidth',2},1);
C=shadedErrorBar([1:size(cors_Back_identities,1)],nanmean(cors_Back_identities'),nanstd(cors_Back_identities')./sqrt(size(cors_Back_identities,2)),{'color',[0.8 0.1 0.1],'LineWidth',2},1);
D=shadedErrorBar([1:size(cors_Back_task,1)],nanmean(cors_Back_task'),nanstd(cors_Back_task')./sqrt(size(cors_Back_task,2)),{'color',[0.6 0.1 0.1],'LineWidth',2},1);
legend([A.mainLine,B.mainLine,C.mainLine,D.mainLine],'Front-ident','Front-task','Back-ident','Back-task')

