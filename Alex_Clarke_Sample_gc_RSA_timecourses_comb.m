clear

%% Parameters
% Data parameters
cd('/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/MEG/RSA_ROI/RSA_timecourses_coord_RSA_sim_modelRDMs_DNN_summary_final_1_Mval4_TFbands_phz_60ms_circ');
mods = {'alexnet_convpca_all' 'mikenet_pca_allticks'};
ROIs = {'Occip' 'LpVTC' 'LATL' 'RpVTC' 'RATL'};
in = [2 4:5 7:15]; % subjects in
tw = [29:101]; % time-window to use

% GC parameters
ntrials   = 12;     % number of trials/subjects
nobs      = length(tw);   % number of observations per trial
regmode   = [];  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = [];  % information criteria regression mode ('OLS', 'LWR' or empty for default)
mordert    = '';%'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 8; %10;     % maximum model order for model order estimation
acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
fs        = 50;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)
seed      = 0;      % random seed (0 for unseeded)

morderI = [2];

% This part collects the data in X, then averages over freq band
for freq = 1:2  % freqband
    ii = 0;
    for mod = 1:length(mods)  % visual or semantic        
        %% Get data
        for roi = 1:length(ROIs)
            ii = ii+1;

            load(['TF_bands_phz_RS_timecourse_' ROIs{roi} 'Coord_source_80ms_sTW.mat']);
            a = rsa_out.(mods{mod})(in,:,:); a(isnan(a)) = 0;
            
            for sub = 1:size(a,1)
                X(ii,:,sub,freq) = permute(a(sub,freq,tw),[2 3 1]);  % ['region',time,subjects,freqband]
            end
        end
    end
end
    
X = mean(X,4); % average data across frequencies


%% GC part

    rng_seed(seed);
    nvars = size(X,1); % number of variables
    
    %% Model order estimation (<mvgc_schema.html#3 |A2|>)
    % Calculate information criteria up to specified maximum model order.
    ptic('\n*** tsdata_to_infocrit\n');
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
    ptoc('*** tsdata_to_infocrit took ');
    
    % Plot information criteria.
    figure; clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');
    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);
    
    if strcmp(mordert,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(mordert,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        morder = morderI
        fprintf('\nusing specified model order = %d\n',morderI);
    end
    
    %% VAR model estimation (<mvgc_schema.html#3 |A2|>)
    % Estimate VAR model of selected order from data.
    
    ptic('\n*** tsdata_to_var... ');
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    ptoc;
    
    % Check for failed regression
    assert(~isbad(A),'VAR estimation failed');
    
    %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)
    % The autocovariance sequence drives many Granger causality calculations (see
    % next section). Now we calculate the autocovariance sequence G according to the
    % VAR model, to as many lags as it takes to decay to below the numerical
    % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
    
    ptic('*** var_to_autocov... ');
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    ptoc;
    
    % The above routine does a LOT of error checking and issues useful diagnostics.
    % If there are problems with your data (e.g. non-stationarity, colinearity,
    % etc.) there's a good chance it'll show up at this point - and the diagnostics
    % may supply useful information as to what went wrong. It is thus essential to
    % report and check for errors here.
    
    var_info(info,true); % report results (and bail out on error)
    
    %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)
    
    ptic('*** autocov_to_pwcgc... ');
    F = autocov_to_pwcgc(G);
    ptoc;
    
    % Check for failed GC calculation
    assert(~isbad(F,false),'GC calculation failed');
    
    % Significance test using theoretical null distribution, adjusting for multiple
    % hypotheses.
    pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
    sig  = significance(pval,alpha,mhtc);            
    
    % Using permutations
    ptic('\n*** tsdata_to_mvgc_pwc_permtest\n');  
    FP = permtest_tsdata_to_pwcgc(X,morder,morder,5000);
    ptoc('*** tsdata_to_mvgc_pwc_permtest took ',[],1);              
    pval_p = empirical_pval(F,FP);
    sig_p  = significance(pval_p,alpha,mhtc);
      
    % Unbiased GC
    ubF = F - squeeze(mean(FP,1));
    pval_ub = mvgc_pval(ubF,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
    sig_ub  = significance(pval,alpha,mhtc);       
    
    % Plot time-domain causal graph, p-values and significance (theoretical)
    figure; clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])
    
    % Plot time-domain causal graph, p-values and significance (permutation)
    figure; clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval_p);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig_p);
    title(['Significant at p = ' num2str(alpha)])
    
    % Plot time-domain causal graph, p-values and significance (unbiased GC)
    figure; clf;
    subplot(1,3,1);
    plot_pw(ubF);
    title('Unbiased Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval_ub);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig_ub);
    title(['Significant at p = ' num2str(alpha)])
    
    % Collect results
    all_F = F;
    all_sig = sig;
    all_pval = pval;
    
    % For good measure we calculate Seth's causal density (cd) measure - the mean
    % pairwise-conditional causality. We don't have a theoretical sampling
    % distribution for this.
    cd = mean(F(~isnan(F)));
    
    fprintf('\ncausal density = %f\n',cd);
    
    %% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)
    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution - again, this only requires the autocovariance sequence.
    
    ptic('\n*** autocov_to_spwcgc... ');
    f = autocov_to_spwcgc(G,fres);
    ptoc;
    
    % Check for failed spectral GC calculation
    assert(~isbad(f,false),'spectral GC calculation failed');
    
    % Plot spectral causal graph.
    %figure; clf;
    %plot_spw(f,fs);
    
    %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)
    % Check that spectral causalities average (integrate) to time-domain
    % causalities, as they should according to theory.
    
    fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
    Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
    mad = maxabs(F-Fint);
    madthreshold = 1e-5;
    if mad < madthreshold
        fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
    else
        fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
    end
    