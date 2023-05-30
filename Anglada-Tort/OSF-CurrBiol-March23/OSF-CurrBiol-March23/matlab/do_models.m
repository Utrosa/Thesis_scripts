function [sim_data,HAS_sim]=do_models(todo_letter,p1,p2,data1,NK1,xxx,miu_mapi,sig_mapi,pyxi,p0,KernelWidth)

switch todo_letter
    case 'M' % interval size model
        NOISE_COMBINED_REDUCE_FACTOR=p1; %param
        BIAS_COMBINED_REDUCE_FACTOR=p2; % param
        [sim_data,HAS_sim]=sim_aditive_model1(data1(:,1,:),NK1,xxx,miu_mapi*BIAS_COMBINED_REDUCE_FACTOR,sig_mapi*NOISE_COMBINED_REDUCE_FACTOR,KernelWidth); % simualte motor (adding the possiblity to sloghtly deviate from the original values)
    case 'P'% pref model
        gamma=p1; % params
        sigma_sen=p2; % params
        pUtil=exp(gamma*pyxi); % utility from Harrison et al. 2020 -> exponent of the rating.
        pUtil=pUtil/sum(pUtil); % normalized
        pR=pUtil';  % prior = exp (gamma*rating)
        tic()
        [pRR,~]=compute_bayes_prior_model(pR,sigma_sen,p0,xxx,NK1); % analytical compuation of the transition probabilties of the model

        [sim_data,HAS_sim]=sim_prob_model1(data1(:,1,:),NK1,xxx,pRR,KernelWidth); % actual prference model
    otherwise
        fprintf('Wrong todo value\n');
        assert(1==0);
end