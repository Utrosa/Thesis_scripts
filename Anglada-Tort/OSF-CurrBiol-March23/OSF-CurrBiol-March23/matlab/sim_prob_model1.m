function [sim_data,HAS_sim]=sim_prob_model1(intial_data,NK,xx,pRR,KernelWidth)
%fprintf('computing simulation data (probabilitic model)...')
UN=size(intial_data,1);
NI=size(intial_data,3);
sim_data=nan(UN,NK,NI);
sim_data(:,1,:)=intial_data; % get the data from first iteartion as initial start

for K=2:NK
    %fprintf('.')
    for I=1:NI
        myx=sim_data(:,K-1,I); % get data from previous iteration
        idx=find_indx_in_array(myx,xx); % find the index of the point in xx
        
        myy=nan(size(myx));

        for ll=1:length(idx)
            ridx=roulette_distribution(pRR(:,idx(ll)),1); % randomized an index based on the conditional distribution
            rval=xx(ridx); % compute value from index
            myy(ll)=rval; % save
        end
        
        sim_data(:,K,I)=myy; % save final vector

    end
end
%fprintf('\n');

% computing the marginals of the simulations
%fprintf('computing smoothed kernels...')
HAS_sim=cell(NK,1);

for K=1:NK
    %fprintf('.');
    HAS_sim{K}=compute_marginal(sim_data(:,K,:),xx,KernelWidth);
end






%%



