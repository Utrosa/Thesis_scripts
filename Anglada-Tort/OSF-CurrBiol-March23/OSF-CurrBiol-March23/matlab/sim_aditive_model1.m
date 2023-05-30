function [sim_data,HAS_sim]=sim_aditive_model1(intial_data,NK,xx,miu_map,sig_map,KernelWidth)

%%% simulate the aditive model that depends only on interval size (miu_map)
%%% with noise depend on interval size (sig_map)

UN=size(intial_data,1);
NI=size(intial_data,3);
sim_data=nan(UN,NK,NI);
sim_data(:,1,:)=intial_data; % only first iteration

for K=2:NK
    for I=1:NI
        myx=sim_data(:,K-1,I); %use previous iteratio

        miu=interp1(xx,miu_map,myx); % convert to the current units
        sig=interp1(xx,sig_map,myx); % convert to the current units

        myy=miu+myx+randn(size(myx)).*sig; % run the model
        sim_data(:,K,I)=myy; % save data

    end
end

% compute marginals:
HAS_sim=cell(NK,1);
for K=1:NK
    HAS_sim{K}=compute_marginal(sim_data(:,K,:),xx,KernelWidth);
end
