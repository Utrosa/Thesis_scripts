function [sim_data,HAS_sim]=sim_combined_model1(intial_data,NK,xx,pRR,miu_map,sig_map,BIAS_COMBINED_REDUCE_FACTOR,NOISE_COMBINED_REDUCE_FACTOR,KernelWidth)

%fprintf('computing simulation data (combined model)...')
UN=size(intial_data,1);
%NK=%
NI=size(intial_data,3);
sim_data=nan(UN,NK,NI);
sim_data(:,1,:)=intial_data;


for K=2:NK
   % fprintf('.')
    for I=1:NI
        myx=sim_data(:,K-1,I);

       

        idx=find_indx_in_array(myx,xx);
        
        myy1=nan(size(myx));

        for ll=1:length(idx)
            ridx=roulette_distribution(pRR(:,idx(ll)),1);
            rval=xx(ridx);
            myy1(ll)=rval;
        end
        
        myx1=myy1;

        miu=interp1(xx,miu_map,myx1);
        sig=interp1(xx,sig_map,myx1);

        myy=myx1+BIAS_COMBINED_REDUCE_FACTOR*miu+randn(size(myx1)).*sig*NOISE_COMBINED_REDUCE_FACTOR;
        sim_data(:,K,I)=myy;

    end
end
%fprintf('\n');

HAS_sim=cell(NK,1);

for K=1:NK
 %   fprintf('.');
    HAS_sim{K}=compute_marginal(sim_data(:,K,:),xx,KernelWidth);
end
%fprintf('\n')







%%



