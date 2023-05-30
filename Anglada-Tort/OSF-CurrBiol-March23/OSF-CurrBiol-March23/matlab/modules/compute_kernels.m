function [FA,FAS]=compute_kernels(data,xgrid,ygrid,sigma_2d)

[X,Y]=meshgrid(xgrid,ygrid);


NU=size(data,1);
NK=size(data,2);
NI=size(data,3);

FA=zeros(size(X));
FAS=cell(NK,1);
for K=1:NK
    FAS{K}=zeros(size(X));
end
fprintf('Computing kernels...')
Sigma=eye(2)*(sigma_2d.^2);
for K=1:NK
    fprintf('.');
    for U=1:NU

        for I=1:(NI-1)
            mu=[data(U,K,I),data(U,K,I+1)];
            if isnan(sum(mu))
                continue
            end

            F = mvnpdf([X(:) Y(:)],mu,Sigma);

            F=F/sum(F);
            F = reshape(F,length(ygrid),length(xgrid));
            FA=FA+F;
            FAS{K}=FAS{K}+F;
        end

    end

end
fprintf('!\n');
%%
for K=1:NK
    FAS{K}=FAS{K}/sum(sum(FAS{K}) );
end
FA=FA/sum(FA(:));