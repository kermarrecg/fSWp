function [v1s,v2s]= f_fSW_BoxCox(vr, dLambda, cp, mXreg,vv2)

vv1= [10^(-10)  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ];

%vv1= [10^(-10)  0.25   0.5  0.75 ];

cnv1=length(vv1);
cnv2=length(vv2);
 
mloglik =NaN(cnv2,cnv1);
opts = optimset('Display','Iter','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter',400, 'MaxFunEvals', 400,...
    'LargeScale', 'off' );
ck=length(dLambda);

 
for ii=1:cnv2
    v2 = vv2(ii)
    for jj=1:cnv1  
        
        v1 = vv1(jj)
        [vP0,~]=fgrid_InitCondsBNJACOB(vr, dLambda , mXreg,cp,v1,v2  );
        fU = @(vP) fSsfLogLik_HarveyJACOB(vP, vr, dLambda, cp, mXreg,v1,v2);
        [~, dLogLik, ~, ~] = fminsearch(fU, vP0, opts);
        mloglik(ii,jj)=dLogLik
 
    end
end

[ind1,ind2]=find(mloglik==min(min(mloglik)));

v2s=vv2(ind1);
v1s=vv1(ind2);

end