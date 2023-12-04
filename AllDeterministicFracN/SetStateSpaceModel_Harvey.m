function [mZt, mGG, mT, mHH, va, mP,mH1] = SetStateSpaceModel_Harvey(vP,cn,cp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign hyperparameters  %%vP=randn(2*k+2,1)
 
dsigma_eta = exp(2*vP(1 )).^0.5;
dd =  exp(vP(2))./(1+exp(vP(2)));

%-----------------------------------------------------------------------

mH1= dsigma_eta;

opts = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter', 1000, 'MaxFunEvals', 1000,...
    'LargeScale', 'off' );

 
vP10 = zeros(2*cp,1);

f    = @(vP1)fARMA_Approx(dd, cn, cp, vP1); % function
[vP1, ~, ~, ~] = fminunc(f, vP10, opts );
vtAR    = fFisherInvTransform(vP1(1:cp));
vtMA    = fFisherInvTransform(vP1(cp+1:2*cp));
[mZ0, mT0, mHH0, va0, ~] = fSsfARMA(fReparAR(vtAR), fReparMA(vtMA), dsigma_eta^2);
mT = mT0; mHH = mHH0; mP = mHH0; va = va0;

 



mGG = 0;% exp(2*vP(3)).^0.5;
 
mZt = ones(cn,1)*mZ0 ;
%mZt = mZSeas; 
end
