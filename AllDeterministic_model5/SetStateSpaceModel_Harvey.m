function [mZt, mGG, mT, mHH, va, mP,mH1] = SetStateSpaceModel_Harvey(vP,cn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign hyperparameters  %%vP=randn(2*k+2,1)
 
dsigma_eta = exp(2*vP(1 )).^0.5;
dphi =  vP(2)/sqrt(1 + vP(2)^2);

%-----------------------------------------------------------------------

mZ0=1;
mT  = dphi ;     mHH =  dsigma_eta^2 ; mH1= dsigma_eta;
va =0;  
mP = mHH;

  
mGG = 0;% dsigma_eta^2 ;%exp(2*vP(3+2*ck)).^0.5;
 
mZt = ones(cn,1)*mZ0 ;
%mZt = mZSeas; 
end
