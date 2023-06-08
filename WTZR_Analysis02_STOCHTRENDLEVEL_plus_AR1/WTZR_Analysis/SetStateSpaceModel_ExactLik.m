function [mZt, mGG, mT, mHH, va, mP ] = SetStateSpaceModel_ExactLik(vP,iTrend, mZSeas,mXreg,vBeta0,cp)


[mZt_2, mGG_2, mT_2, mHH_2, mW,va_2, mP_2 , mW0] = fSetSSModel(vP, iTrend, mZSeas,cp);
ckd = size(mW0,2); 
vBeta=vBeta0(ckd+1:end);
ck=length(vBeta);
mZt=[mZt_2, mXreg(:, 1:end)];
mGG=mGG_2;
mT= blkdiag(mT_2, eye(ck,ck));
mHH=  blkdiag(mHH_2, zeros(ck,ck));
va = [va_2; vBeta];
mP= blkdiag(mP_2, zeros(ck,ck));

end