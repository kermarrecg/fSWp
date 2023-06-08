function [mZt, mGG, mT, mHH, va, mP ] = SetStateSpaceModel_ExactLik(vP,dLambda,mX_reg,vBeta,cn,cp)

mk=length(dLambda);
[mZt_2, mGG_2, mT_2, mHH_2, va_2, mP_2 ] = SetStateSpaceModel_Harvey(vP,dLambda,cn,cp);
 
ck=length(vBeta);
mZt=[mZt_2, mX_reg];
mGG=mGG_2;
mT= blkdiag(mT_2, eye(ck,ck));
mHH=  blkdiag(mHH_2, zeros(ck,ck));
va = [va_2; vBeta];
mP= blkdiag(mP_2, zeros(ck,ck));

end