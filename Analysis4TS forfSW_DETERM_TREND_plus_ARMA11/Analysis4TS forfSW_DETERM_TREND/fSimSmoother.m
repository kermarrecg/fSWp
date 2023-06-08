function [vTrend,v1stCyc,v2ndCyc,vTrend_var,v1stCyc_var,v2ndCyc_var]=fSimSmoother(vy,vP_hat,dLambda,mX,vBeta,cp,v1,v2,MC)
 

cT=length(vy);
ck=length(dLambda);
 

[mZt, mGG, mT, mHH, va, mP,mH] = SetStateSpaceModel_Harvey(vP_hat,dLambda,cT,cp);
[mAlpha_hat, ~, ~, ~, ~] = fAugmentedStateSmoother_Harvey(vy, mZt, mGG, mT, mHH, va, mP, [], mX);

 
vTrend_MC  = zeros( cT,1);
v1stCyc_MC = zeros( cT,1);
v2ndCyc_MC = zeros( cT,1);
vTrend_MC2  = zeros( cT,1);
v1stCyc_MC2 = zeros( cT,1);
v2ndCyc_MC2 = zeros( cT,1); 

for mc=1:MC
    vy_plus=NaN(cT,1);
    mAlpha_plus= NaN(size(mAlpha_hat,1),cT);
    mAlpha_plus(:,1)= mAlpha_hat(:,1);
    
    mEta=randn(2,cT);
    mEta0=randn(2*ck,cT+cp+1);
    for k=1:2*ck
        for i=1:cp+1
            mEta=[mEta;mEta0(k,cp-i+2:cT+cp-i+1)];
        end
    end
    for t=1:cT
        %mEta= kron(randn(1+2*ck,1),ones(cp+1,1));
        vy_plus(t)=mX(t,:)*vBeta + mZt(t,:)* mAlpha_plus(:,t);
        if t==cT;break;end
        mAlpha_plus(:,t+1)=mT* mAlpha_plus(:,t)+mH.*mEta(:,t);
    end

%     figure();plot(vtime,vy_plus)
%     figure();plot(vtime,vy )

    [mAlpha_plushat, ~, ~, ~, ~] = fAugmentedStateSmoother_Harvey(vy_plus, mZt, mGG, mT, mHH, va, mP, [], mX);

    mAlpha_tilde=mAlpha_hat-mAlpha_plushat+mAlpha_plus;
    SmoothEst=mAlpha_tilde([1,2:(cp+1):( 1+2*ck)*(cp)],:);
 
    mX2=mX(:,3:end);
    SmoothEst2=SmoothEst(2:end,:);
    vBeta2=vBeta(3:end);
if v1>0.001
    vTrend_MC= vTrend_MC       + (1 + v1*(mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)')).^(1/v1) - v2;
    vTrend_MC2= vTrend_MC2       + (1 + v1*(mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)')).^(2/v1) - v2;
    v1stCyc_MC =  v1stCyc_MC   + (1 + v1*(  mX(:,1:4)*vBeta(1:4) +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) )).^(1/v1) - v2;
    v1stCyc_MC2 =  v1stCyc_MC2   + (1 + v1*( mX(:,1:4)*vBeta(1:4) +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) )).^(2/v1) - v2;
    v2ndCyc_MC =  v2ndCyc_MC   + (1 + v1*( mX(:,[1,2,5,6])*vBeta([1,2,5,6])  +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2))).^(1/v1) - v2;
    v2ndCyc_MC2 =  v2ndCyc_MC2   + (1 + v1*( mX(:,[1,2,5,6])*vBeta([1,2,5,6])   +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2))).^(2/v1) - v2; 
else
    vTrend_MC=    vTrend_MC    + exp( mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)')  - v2;
    vTrend_MC2=    vTrend_MC2    + exp( mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)').^2  - v2;
    v1stCyc_MC =  v1stCyc_MC   + exp( (mX(:,1:4)*vBeta(1:4)  +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) ))  - v2;
    v1stCyc_MC2 =  v1stCyc_MC2   + exp( ( mX(:,1:4)*vBeta(1:4)  +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) )).^2  - v2;
    v2ndCyc_MC =  v2ndCyc_MC   + exp( (mX(:,[1,2,5,6])*vBeta([1,2,5,6])  +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2)))  - v2;
    v2ndCyc_MC2 =  v2ndCyc_MC2   + exp( ( mX(:,[1,2,5,6])*vBeta([1,2,5,6])  +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2))).^2  - v2;
end

mc
end
vTrend =vTrend_MC./MC;
v1stCyc =v1stCyc_MC./MC;
v2ndCyc =v2ndCyc_MC./MC;
vTrend_var=  vTrend_MC2/MC   - vTrend.^2;
v1stCyc_var=  v1stCyc_MC2/MC  -v1stCyc.^2;
v2ndCyc_var=  v2ndCyc_MC2/MC - v2ndCyc.^2;

end