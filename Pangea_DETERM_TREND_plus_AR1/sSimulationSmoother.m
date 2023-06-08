
MC=1000; v1=0;

cT=length(vy);
ck=length(dLambda);


[mZt, mGG, mT, mHH, va, mP,mH] = SetStateSpaceModel_Harvey(vP_hat,dLambda,cT,cp);
[mAlpha_hat, ~, ~, ~, ~] = fAugmentedStateSmoother_Harvey(vy, mZt, mGG, mT, mHH, va, mP, [], mX);

vpseudoy= zeros( cT,1);
vTrend_MC  = zeros( cT,1);
v1stCyc_MC = zeros( cT,1);
v2ndCyc_MC = zeros( cT,1);
v1stAmp_MC = zeros( cT,1);
v1stPh_MC  = zeros( cT,1);
v2ndAmp_MC = zeros( cT,1);
v2ndPh_MC  = zeros( cT,1);

for mc=1:MC
    vy_plus=NaN(cT,1);
    mAlpha_plus= NaN(1+2*ck*(cp+1),cT);
    mAlpha_plus(:,1)= mAlpha_hat(:,1);
    
    mEta=randn(1,cT);
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
    SmoothEst=mAlpha_tilde([1,2:(cp+1):( 1+2*cCyc)*(cp)],:);
 
    mX2=mX(:,3:end);
    SmoothEst2=SmoothEst(2:end,:);
    vBeta2=vBeta(3:end);
if v1>0.001
    vTrend_MC= vTrend_MC       + (1 + v1*(mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)')).^(1/v1);
    v1stCyc_MC =  v1stCyc_MC   + (1 + v1*( - mX2(:,1:2)*vBeta2(1:2) +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) )).^(1/v1);
     vpseudoy  =   vpseudoy  + (1 + v1*( mX*vBeta + SmoothEst(1,:)'+ sum(mX(:,3:end).* SmoothEst(2:end,:)',2))).^(1/v1);
    v2ndCyc_MC =  v2ndCyc_MC   + (1 + v1*( - mX2(:,3:4)*vBeta2(3:4) +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2))).^(1/v1);
    v1stAmp_MC =  v1stAmp_MC   + (1 +  v1*sqrt((SmoothEst2(1,:)-vBeta2(1)).^2  + (SmoothEst2(2,:)-vBeta2(2)).^2) ).^(1/v1);
    v1stPh_MC =  v1stPh_MC     + (1 +  v1*(atan((SmoothEst2(2,:)-vBeta2(2))./(SmoothEst2(1,:)-vBeta2(1)))) ).^(1/v1);
    v2ndAmp_MC =  v2ndAmp_MC   + (1 +  v1*(sqrt((SmoothEst2(3,:)-vBeta2(3)).^2  + (SmoothEst2(4,:)-vBeta2(4)).^2)) ).^(1/v1);
    v2ndPh_MC =  v2ndPh_MC     + (1 +  v1*(atan((SmoothEst2(4,:)-vBeta2(4))./(SmoothEst2(3,:)-vBeta2(3)))) ).^(1/v1);
else
    vTrend_MC=    vTrend_MC    + exp( mX(:,1:2)*vBeta(1:2) + SmoothEst(1,:)') ;
    v1stCyc_MC =  v1stCyc_MC   + exp( ( - mX2(:,1:2)*vBeta2(1:2) +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) )) ;
    vpseudoy  =   vpseudoy  +   exp( mX*vBeta + SmoothEst(1,:)'+ sum(mX(:,3:end).* SmoothEst(2:end,:)',2)) ;                   
    v2ndCyc_MC =  v2ndCyc_MC   + exp( ( -mX2(:,3:4)*vBeta2(3:4) +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2))) ;
    v1stAmp_MC =  v1stAmp_MC   + exp(  sqrt((SmoothEst2(1,:)-vBeta2(1)).^2  + (SmoothEst2(2,:)-vBeta2(2)).^2) )' ;
    v1stPh_MC =   v1stPh_MC    + exp( (atan((SmoothEst2(2,:)-vBeta2(2))./(SmoothEst2(1,:)-vBeta2(1)))) )';
    v2ndAmp_MC =  v2ndAmp_MC   + exp(  (sqrt((SmoothEst2(3,:)-vBeta2(3)).^2  + (SmoothEst2(4,:)-vBeta2(4)).^2)) )' ;
    v2ndPh_MC =   v2ndPh_MC    + exp( (atan((SmoothEst2(4,:)-vBeta2(4))./(SmoothEst2(3,:)-vBeta2(3)))) )' ;
end

mc
end

vTrend =vTrend_MC./MC;
v1stCyc =v1stCyc_MC./MC;
v2ndCyc =v2ndCyc_MC./MC;
v1stAmp =v1stAmp_MC./MC;
v1stPh =v1stPh_MC./MC;
v2ndAmp =v2ndAmp_MC./MC;
v2ndPh =v2ndPh_MC./MC;

figure()
% plot(vtime,vr,'b')
% hold on
plot(vtime,vr,'b')
hold on
plot(vtime,vTrend ,'r')

figure()
plot(vtime,vr,'b')
hold on
plot(vtime,vpseudoy/MC,'r--')

figure()
plot(vtime,vr,'b')
hold on
plot(vtime,vTrend +v1stCyc + v2ndCyc,'r--')


figure()
plot(vtime,v1stCyc)
figure()
plot(vtime,v2ndCyc)

figure()
plot(vtime,v1stAmp)
figure()
plot(vtime,v2ndAmp)

