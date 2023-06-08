%% APPLICATION


clear all
close all


%% DO YOU WANT BOX COX ? YES BC=1, NO BC=0
BC=0;

% econometrics toolbox required for autocorr
%% Graphic parameters
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
%% Load the GPS daily time series vy2
%replace "DRAO_NTAL_SR.txt" with the name of the file you want to process
LBpvalues=NaN(5,1); R2=NaN(5,1); SumCorr2=NaN(5,1);
mPar=NaN(6,5);

for jjj= 1

    if jjj==1; structData = readtable(['DRAO_NTAL_SR.txt']); end
 

mData = table2array(structData);
vtime0 =  datetime(mData(:,1),'convertfrom','modifiedjuliandate');
vy0    = mData(:,2) ;

vtime=(vtime0(1):vtime0(end))';
cn    = length(vtime);
vy2=NaN(cn,1); i=1;
for t=1:cn
    if vtime(t)==vtime0(i)
        vy2(t)= (vy0(i)); i=i+1;
    end
end
 
vr=vy2; 
 
 

%% Plot the series


g1 = figure('Name','Original Series');
subplot(2,1,1)
plot(vtime,  vr);
xlim([vtime(1)  vtime(end)])
%ylabel('Trading Volume', 'rotation', 90,'Interpreter','latex');
xlabel(' ', 'Interpreter','latex'); ylabel('',   'Interpreter','latex');
title( strcat('Series',string(jjj)), 'Fontweight', 'Normal', 'FontSize', 20,  'Interpreter','latex');
subplot(2,1,2)
[fxx,ww]=plomb(vr,[],2*pi*(0:cn-1)/cn);
plot(ww,fxx,'b');
%ylim([0 10]);
xlim([0, pi/6]);
xlabel(' ', 'Interpreter','latex'); ylabel('',   'Interpreter','latex');
title('GPS Periodogram ', 'Fontweight', 'Normal', 'FontSize', 20,  'Interpreter','latex');
orient(g1,'landscape')
% print(g1,strcat('Series',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')



%% Specification set periodicity and ARMA approximation of FN components

harm=2;
vPeriod  = [365.25];
dLambda  = (2*pi./vPeriod)*(1:harm); 
cCyc=length(dLambda);

cp=3; %ARMA(cp,cp) approximation of each FN components (see Maddanu Proietti 2022, Trends in ethane)

mX =[ones(cn,1) (1:cn)']; %Set regression components
for jk=1:cCyc
    mX = [mX , cos(dLambda(jk)*(1:cn)'), sin(dLambda(jk)*(1:cn)')] ;
end


%% BOX COX TRANSFORM


   if BC==1
   [v1,v2]= f_fSW_BoxCox(vr, dLambda, cp, mX,-1.1*min(vr));
   %v1=vVV(jjj);
   if v1<0.001
       vy=log(vr);
   else
       vy=((vr+v2).^v1 - 1)/v1;
   end
   else
       v1=1;
       vy=vr;
   end




%% Estimation
 

%Initial consitions
%[vP0,vLL]= fgrid_InitConds(vy, dLambda , mX,cp  );
%[vP0,vLL]= fgrid_InitCondsLC(vy , mX, dLambda,cp ) ;
[vP0,vLL]= fgrid_InitCondsRN(vy, dLambda , mX,cp  );

% maximization of the log likelihood

fU  = @(vP)  fSsfLogLik_Harvey(vP, vy, dLambda, cp, mX); % function
opts = optimset('Display','Iter','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter',500, 'MaxFunEvals', 500,...
    'LargeScale', 'off' );
[vP_hat, dDeviance, exitflag, ~] = fminsearch(fU, vP0, opts);
 
vd_hat          = exp(vP_hat(1:cCyc))./(1+exp(vP_hat(1:cCyc)))
vSigma2_eta_hat = exp(2*vP_hat(cCyc+1:2*cCyc))
dsigma2_eta0_hat = exp(2*vP_hat(1+2*cCyc)) 
dphi_hat =  vP_hat(2*cCyc+2)/sqrt(1 + vP_hat(2*cCyc+2)^2)
dpsi_hat =  vP_hat(2*cCyc+3)/sqrt(1 + vP_hat(2*cCyc+3)^2)

mPar(:,jjj)=[ vd_hat'  vSigma2_eta_hat'  dsigma2_eta0_hat dphi_hat   ]    ;

% Exatract all the components 
[mZt, mGG, mT, mHH, va, mP ] = SetStateSpaceModel_Harvey(vP_hat ,dLambda,cn,cp);
[mArt_smooth, mPrt_smooth, vBeta, mVarBeta, aP] = fAugmentedStateSmoother_Harvey(vy, mZt, ....
    mGG, mT, mHH, va, mP, [], mX);

[mZts, mGGs, mTs, mHHs, vas, mPs ] = SetStateSpaceModel_ExactLik(vP_hat,dLambda,mX,vBeta,cn,cp);
[vRes, vVarRes, mStatePred, mCovStatePred,dLogLik_exact ] = KalmanFilter_ExactLik(vy,mZts , mGGs, mTs, mHHs, vas, mPs);


% SmoothEst=mArt_smooth([1:(cp+1):( 2*cCyc)*(cp+1)],:);
% SmoothEst_var=mPrt_smooth([1:(cp+1):( 2*cCyc)*(cp+1)],:);
SmoothEst=mArt_smooth([1,3:(cp+1):( 1+2*cCyc)*(cp)],:);
SmoothEst_var=mPrt_smooth([1,3:(cp+1):( 1+2*cCyc)*(cp)],:);

vy_hat=mX*vBeta +  SmoothEst(1,:)' + sum(mX(:,[3:end]).* SmoothEst(2:end,:)',2);
 
g3=figure()
plot(vtime,vy,'b')
hold on
plot(vtime,vy_hat,'r--','LineWidth',1.7)
xlim([vtime(1) vtime(end)])
title('Kalman Filter smooth fit of the series')
orient(g3,'landscape')
%print(g3,'SeriesFit','-dpdf', '-fillpage', '-painters','-r600')

 
if BC==1
    [vTrend,v1stCyc,v2ndCyc,vTrend_var,v1stCyc_var,v2ndCyc_var]=fSimSmoother(vy,vP_hat,dLambda,mX,vBeta,cp,v1,v2,1000);
end

g1=figure()
plot(vtime,  vr  ,'b:')
hold on
if BC==1
    if v1>0.01
        plot(vtime, (1 +v1*mX(:,1:2)*vBeta(1:2)).^(1/v1) - v2 ,'r','LineWidth',2)
    else
        plot(vtime, exp( mX(:,1:2)*vBeta(1:2)) -v2 ,'r','LineWidth',2)
    end
else
    plot(vtime,  ( mX(:,1:2)*vBeta(1:2))  ,'r','LineWidth',2)
end
xlim([vtime(1) vtime(end)])
title('Deterministic Trend component')
orient(g1,'landscape')
print(g1,strcat('Trend',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')

g1=figure()
plot(vtime, SmoothEst(1,:),'r','LineWidth',2)
xlim([vtime(1) vtime(end)])
title('Arma(1,1) component')
orient(g1,'landscape')

mX2=mX(:,3:end);
SmoothEst2=SmoothEst(2:end,:);
SmoothEst_var2=SmoothEst_var(2:end,:);
vBeta2=vBeta(3:end);
mVarBeta2=diag(mVarBeta(3:end,3:end));
vCycVar=(mX2(1).^2).*mVarBeta2(1)+ (mX2(2).^2).*mVarBeta2(2) + sum((mX2(:,1:2).^2).* (SmoothEst_var2(1:2,:))',2) ;
% % mX2=mX(:,2:end);
% % SmoothEst2=SmoothEst(2:end,:);
% % vBeta2=vBeta(2:end);


g1=figure()
if BC==1
    plot(vtime, vr ,'b:')
    hold on
    plot(vtime, v1stCyc ,'r','LineWidth',2) 
    else

    plot(vtime,  - mX2(:,1:2)*vBeta2(1:2) +  sum(mX2(:,1:2).* SmoothEst2(1:2,:)',2) ,'r')
end
title('1st Cyclical component')
xlim([vtime(1) vtime(end)])
orient(g1,'landscape') 
print(g1,strcat('CyC',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')

if harm>1
    g1=figure()
    if BC==1
        plot(vtime, vr ,'b:')
        hold on
        plot(vtime, v2ndCyc ,'r','LineWidth',2)
    else
        plot(vtime, - mX2(:,3:4)*vBeta2(3:4) +  sum(mX2(:,3:4).* SmoothEst2(3:4,:)',2) ,'r')
    end
    title('2nd Cyclical component')
    xlim([vtime(1) vtime(end)])
    orient(g1,'landscape')
    print(g1,strcat('2CyC',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')
end

g1=figure()
subplot(2,1,1)
plot(vtime,sqrt((SmoothEst2(1,:)-vBeta2(1)).^2  + (SmoothEst2(2,:)-vBeta2(2)).^2) ,'r')
title('Amplitude 1st Cyclical component')
xlim([vtime(1) vtime(end)])
subplot(2,1,2)
plot(vtime,atan((SmoothEst2(2,:)-vBeta2(2))./(SmoothEst2(1,:)-vBeta2(1))) ,'r')
title('Phase 1st Cyclical component')
xlim([vtime(1) vtime(end)])
orient(g1,'landscape') 
print(g1,strcat('Amp',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')

if harm>1
g1=figure()
subplot(2,1,1)
plot(vtime,sqrt((SmoothEst2(3,:)-vBeta2(3)).^2  + (SmoothEst2(4,:)-vBeta2(4)).^2) ,'r')
title('Amplitude 2nd Cyclical component')
xlim([vtime(1) vtime(end)])
subplot(2,1,2)
plot(vtime,atan((SmoothEst2(4,:)-vBeta2(4))./(SmoothEst2(3,:)-vBeta2(3))) ,'r')
title('Phase 2nd Cyclical component')
xlim([vtime(1) vtime(end)])
orient(g1,'landscape') 
print(g1,strcat('2Amp',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')

end

g4 = figure('Name','res');
subplot(3,1,1);
plot(vtime, vRes./sqrt(vVarRes))
xlim([vtime(1) vtime(end)])
title('Residuals: time series plot, periodogram, correlogram')
subplot(3,1,2);
[fxxRES,wwRES]=plomb(vRes./sqrt(vVarRes) ,[],2*pi*(0:cn-1)/cn);
plot(wwRES,fxxRES,'b');xlim([0, 0.1]);title('Lomb-Scargle periodogram','Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex')
%[ fxxRES, wwRES ] = fPeriodogram(vRes./sqrt(vVarRes));
xlabel('$\omega$', 'Interpreter','latex');
ylabel('$I_{\epsilon}(\omega)$',   'Interpreter','latex');
xlim([0, pi ]);
subplot(3,1,3);
Lag=120;
acfR =NaN(Lag+1,1);
for k=0:Lag
    acfR (k+1)= real( (1/cn) * sum(fxxRES(1:end).*exp(1i*k *wwRES(1:end)))) ;
end
SE=real(sqrt((1/cn)*(1+2*sum(real(acfR(1:51)/acfR(1)).^2))));
% % [acf,lags,bounds] = autocorr(vRes./sqrt(vVarRes),cn-1 );
% % bar(0:50, acf(1:51),'r')
bar(0:50,acfR(1:51)/acfR(1),'r')
% yline(bounds(1)  , 'b--')
% yline( bounds(2)  , 'b--')
yline(2*SE , 'b--')
yline(-2*SE , 'b--')
title('Sample autocorrelation', 'Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
xlabel('k','Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
ylabel('$\hat{\rho}_\epsilon(k)$','Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
orient(g4,'landscape')
print(g4,strcat('Res',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')
[h,pValue,stat,cValue] = lbqtest(vRes./sqrt(vVarRes));
LBpvalues(jjj)=pValue;
R2(jjj)  = 1 - sum(vRes.^2)/sum((vy-mean(vy)).^2);
SumCorr2(jjj) = sum((acfR(2:end)/acfR(1)).^2);
jjj
%close all
end


 
