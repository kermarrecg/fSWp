%% Model Specification

%y_t= nu_t + s_t + u_t
%mu_t= b_0 + b_1 t
%s_t= a_1 cos(wt) + a*_1 sin(wt) + a_2 cos(2wt) + a*_2 sin(2wt) 
%u_t= (1-L)^d e_t
%e_t is a WN(0,sigma^2)

close all
clear all

%% DO YOU WANT BOX COX ? YES BC=1, NO BC=0
BC=0;

%% Graphic parameters
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
%% Load the GPS daily time series vy2
LBpvalues=NaN(10,1); R2=NaN(10,1); SumCorr2=NaN(10,1);
mPar=NaN(8,1);Miss=NaN(10,1); vV=NaN(8,1);

 

for jjj=  1 : 4

 
    if jjj==1; structData = readtable(['DRAO_PWV.txt']); ID=  'DRAO';end
    if jjj==2; structData = readtable(['DRAO_NTAL_SR.txt']); ID=  'DRAO';end
    if jjj==3; structData = readtable(['DRAO_HYDL_SR.txt']); ID=  'DRAO';end
    if jjj==4; structData = readtable(['DRAO_IGS.txt']); ID=  'DRAO';end
    
    if jjj==1
        vtime0 = datetime(structData.Date_Time,"InputFormat",'yyyy-MM-dd''T''HH:mm');
        vy0    =structData.PWV_mm_;
        step=duration('06:00:00');
        vtime1=vtime0(1):step:vtime0(end);
        cn0    = length(vtime1);
        vy2=NaN(cn0,1); i=1;
        for t=1:cn0
            if vtime1(t)==vtime0(i)
                vy2(t)= (vy0(i)); i=i+1;
            end
        end

        TT=timetable(vtime1',vy2);
        TTdailyMean = retime(TT,'daily','mean');

        vy0=(TTdailyMean.vy2);

        vtime0=TTdailyMean.Time;

    else
        mData = table2array(structData);
        vtime0 =  datetime(mData(:,1),'convertfrom','modifiedjuliandate');
        vy0    = mData(:,2) ;
   end

vtime=(vtime0(1):vtime0(end))';
cn    = length(vtime);
vy2=NaN(cn,1); i=1;
for t=1:cn
    if vtime(t)==vtime0(i)
        vy2(t)= (vy0(i)); i=i+1;
    end
end
 
vr=vy2; 
 

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
   
 
    if BC==1
        vr(vr<=0)=NaN; sum(vr<=0)
        [v1,v2]= f_PinkNoise_BoxCox(vr,   mX,cp);
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

   vV(jjj)=v1;
 

    Miss(jjj)=sum(isnan(vy))/cn


    %% Plot the series


    g1 = figure('Name','Original Series');
    subplot(2,1,1)
    plot(vtime,  vr);
    xlim([vtime(1)  vtime(end)])
    %ylabel('Trading Volume', 'rotation', 90,'Interpreter','latex');
    xlabel(' ', 'Interpreter','latex'); ylabel('',   'Interpreter','latex');
    title( strcat(ID,' Series ',string(jjj)), 'Fontweight', 'Normal', 'FontSize', 20,  'Interpreter','latex');
    subplot(2,1,2)
    [fxx,ww]=plomb(vr,[],2*pi*(0:cn-1)/cn);
    plot(ww,fxx,'b');
    % [ fxx, ww ] = fPeriodogram(vy);
    %ylim([0 10]);
    xlim([0, pi/6]);
    xlabel(' ', 'Interpreter','latex'); ylabel('',   'Interpreter','latex');
    title('GPS Periodogram ', 'Fontweight', 'Normal', 'FontSize', 20,  'Interpreter','latex');
    orient(g1,'landscape')
    print(g1,strcat(ID,' Series ',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')



    %% Estimation


    %Initial conditions
%    dd0=exp(0.3)./(1+exp(0.3));
% %  vP0=[ log( nanvar(vy)*((gamma(1-dd0)^2)/(gamma(1-2*dd0))))*0.5 ;0.3  ];
%     vP0=[ log( nanvar(vy) )*0.5 ; 0.3  ];
[vP0,vLL]= fgrid_InitCondsPN(vy , mX,cp  );

    % maximization of the log likelihood

    fU  = @(vP)  fSsfLogLik_Harvey(vP, vy, mX,cp); % function
    opts = optimset('Display','Iter','TolX',1e-4,'TolFun',1e-4,...
        'Diagnostics','off', 'MaxIter',1500, 'MaxFunEvals', 1500,...
        'LargeScale', 'off' );
    [vP_hat, dDeviance, exitflag, ~] = fminsearch(fU, vP0, opts);

 
    dsigma_eta_hat = exp(2*vP_hat(1 )).^0.5;
    dd_hat =  exp(vP_hat(2))./(1+exp(vP_hat(2)));
 

    % Exatract all the components
    [mZt, mGG, mT, mHH, va, mP ] = SetStateSpaceModel_Harvey(vP_hat ,cn,cp);
    [mArt_smooth, mPrt_smooth, vBeta, mVarBeta, aP] = fAugmentedStateSmoother_Harvey(vy, mZt, ....
        mGG, mT, mHH, va, mP, [], mX);

    mPar(:,jjj)=[ vBeta'   dsigma_eta_hat dd_hat]    ;

    [mZts, mGGs, mTs, mHHs, vas, mPs ] = SetStateSpaceModel_ExactLik(vP_hat,mX,vBeta,cn,cp);
    [vRes, vVarRes, mStatePred, mCovStatePred,dLogLik_exact ] = KalmanFilter_ExactLik(vy,mZts , mGGs, mTs, mHHs, vas, mPs);


  SmoothEst=mArt_smooth  ;
  SmoothEst_var=mPrt_smooth ;
  vy_hat=mX*vBeta + SmoothEst(1,:)' ;
%     vy_hat=mX*vBeta +sum(SmoothEst(1:2,:)',2) + sum(mX(:,[2:end]).* SmoothEst(3:end,:)',2);

    g3=figure()
    plot(vtime,vy,'b')
    hold on
    plot(vtime,vy_hat,'r--','LineWidth',1.7)
    xlim([vtime(1) vtime(end)])
    title('Kalman Filter smooth fit of the series')
    orient(g3,'landscape')
    %print(g3,'SeriesFit','-dpdf', '-fillpage', '-painters','-r600')

    figure()
    plot(vtime,  SmoothEst(1,:)'  ,'r','LineWidth',2)
    
    if BC==1
        if v1>0.001
            v1stCyc =   (1 + v1* mX(:,1:4)*vBeta(1:4) ).^(1/v1);
            v2ndCyc =   (1 + v1*( mX(:,[1,2,5,6])*vBeta([1,2,5,6]) )).^(1/v1);
        else
            v1stCyc  =    exp( (mX(:,1:4)*vBeta(1:4)  )) ;
            v2ndCyc  =   exp( (mX(:,[1,2,5,6])*vBeta([1,2,5,6])  )) ;

        end
    else 
        v1stCyc  =     ( (mX(:,1:4)*vBeta(1:4)  )) ;
        v2ndCyc  =    ( (mX(:,[1,2,5,6])*vBeta([1,2,5,6])  )) ;
    end

    g1=figure()
    plot(vtime,  vr  ,'b:')
    hold on
%   
    if BC==1
    
    
    if v1>0.01
        plot(vtime, (1 +v1*mX(:,1:2)*vBeta(1:2)).^(1/v1) ,'r','LineWidth',2)
    else
        plot(vtime, exp( mX(:,1:2)*vBeta(1:2))  ,'r','LineWidth',2)
    end
    else
        plot(vtime, ( mX(:,1:2)*vBeta(1:2))  ,'r','LineWidth',2)
    end
    xlim([vtime(1) vtime(end)])
    title('Trend component of the original series')
    orient(g1,'landscape')
    print(g1,strcat('ID','_Trend_',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')
    

    g1=figure()
    plot(vtime,  SmoothEst(1,:)'  ,'r','LineWidth',1.2)
    xlim([vtime(1) vtime(end)])
    title('Pink Noise component ')
    orient(g1,'landscape')
    % print(g1,strcat('ID','_PinkNoise_',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')


    g1=figure()
    plot(vtime, vr ,'b:')
    hold on
    plot(vtime, v1stCyc ,'r','LineWidth',2) 
    title('1st Cyclical component of the orignal series')
    xlim([vtime(1) vtime(end)])
    orient(g1,'landscape')
    print(g1,strcat('ID','_CyC_',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')

    if harm>1
        g1=figure()
        plot(vtime, vr ,'b:')
        hold on
        plot(vtime, v2ndCyc ,'r','LineWidth',2)
        hold on
        title('2nd Cyclical component of the original series')
        xlim([vtime(1) vtime(end)])
        orient(g1,'landscape')
        print(g1,strcat('ID','_2ndCyC_',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')
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
    SE=real(sqrt((1/cn)*(1+2*sum(real(acfR(1:51)/acfR(1)).^2))));bounds=SE*[-1 2];
    %[acf,lags,bounds] = autocorr(vRes./sqrt(vVarRes),cn-1 );
    % bar(0:50, acf(1:51),'r')
    bar(0:50,acfR(1:51)/acfR(1),'r')
    yline(bounds(1)  , 'b--')
    yline( bounds(2)  , 'b--')
    title('Sample autocorrelation', 'Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
    xlabel('k','Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
    ylabel('$\hat{\rho}_\epsilon(k)$','Fontweight', 'Normal', 'FontSize', 15,  'Interpreter','latex');
    orient(g4,'landscape')
    print(g4,strcat('ID','_Res_',string(jjj)),'-dpdf', '-fillpage', '-painters','-r600')
    [h,pValue,stat,cValue] = lbqtest(vRes./sqrt(vVarRes));
    LBpvalues(jjj)=pValue;
    R2(jjj)  = 1 - sum(vRes.^2)/sum((vy-mean(vy)).^2);
    SumCorr2(jjj) = sum((acfR(2:end)/acfR(1)).^2);
    jjj
    disp(['The uncertency of the trend is +-',string(sqrt(mVarBeta(2,2)))])
    %close all
end
 

 