%% Monte Carlo experiment

clear all 
close all

% please set the experiment

%% Number of iterations MC experiment

MC=10;

%% Set the sample size of your time series

vn=[250 1000];

%% Set the trend as a Random Walk with variance

dsigma2_eps= 0.05  ;


%% Set the parameters of the fSW process and its frequency

vd=[0.35  0.45];  
dsigma2=0.7;
dLambda=2*pi/365.25;


%% Run the Monte Carlo experiment



opts = optimset('Display','Off','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter',500, 'MaxFunEvals', 500,...
    'LargeScale', 'off' );


iTrend=1;
cp=3;
cnd=length(vd);
cnn=length(vn);
mPar=NaN(cnd, 2+iTrend ,MC,vn(cnn));
mMSEi=NaN(MC,cnn);
mMSE =NaN(cnd,cnn);

vmu0=NaN(vn(cnn),1);
mX0 =[ones(vn(cnn),1) ]; %Set regression components
cCyc=length(dLambda);

for jk=1:cCyc
    mX0 = [mX0 , cos(dLambda(jk)*(1:vn(cnn))'), sin(dLambda(jk)*(1:vn(cnn))')] ;
end
    
for ii=1:cnd
    dd=vd(ii);
    for mc=1:MC

        for nn=1:cnn
            cn=vn(nn);
            
            [vs0]=fSim_fSW_typeII(dd,dLambda,dsigma2,vn(cnn));
            vs=vs0(1:cn);
            veps=sqrt(dsigma2_eps)*randn(vn(cnn)+5000,1);
            %e= randn(vn(cnn)+5000,1);
            
            for t=1:vn(cnn)+5000
                vmu0(t) = sum(veps(1:t));
            end
            vmu =  vmu0(5001:5000+cn)  ;
            vy = vmu + vs; 
            
            %figure();plot(1:cn,vy);
            
            mX=mX0(1:cn,:);
            [vP0,vLL]= fgrid_InitCondsBN(vy, dLambda , mX,cp  );

            fU  = @(vP)  fSsfLogLik_Harvey(vP, vy, dLambda, cp, mX); % function

            [vP_hat, dDeviance, exitflag, ~] = fminsearch(fU, vP0, opts);

            vd_hat          = exp(vP_hat(1:cCyc))./(1+exp(vP_hat(1:cCyc)));
            vSigma2_eta_hat = exp(2*vP_hat(cCyc+1:2*cCyc));
            dsigma2_eta0_hat = exp(2*vP_hat(1+2*cCyc)) ;

            mPar(ii, :,mc,nn)=[ vd_hat'  vSigma2_eta_hat'  dsigma2_eta0_hat   ]    ;

            % Exatract all the components
            [mZt, mGG, mT, mHH, va, mP ] = SetStateSpaceModel_Harvey(vP_hat ,dLambda,cn,cp);
            [mArt_smooth, mPrt_smooth, vBeta, mVarBeta, aP] = fAugmentedStateSmoother_Harvey(vy, mZt, ....
                mGG, mT, mHH, va, mP, [], mX);
            %SmoothEst=mArt_smooth([1,2:(cp+1):( 1+2*cCyc)*(cp)],:);
            vmu_hat=vBeta(1) + mArt_smooth(1,:)'  ;
            %figure();plot(1:cn,vmu);hold on;plot(1:cn,vmu_hat)
            mMSEi(mc,nn)=mean((vmu-vmu_hat).^2);
        end
        [ii mc]
    end
    mMSE(ii,:)=mean(mMSEi);
end


mMSE  





