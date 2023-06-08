function [mZt, mGG, mT, mHH, mW, va, mP, mW0] = fSetSSModel(vP, iTrend, mZSeas,cp)
% iTrend = 1 for LLM iTrend = 2 for LLTM
%% Hyperparameters
dSeta  = exp(2*vP(1));
if iTrend == 2
    dSzeta = exp(2*vP(2));
end
[cn, cms] = size(mZSeas);
cs        = ceil(cms/2);
vd        = exp(vP(3:2+cs))./(1+exp(vP(3:2+cs))); % memory parameters
vSomega   =  exp(2*vP(3+cs:2+2*cs)) ;   % variance parameters
%% measurement equation
if iTrend == 2
    mZt =  (ones(cn,1) * [1, 0, 1]) ;
else  mZt =  ones(cn,1) ;
end
mGG =0;dSAR = exp(2*vP(end-1)); dphi=  vP(end)/sqrt(1 + vP(end)^2);%0; % G_t * G_t'
%% transition equation
if iTrend == 2
    mT  =[1 1 0; 0 1 0; 0 0 dphi];     mHH = diag([dSeta, dSzeta, dSAR]);
else mT = 1; mHH = dSeta; 
end
%cp  = 2; % Hartl MA and AR order
vp0 = zeros(2*cp,1);
 
opts = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter', 1000, 'MaxFunEvals', 1000,...
    'LargeScale', 'off' );

for jd = 1:cs
    f    = @(vp1)fARMA_Approx(vd(jd), cn, cp, vp1); % function
    [vp1, ~, ~, ~] = fminunc(f, vp0,opts );
    vtAR    = fFisherInvTransform(vp1(1:cp));
    vphi   =  fReparAR(vtAR) ;
    vtMA    = fFisherInvTransform(vp1(cp+1:2*cp));
    vtheta  = fReparMA(vtMA) ;
    [mZarma, mTarma, mHHarma, ~, ~,~] = fSsfARMA(vphi, vtheta,vSomega(jd));
    mZt  = [mZt,  mZSeas(:, 2*jd-1) * mZarma, mZSeas(:, 2*jd ) * mZarma];
    mT   = blkdiag(mT, kron(eye(2),mTarma));
    mHH  = blkdiag(mHH,kron(eye(2),mHHarma));
    vp0 = vp1;
end
cm  = size(mT,2);
mW  = 0;  
%% Initial conditions 
va  = zeros(cm, 1); 
mP  = mHH;
if iTrend == 2
    mW0 = mT(:,1:3);
else mW0 = [1; zeros(cm-1,1)];
end
end