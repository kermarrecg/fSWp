function [mZt, mGG, mT, mHH, va, mP,mH1] = SetStateSpaceModel_Harvey(vP,dLambda,cn,cp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ck=length(dLambda);
 
% Assign hyperparameters  %%vP=randn(2*k+2,1)
dd = exp(vP(1:ck))./(1+exp(vP(1:ck))) ;
dsigma  = exp(2*vP(ck+1:2*ck)).^0.5;
dsigma_eta = exp(2*vP(1+2*ck)).^0.5;
%dd_trend=exp(vP(2*ck+2))./(1+exp(vP(2*ck+2))) ;
dphi =  vP(2*ck+2)/sqrt(1 + vP(2*ck+2)^2);
%dsigma_zeta=exp(2*vP(3+2*ck)).^0.5;
%-----------------------------------------------------------------------


opts = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter', 1000, 'MaxFunEvals', 1000,...
    'LargeScale', 'off' );

 
vP10 = zeros(2*cp,1);
% 
%cm = max(cp,cp+1);     % number of state elements
% mZ0 = [1, zeros(1, cm-1)];    
% mT = zeros(cm,cm); mT(1,1) = 1;
% mT(:,2:end) = [eye(cm-1); zeros(1, cm-1)];%  mT(1,1)=dphi;
% %mT(2,2) = dphi;
% mH = [dsigma_eta ;zeros(cm-1,1)]; mH(2:2+cp-1) = [0; zeros(cp-1,1)];
% %mH = [dsigma_eta; dsigma_zeta ;zeros(cm-2,1)]; mH(3:2+cp-1) = [ zeros(cp-1,1)];
% mHH = (mH*mH'); % H_t * H_t'  
% va = zeros(cm,1);     
% mP = mHH;

% mZ0=[1 , 1];
% mT  =[1 0; 0 dphi];     mHH = diag([dsigma_eta^2, dsigma_zeta^2]);
% va = zeros(2,1);   
mZ0=1;
mT  = dphi ;     mHH =  dsigma_eta^2 ; mH1= dsigma_eta;
va =0;  
mP = mHH;

% f    = @(vP1)fARMA_Approx(dd_trend, cn, cp, vP1); % function
% [vP1, ~, ~, ~] = fminunc(f, vP10, opts );
% vtAR    = fFisherInvTransform(vP1(1:cp));
% vtMA    = fFisherInvTransform(vP1(cp+1:2*cp));
% [mZ0, mT0, mHH0, va0, ~] = fSsfARMA(fReparAR(vtAR), fReparMA(vtMA), dsigma_eta^2);
%mT = mT0; mHH = mHH0; mP = mHH0; va = va0;
%mT = []; mHH = []; mP = []; va = [];

for jd=1:ck

    f    = @(vP1)fARMA_Approx(dd(jd), cn, cp, vP1); % function
    [vP1, ~, ~, ~] = fminunc(f, vP10, opts );
    vtAR    = fFisherInvTransform(vP1(1:cp));
    vtMA    = fFisherInvTransform(vP1(cp+1:2*cp));
    %[vphi, vtheta]
    [~, mTs, mHHs, vas,~, mH] = fSsfARMA(fReparAR(vtAR), fReparMA(vtMA), dsigma(jd)^2);
    mH1=[mH1; mH;mH];
    mT  = blkdiag(mT, kron(eye(2),mTs));
    mHH = blkdiag(mHH, kron(eye(2), mHHs));
    %mP = blkdiag(mP, kron(eye(2), mHHs));     
    mP = blkdiag(mP, kron(eye(2), mHHs));

    va = [va; vas;vas];
    
end

mGG = 0;% dsigma_eta^2 ;%exp(2*vP(3+2*ck)).^0.5;
vt    = linspace(1,cn,cn)';
mZSeas = zeros(cn, 2*ck*(cp+1));
for k = 1:ck
    mZSeas(:,2*(k-1)*(cp+1)+1) = cos(dLambda(k)*vt);
    mZSeas(:,(2*k-1)*(cp+1)+1) = sin(dLambda(k)*vt);
end
mZt = [ones(cn,1)*mZ0, mZSeas];
%mZt = mZSeas; 
end
