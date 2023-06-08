function [mZt, mGG, mT, mHH, va, mP] = SetStateSpaceModel_Harvey(vP,dLambda,cn,cp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ck=length(dLambda);
 
% Assign hyperparameters  %%vP=randn(2*k+2,1)
dd = exp(vP(1:ck))./(1+exp(vP(1:ck))) ;
dsigma  = exp(2*vP(ck+1:2*ck)).^0.5;
dsigma_eta = exp(2*vP(1+2*ck)).^0.5;
%dd_trend=exp(vP(2*ck+2))./(1+exp(vP(2*ck+2))) ;
dphi =  vP(2*ck+2)/sqrt(1 + vP(2*ck+2)^2);
dtheta =  vP(2*ck+3)/sqrt(1 + vP(2*ck+3)^2);
%dsigma_zeta=exp(2*vP(3+2*ck)).^0.5;
%-----------------------------------------------------------------------


opts = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,...
    'Diagnostics','off', 'MaxIter', 1000, 'MaxFunEvals', 1000,...
    'LargeScale', 'off' );

 
vP10 = zeros(2*cp,1);

% mZ0=[1 , 0];
% mT  =[1 1; 0 1];     mHH = diag([dsigma_eta^2, dsigma_zeta^2]);
% va = zeros(2,1);   

%mZ0=1;
% mT  =dphi;     mHH =  dsigma_eta^2 ;
% va = 0;  
% mP = mHH; 

[mZ0, mT, mHH, va, mP,mH1] = fSsfARMA(dphi, dtheta,  dsigma_eta^2);
%mT = []; mHH = []; mP = []; va = [];

for jd=1:ck

    f    = @(vP1)fARMA_Approx(dd(jd), cn, cp, vP1); % function
    [vP1, ~, ~, ~] = fminunc(f, vP10, opts );
    vtAR    = fFisherInvTransform(vP1(1:cp));
    vtMA    = fFisherInvTransform(vP1(cp+1:2*cp));
    %[vphi, vtheta]
    [~, mTs, mHHs, vas, ~] = fSsfARMA(fReparAR(vtAR), fReparMA(vtMA), dsigma(jd)^2);
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
