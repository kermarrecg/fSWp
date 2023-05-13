function [mZt, mGG, mT, mHH, mW, va, mP, mW0] = fTrendPlusTrigoSeasCycles(dSeta, dSzeta, ....
                                                      mZSeas, mHHSeas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cm  = 2+size(mZSeas,2);
cn  = size(mZSeas,1);
cK  = cm;              % number of diffuse elements
%% system matrices
mZt = [(ones(cn,1) * [1, 0])  mZSeas];
mT  = blkdiag([1 1; 0 1], eye(cm-2)); 
mGG = 1; % G_t * G_t'
mHH =  blkdiag(diag([dSeta; dSzeta]), mHHSeas);   % H_t * H_t'  
mW  = zeros(cm, cK);  
%% Initial conditions 
va  = zeros(cm, 1); 
mP  = mHH;
mW0 = mT ;
end
