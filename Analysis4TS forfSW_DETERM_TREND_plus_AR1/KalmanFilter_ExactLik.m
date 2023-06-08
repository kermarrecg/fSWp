function [vInnovations, vVarInnovations, mStatePred, mCovStatePred,...
    dLogLik ] = KalmanFilter_ExactLik(vy,mZt , mGG, mT, mHH, va, mP)
%----------------------------------------------------------------------
cm = length(va); % n. of state elements 
cn = length(vy);
%cN=sum(isnan(vy)==0);
%----------------------------------------------------------------------
 
% Initialisation of matrices and scalars
dLogf = 0; dSumSquares = 0;     
vInnovations = NaN(cn,1);  % stores the KF innovations
vVarInnovations = NaN(cn,1); 
mStatePred = NaN(cm, cn); % stores the states predictions
mCovStatePred = NaN(cm, cn);
% vd=nanmean(vy);
cN=sum(isnan(vy)== 0);
for i = 1:cn
    if isnan(vy(i))== 0
    mZ = mZt(i,:);
	dv = vy(i) -  mZ * va  ; 		dF = (mZ) * mP * mZ' + mGG;
                                    vK = (mT * mP * mZ' )/dF;
	va = mT * va  + vK * dv; 	mP = mT * mP * mT' + mHH - vK *  dF * vK' ;
    
    vInnovations(i) = dv;       vVarInnovations(i) = dF;
    dLogf = dLogf + log(dF);    dSumSquares = dSumSquares + dv^2 /dF; %in this line we compute sum(logF) and sum(v*F^(-1)*v')
    else
    va = mT * va   ; 	mP = mT * mP * mT' + mHH;
    end
    % store quantities
    
    mStatePred(:,i) = va;       mCovStatePred(:,i) = diag(mP);
end 
dLogLik = - 0.5 * (cN * log(2 * pi) + dLogf + dSumSquares ); %here we compute the LogLik
 
end
 