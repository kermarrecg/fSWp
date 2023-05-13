function [vInnovations, vVarInnovations, mStatePred, mCovStatePred,...
          dSigma2, dLogLikC, vBeta, mVarBeta] = fAugmentedKalmanFilterTpS(vy, mZt, mGG, mT, mHH, va, mP, mW0, mXreg)
%----------------------------------------------------------------------
cT = length(vy);  % number of time points 
cm = size(mT,1);
ckd = size(mW0,2);   % ckd is number of diffuse state elements
ckreg = size(mXreg,2);      % explanatory variables in meas eqn
ck = ckd + ckreg;  % number of regression effects
%--------------------------------------------------------------------------
% Initialisation of matrices and scalars
dLogf = 0; dSumSquares = 0; vs = 0; mS = 0;    
cn = 0;  % observations counter
vInnovations = NaN(1, cT);  % stores the KF innovations
vVarInnovations = NaN(1, cT); 
mStatePred = NaN(cm, cT); % stores the states predictions
mCovStatePred = NaN(cm, cT);
% -------------------------------------------------------------------------
if ( ckreg == 0) mA = -mW0;
    else  mA = [-mW0 zeros(cm, ckreg)];
end
mW = 0; 
for i = 1:cT
    %disp(['Filtering observation: ', num2str(i)]);
    if ( ckreg == 0) mX = zeros(1, ckd);
    else  mX = [zeros(1, ckd) mXreg(i,:)];
    end
    mZ =mZt(i,:);
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv = vy(i) - mZ * va;           vV = mX - mZ * mA;			
        %disp(dv);
		df= mZ * mP * mZ' + mGG;
		vK = mT * mP * mZ' / df;
		va = mT * va + vK * dv; 		mA = mW + mT * mA + vK * vV;
		mP = mT * mP * mT' + mHH - vK * vK' * df ;
		dLogf = dLogf + log(df);        
		dSumSquares = dSumSquares + dv^2 /df;
	 	vs = vs + dv * vV' / df ;
		mS = mS + vV' * vV/df;
		cn = cn + 1;  % observation counter
	else % y(i) is   missing
		va = mT * va;
		mA =  mW + mT * mA;
		mP = mT * mP * mT' + mHH; 
    end 
    if (cn == ck) 
%        mVarBeta = inv(mS);  	vBeta = mS\vs;
        [mQ,mR] = qr(mS);
        opts.UT = true;
        mRinv = linsolve(mR,eye(size(mR)),opts);
        mVarBeta = mRinv*mQ';  	vBeta = mVarBeta * vs;
        vaf = va - mA * vBeta; 		mPf = mP + mA * mVarBeta * mA';
    end
    if cn>ck
        vInnovations(i) = dv - vV  * vBeta;  vVarInnovations(i) = df + vV  * mVarBeta * vV';
		mStatePred(:,i) = vaf;		mCovStatePred(:,i) = diag(mPf);
        [mQ,mR] = qr(mS);
        opts.UT = true;
        mRinv = linsolve(mR,eye(size(mR)),opts);
        mVarBeta = mRinv*mQ';  	vBeta = mVarBeta * vs;
%         mVarBeta = inv(mS);  	vBeta = mVarBeta * vs;
        vaf = va - mA * vBeta; 		mPf = mP + mA * mVarBeta * mA';
    end
end
dSigma2 = (dSumSquares-vs'*vBeta)/(cn-ck);
dLogLikC = 0.5 * ( (cn - ck) * ( log(2*pi) +1 + log(dSigma2)) + dLogf + log(det(mS)));  
dPev =  dSigma2 * (df + vV  * mVarBeta * vV');
end
 