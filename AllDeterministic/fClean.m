function [ vyclean, vWeights] = ...
 fClean(vy, mZt, mGG, mT, mHH, mW, va_star, mP_star, mW0, mXreg, dSigma2, dthrshld, cadd)
%*********************************************************************************
% Extends Martin and Thomson (1982) to general nonstationary state space form 
% Robust Augmented Filter - returns a 'clean' series
% cc is Tukey's biweight threshold 
%----------------------------------------------------------------------
cT = length(vy);  % number of time points 
cm = size(mT,1);
ckd = size(mW0,2);
%[cm, ckd] = size(mW0);   % ckd is number of diffuse state elements
ckreg = size(mXreg,2);      % explanatory variables in meas eqn
ck = ckd + ckreg;  % number of regression effects
ckadd = ck + cadd;

%--------------------------------------------------------------------------
cneff = 0;  % observations counter
i = 0;
vs = 0; mS = 0;  
vWeights = NaN(1, cT); vyclean = NaN(1, cT); 
if ckreg == 0 
    mA = -mW0;
else  mA = [-mW0 zeros(cm, ckreg)]; 
end

%--------------------------------------------------------------------------
while cneff < ckadd
    i = i + 1;
%     disp(['Filtering observation: ', num2str(i)]);
    if  ckreg == 0
        mX = zeros(1, ckd);    
    else
        mX = [zeros(1, ckd) mXreg(i,:)];
    end
    mZ  = mZt(i,:); 
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv_star = vy(i) - mZ * va_star;           vV = mX - mZ * mA;			
		df_star= mZ * mP_star * mZ' + mGG;        
        vK = (mT * mP_star * mZ') / df_star;
        cneff = cneff + 1;  % observation counter
   		va_star = mT * va_star + vK * dv_star; 		mA = -mW + mT * mA + vK * vV;
		mP_star = mT * mP_star * mT' + mHH - vK * vK' * df_star ;
	 	vs = vs + dv_star * vV' / df_star ;
		mS = mS + vV' * vV/df_star;
	else % y(i) is   missing
		va_star = mT * va_star; mA =  -mW + mT * mA; 
        mP_star = mT * mP_star * mT' + mHH; 
    end
end 
% disp('end of first run');
[mQ,mR] = qr(mS);
opts.UT = true;
mRinv = linsolve(mR,eye(size(mR)),opts);
% mVarBeta = inv(mS);  	vBeta = mVarBeta * vs;
mVarBeta = mRinv*mQ'; 	vBeta = mVarBeta * vs;


for i = i+1:cT
%     disp(['Filtering observation: ', num2str(i)]);
    if  ckreg == 0
        mX = zeros(1, cm);     
    else  mX = [zeros(1, ckd) mXreg(i,:)];    
    end
    mZ  = mZt(i,:); 
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv_star = vy(i) - mZ * va_star;           vV = mX - mZ * mA;			
		df_star = mZ * mP_star * mZ' + mGG;        
        vK = ( mT * mP_star * mZ' ) / df_star;
        cneff = cneff + 1;  % observation counter
        dv = dv_star - vV  * vBeta; df = df_star + vV * mVarBeta * vV';
        
        % standardized innovations
        dt = dv/sqrt(dSigma2 * df);  
        
        % transform standardized innovations according to the chosen
        % influence function
        %dpsi =  fBiweight(dt,dthrshld);
        dpsi = fBiif(dt, dthrshld);
        dpsis = dpsi * sqrt(dSigma2);
        
        % fixed effects updated estimates
        vBeta = vBeta + mVarBeta * vV' * dpsis/sqrt(df);
        mVarBeta = mVarBeta - (dpsi / dt) * mVarBeta * (vV'*vV) * mVarBeta/df ;
        
        % updated state estimates
        vart = va_star - mA * vBeta + mP_star * mZ' * dpsis/sqrt(df);  
        mPrt = mP_star - (dpsi / dt)* mP_star * (mZ'*mZ) * mP_star/df_star + mA*mVarBeta*mA';
        
        % clean series
        vyclean(i) = mZ * vart + mX * vBeta + mGG * dpsis/sqrt(df);
        vWeights(i) = (dpsi / dt);
        
        % prediction step 
        va_star = mT * va_star + vK * (dpsi / dt) * dv_star; 		mA = -mW + mT * mA + vK * (dpsi / dt)* vV;
		mP_star = mT * mP_star * mT' + mHH - (dpsi / dt) * (vK * vK') * df_star;

    else % y(i) is   missing
		va_star = mT * va_star; mA =  -mW + mT * mA; 
        mP_star = mT * mP_star * mT' + mHH; 
    end 
end
vWeights(1:ckadd) = 1;
vyclean(1:ckadd) = vy(1:ckadd);
end
 