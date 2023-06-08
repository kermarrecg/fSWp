function [ dLogLik, valpha, mPalpha, vRes, vResVar] = fKalmanFilterReg( vy, mXreg,mZt, mGG, mT, mHH, va, mP )
%----------------------------------------------------------------------
cT = length(vy);  % number of time points 


%mHGt=mH*mG';
mW0=[];
cm  = size(mT,1);
ckd = size(mW0,2);   % ckd is number of diffuse state elements
ckreg = size(mXreg,2);      % explanatory variables in meas eqn
ck = ckd + ckreg;  % number of regression effects
%----------------------------------------------------------------------
% Initialisation of matrices and scalars
valpha=NaN(length(va),cT);
mPalpha=NaN(length(va),cT);
vRes=NaN(cT,1);
vResVar=NaN(cT,1);
dLogf = 0; dSumSquares = 0; vs = 0; mS = 0;    
cn = 0;  % observations counter
if ( ckreg == 0) mA = -mW0;
    else  mA = [-mW0 zeros(cm, ckreg)];
end
mW = 0;
for i = 1:cT
    %disp([' Filtering observation: ', num2str(i)]);
    if ( ckreg == 0) mX = zeros(1, ckd);
    else  mX = [zeros(1, ckd) mXreg(i,:)];
    end
    mZ = mZt(i,:);
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv = vy(i) - mZ * va;           vV = mX - mZ * mA;			
		df= mZ * mP * mZ' + mGG;
 		%vK = (mT * mP * mZ' + mHGt) / df;
        vK = (mT * mP * mZ') / df;
		va = mT * va + vK * dv; 		mA = mW + mT * mA + vK * vV;
		mP = mT * mP * mT' + mHH - vK * vK' * df ;
		dLogf = dLogf + log(df);        dSumSquares = dSumSquares + dv^2 /df;
	 	vs = vs + dv * vV' / df ;		mS = mS + vV' * vV/df;
		cn = cn + 1;  % observation counter
        valpha(:,i)=va;
        mPalpha(:,i)=diag(mP);
        vRes(i,1)=dv;
        vResVar(i,1)=df;
	else % y(i) is   missing
		va = mT * va;
		mA =  mW + mT * mA;
		mP = mT * mP * mT' + mHH; 
        valpha(:,i)=va;
        mPalpha(:,i)=diag(mP);
    end 
end
%
% [mQ,mR] = qr(mS);
% opts.UT = true;
% mRinv = linsolve(mR,eye(size(mR)),opts);
% mVarBeta = mRinv*mQ';  	vBeta = mVarBeta * vs;

mVarBeta = inv(mS);  	vBeta = mVarBeta * vs;
dLogLik = 0.5 * ((cn - ck) * log(2*pi)+dLogf+log(det(mS))+dSumSquares-vs'*vBeta);  