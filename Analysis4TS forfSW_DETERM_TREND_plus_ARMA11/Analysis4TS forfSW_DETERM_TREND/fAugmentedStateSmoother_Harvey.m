function [mAs, mPs, vBeta, mVarBeta,aP ] = fAugmentedStateSmoother_Harvey(vy, mZt, ....
    mGG, mT, mHH, va, mP, mW0, mXreg)
%----------------------------------------------------------------------
cT = length(vy);  % number of time points 
cm = size(mT,1);
ckd = size(mW0,2);   % ckd is number of diffuse state elements
ckreg = size(mXreg,2);      % explanatory variables in meas eqn
ck = ckd + ckreg;  % number of regression effects
%----------------------------------------------------------------------
% Initialisation of matrices and scalars
av  = NaN([1 1 cT]);     aV = NaN([1 ck cT]);
af  = NaN([1 1 cT]);     aK = NaN([cm 1 cT]);
aaf = NaN([cm 1 cT]);    aAf = NaN([cm ck cT]); 
aP  = NaN([cm cm cT]);   vs = 0; mS = 0;    
cn = 0;  % observations counter
% -------------------------------------------------------------------------
if ( ckreg == 0); mA = -mW0;
    else ; mA = [-mW0 zeros(cm, ckreg)];
end
mW = 0 ;
for i = 1:cT
    if ( ckreg == 0); mX = zeros(1, ckd);
    else;  mX = [zeros(1, ckd) mXreg(i,:)];
    end
    mZ = mZt(i,:);
	aaf(:,:,i) = va; 	aAf(:,:,i) = mA; 	 aP(:,:,i) = mP;	  
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv = vy(i) - mZ * va;           vV = mX - mZ * mA;			
       % disp(dv);
		df= mZ * mP * mZ' + mGG;
		%vK =( mT * mP * mZ' +  mHGt)/ df;
        vK =( mT * mP * mZ' )/ df;
        av(:,:,i) = dv; 		aV(:,:,i) = vV;	
		af(:,:,i) = df; 		aK(:,:,i) = vK;
		va = mT * va + vK * dv; 		mA = mW + mT * mA + vK * vV;
		mP = mT * mP * mT' + mHH - vK * vK' * df ;
	 	vs = vs + dv * vV' / df ;
		mS = mS + vV' * vV/df;
 		cn = cn + 1;  % observation counter
%          
	else % y(i) is   missing
		va = mT * va;
		mA = mW + mT * mA;
		mP = mT * mP * mT' + mHH; 
    end 
end

[mQ,mR]  = qr(mS); opts.UT  = true;
mRinv    = linsolve(mR,eye(size(mR)),opts);
mVarBeta = mRinv*mQ';  	 vBeta    = mVarBeta * vs;
% vRes = reshape(av,cT,1) - reshape(aV,ck, cT)'* vBeta ;
% vVarRes = reshape(af,cT,1);
%aStateSmo = repmat(NaN, [cm 1 cT]);
%aCovStateSmo = repmat(NaN, [cm cm cT]);
	
vr = zeros(cm,1);  mR = zeros(cm, ck); mN = zeros(cm,cm); 
mAs = NaN(cm, cT); 	mPs = NaN( cm, cT); 
for i = cT:-1.0:1
	%disp(['Smoothing observation: ', num2str(i)]);
    mZ = mZt(i,:);
	if (isnan(vy(i)) == 0)
		
	 	dfinv =  1 / af(:,:,i);
		mL = mT - aK(:,:,i) * mZ;			
		vr = mZ' * dfinv * av(:,:,i) + mL' * vr;	
		mR = mZ' * dfinv * aV(:,:,i) + mL' * mR;
	 	mN = mZ' * dfinv * mZ + mL' * mN * mL;
    else 
        vr = mT' * vr;
		mR = mT' * mR; 
		mN = mT' * mN * mT; 
    end
	mAs(:,i) = aaf(:,:,i) - aAf(:,:,i)*vBeta + aP(:,:,i) * (vr - mR * vBeta );
	mAstar = aAf(:,:,i) + aP(:,:,i) * mR; 
    mPs( :,i) =  diag(aP(:,:,i) - aP(:,:,i) * mN * aP(:,:,i) + mAstar * mVarBeta * mAstar' );
end
 
end