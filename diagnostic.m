k=length(vP_hat);%+length(vBeta)); k_det=length(mPar);
AIC = 2*k- 2*dLogLik_exact
BIC = k*log(length(vy))- 2*dLogLik_exact

%%%%%%%%%%%%%%%%%for model 5 FULL DETERMINISTIC
disp(['sum2corr',string(SumCorr2(jjj))])
disp(['likelihood',string(dLogLik_exact)])
disp(['The uncertency of the trend from KF is +-',string(sqrt(mVarBeta(2,2)))])
disp(['Trend',string(vBeta(2))])
UNCERTAINTY_paper=  (sqrt(dsigma2_eta0_hat)*sqrt((1+dphi_hat)/(1-dphi_hat)))/ sqrt(sum((juliandate(vtime)-mean(juliandate(vtime))).^2));
disp(['uncertainty of the trend from AR(1)',string( UNCERTAINTY_paper)])
    %close all))


%%%%%%%%%%%%%%%%for other models STOCHASTIC PERIODICAL COMPONENTS
disp(['likelihood',string(dLogLik_exact*1e-4)])
disp(['sum2corr',string(SumCorr2(jjj))])
disp(['Trend',string(vBeta(2)*1e3)])
UNCERTAINTY_paper=  (sqrt(dsigma2_eta0_hat)*sqrt((1+dphi_trend_hat)/(1-dphi_trend_hat)))/ sqrt(sum((juliandate(vtime)-mean(juliandate(vtime))).^2));
disp(['Uncertainty of the trend from AR(1)',string( UNCERTAINTY_paper*1e3)])
disp(['The uncertency of the trend from KF is +-',string(sqrt(mVarBeta(2,2))*1e5)])
 %disp(['Trend',string(vBeta2(2))])

