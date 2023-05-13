function dmse = fARMA_Approx(dd, cn, cp, vP)
% ARMA(p,p) approximation of FN(d)
vPsi    = fFN_MAcoeff(dd, cn);  % FN MA coefficients
vtAR    = fFisherInvTransform(vP(1:cp));
vphi    = fReparAR(vtAR);
vtMA    = fFisherInvTransform(vP(cp+1:2*cp));
vtheta  = fReparMA(vtMA);
vPsi_t  = filter([1; vtheta], [1; -vphi],  [1; zeros(cn,1)])'; 
dmse    = dot(cn:-1:1, (vPsi(2:end)-vPsi_t(2:end)).^2)/cn;
end