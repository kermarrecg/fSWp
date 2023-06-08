function [vphi] = fReparAR(vr)   
% Monahan's reparameterisation enforcing stationarity (1984 Bmtk)
% Maps r_j in (-1,1) to AR $\phi_j$ within the stationary region
cp = length(vr); % order of AR
vphi = zeros(cp, 1);
vphi(1) = vr(1);
for k = 2:cp   
    vphi(1:k-1) = vphi(1:k-1) - vr(k) * flip(vphi(1:k-1));
	vphi(k) = vr(k);
end 
end 