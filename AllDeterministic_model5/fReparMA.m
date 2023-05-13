function [ vtheta ] = fReparMA( vr )
% Monahan's reparameterisation enforcing invertibility (1984 Bmtk)
% Maps r_j in (-1,1) to MA $\phi_j$ within the invertibility region
cq = length(vr); % order of MA
vtheta = zeros(cq, 1);
vtheta(1) = vr(1);
for k = 2:cq   
    vtheta(1:k-1) = vtheta(1:k-1)  + vr(k) * flip(vtheta(1:k-1));
	vtheta(k) = vr(k);
end 
end

