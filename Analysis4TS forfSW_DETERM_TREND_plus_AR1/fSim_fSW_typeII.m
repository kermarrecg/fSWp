function [vy]=fSim_fSW_typeII(dd,dLambda,dsigma2,cn)

K=length(dLambda);
vFN=NaN(1,5000+cn);
vFNs=NaN(1,5000+cn);
e= randn(5000+cn,K);
es= randn(5000+cn,K);
mY=NaN(K,5000+cn);
for j=1:K
   mM = sqrt(dsigma2(j))*[cos(dLambda(j)*(-4999:cn)')  sin(dLambda(j)*(-4999:cn)')];
   vMAcoeff = fFN_MAcoeff(dd(j), 5000+cn-1) ; 
   for t=1:5000+cn
      vFN(t)= vMAcoeff(1:t)*flipud(e(1:t,j));
      vFNs(t)= vMAcoeff(1:t)*flipud(es(1:t,j));
   end
   mY(j,:) = mM(:,1).*vFN' +  mM(:,2).*vFNs';
end
%% simulate y
mY2=mY(:,5001:end);
if K>1
vy = sum(mY2)' ; else; vy=mY2'; end


end