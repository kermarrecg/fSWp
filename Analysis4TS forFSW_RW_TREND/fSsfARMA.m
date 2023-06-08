function [mZ, mT, mHH, va, mP,mH1] = fSsfARMA(vphi, vtheta, dsigma2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARMA MODEL - Harvey representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp = length(vphi);
cq = length(vtheta);
cm = max(cp,cq+1);     % number of state elements
 
% Measurement equation: % y_{t} = Z_t \alpha_{t}  
mZ = [1, zeros(1, cm-1)];    
% Transition equation:  % \alpha_{t+1} = T_t \alpha_{t} + H_t \epsilon_t
mT = zeros(cm,cm); mT(1:cp,1) = vphi;
mT(:,2:end) = [eye(cm-1); zeros(1, cm-1)];     
mH = [1; zeros(cm-1,1)]; mH(2:2+cq-1) = vtheta; 
mHH = dsigma2*(mH*mH'); % H_t * H_t'  
mH1= sqrt(dsigma2)*mH;
%%%%%%%%%%%%%%% Initial conditions               %%%%%%%%%%%%%%%%%%%%%%%%%
va = zeros(cm,1);     

%[mV,mD]=eig(eye(cm^2)-kron(mT,mT));
mP = mHH;%  reshape( inv(eye(cm^2)-kron(mT,mT))*reshape(mHH, cm^2,1) , cm, cm);    %      
 
%{ 
mP=zeros(cm,cm);
vTheta= [1;vtheta(1:end)];
mP(:,cm)=dsigma2*(vTheta *vtheta(end));
mP(cm,:)=dsigma2*(vTheta'*vtheta(end));

for row=1:cm-1
    for col=1:cm-1
       ma_coef_1 = vTheta(row)*[1; arma2ma(vphi(row:end)./vTheta(row),vtheta(row:end)./vTheta(row),100)]; 
       ma_coef_2 = vTheta(col)*[1; arma2ma(vphi(col:end)./vTheta(col),vtheta(col:end)./vTheta(col),100)]; 
       mP(row,col)= dsigma2*sum(ma_coef_1.*ma_coef_2);
    end
end
%}
end


