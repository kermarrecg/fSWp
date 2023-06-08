function [vP0,vLL]= fgrid_InitConds(vy, dLambda , mX,cp  )
mm= length( dLambda  )  ;
 

vd0=linspace(.14,.64,4);
vd1=linspace(.14,.64,4);
vr=linspace(1.2,2,3);
vr2=linspace(0.01,0.99,4);
 
if mm==1; [  d1, vr1,r2] = ndgrid( vd1,vr,vr2 );  mD= [ d1(:),vr1(:),r2(:) ]';
elseif mm==2; [  d1,d2 , vr1,r2] = ndgrid( vd0, vd1 ,vr,vr2);  mD= [ d1(:),  d2(:),vr1(:) ,r2(:)]';
elseif mm==3; [  d1,d2,d3 , vr1,r2] = ndgrid( vd0, vd1, vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),vr1(:) ,r2(:)]';
elseif mm==4; [  d1,d2,d3,d4 , vr1,r2] = ndgrid( vd0, vd1, vd1, vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),  d4(:),vr1(:) ,r2(:)]'; 
elseif mm==5; [  d1,d2,d3,d4,d5 , vr1,r2] = ndgrid( vd0, vd1, vd1, vd1,vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:),vr1(:),r2(:)]';
elseif mm==6; [  d1,d2,d3,d4,d5,d6, vr1,r2 ] = ndgrid( vd0,vd1, vd1, vd1, vd1,vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),vr1(:),r2(:)]';
elseif mm==7; [  d1,d2,d3,d4,d5,d6,d7, vr1,r2 ] = ndgrid( vd0,vd1,vd1, vd1, vd1, vd1,vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:),vr1(:),r2(:)]'; 
elseif mm==8; [  d1,d2,d3,d4,d5,d6,d7,d8, vr1,r2 ] = ndgrid( vd0,vd1,vd1,vd1, vd1, vd1, vd1,vd1,vr,vr2);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:),d8(:),vr1(:),r2(:)]'; 
end 
NN=size(mD,2);
vLL=NaN(NN,1);
f_LLik2= @(theta)   fSsfLogLik_Harvey(theta, vy, dLambda, cp, mX)  ;  %%Log likelihood
mS=NaN(mm,NN);
%vthetaphi0=NaN(NN,1);
for i=1:NN

    %       vLL(i)=-f_LLik2(  [log((mD(1,i)*exp(-0.3*(0:(mm-1))'))./(1-mD(1,i)*exp(-0.3*(0:(mm-1))')));...
    %             log((1-mD(2,i)) *nanvar(vy)*exp(-(0:(mm-1))'))*0.5; ...
    %          log(mD(2,i)* nanvar(vy))*0.5  ]   );
    for j=1:mm
        if mD(j,i)<0.5; mS(j,i) = log( ((gamma(1-mD(j,i))^2)/(gamma(1-2*mD(j,i))))* mD(end,i) *nanvar(vy))*0.5; else
            mS(j,i) = log(  mD(end,i) *0.1*nanvar(vy))*0.5; end
    end
    %vthetaphi0(i)=mD(mm+1,i)/sqrt(1+mD(mm+1,i)^2);
%     vLL(i)=-f_LLik2(  [log(0.55./(1-0.55)); log(mD(1:mm,i)./(1-mD(1:mm,i) ));...%log(mD(2,i))./(1-mD(2,i) );...
%      %   log(mD(3,i))./(1-mD(3,i) );  log(mD(4,i))./(1-mD(4,i) ) ;...
%          log( 0.5*nanvar(vy))*0.5; mS(:,i)    ]   );
 vLL(i)=-f_LLik2(  [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i);  log(  mD(end,i) *0.1*nanvar(vy))*0.5 ; mD(end-1,i)]   );
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(end));
i=vi_star(1); 

vP0 =  [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i);  log(  mD(end,i) *0.1*nanvar(vy))*0.5  ; mD(end-1,i)   ]    ;
end