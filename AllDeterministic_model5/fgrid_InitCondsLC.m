function [vP0,vLL]= fgrid_InitCondsLC(vy , mXreg, iTrend, mZSeas,cp,jjj )
mm= size( mZSeas,2  )/2  ;

 
vd0=linspace(.14,.64,4);
vd1=linspace(.14,.64,4);
vr  = linspace(-10,5,4);
vr2 = linspace(-10,5,4);
coef=0.1; 
 
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
f_LLik2= @(theta)  fSsfLogLik(theta, vy, mXreg, iTrend, mZSeas,cp);    %%Log likelihood
mS=NaN(mm,NN);
 
for i=1:NN
     for j=1:mm
        if mD(j,i)<0.5; mS(j,i) = log( ((gamma(1-mD(j,i))^2)/(gamma(1-2*mD(j,i))))* nanvar(vy))*0.5; else
            mS(j,i) = log( nanvar(vy))*0.5; end
    end
 
 vLL(i)=-f_LLik2(  [  mD(end,i)   ;    mD(end-1,i)   ; log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i);  log(coef*nanvar(vy))*0.5 ; 1.3]   );
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(end));
i=vi_star(1); 

vP0 = [  mD(end,i)   ;    mD(end-1,i)   ; log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i); log(coef*nanvar(vy))*0.5 ; 1.3 ] ;
end