function [vP0,vLL]= fgrid_InitCondsBNJACOB(vy0, dLambda , mX,cp,v1,v2  )
mm= length( dLambda  )  ;
 

vd0=linspace(.24,.64,3);
vd1=linspace(.04,.34,3);
vr=linspace(-5,1.8,3);
vr2=linspace(0.1, 0.5,3);
 
if mm==1; [  d1, vr1,r2] = ndgrid( vd0,vr,vr2 );  mD= [ d1(:),vr1(:),r2(:) ]';
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

f_LLik2= @(theta) fSsfLogLik_HarveyJACOB(theta, vy0, dLambda, cp, mX,v1,v2);  %%Log likelihood
mS=NaN(mm,NN);
%vthetaphi0=NaN(NN,1);
for i=1:NN

    %       vLL(i)=-f_LLik2(  [log((mD(1,i)*exp(-0.3*(0:(mm-1))'))./(1-mD(1,i)*exp(-0.3*(0:(mm-1))')));...
    %             log((1-mD(2,i)) *nanvar(vy)*exp(-(0:(mm-1))'))*0.5; ...
    %          log(mD(2,i)* nanvar(vy))*0.5  ]   );
    for j=1:mm
        if mD(j,i)<0.5; mS(j,i) = log( ((gamma(1-mD(j,i))^2)/(gamma(1-2*mD(j,i))))*nanvar(vy0))*0.5; else
            mS(j,i) = log( (1- mD(end,i)) *nanvar(vy0))*0.5; end
    end
    %vthetaphi0(i)=mD(mm+1,i)/sqrt(1+mD(mm+1,i)^2);
 
 vLL(i)=-f_LLik2(  [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i);  log(  mD(end,i) *nanvar(vy0) )*0.5 ;   .7;  0.1]   );
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(end));
i=vi_star(1); 


vP0 =  [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i);  log(  mD(end,i) *nanvar(vy0))*0.5  ;  .7; 0.1   ]    ;
end