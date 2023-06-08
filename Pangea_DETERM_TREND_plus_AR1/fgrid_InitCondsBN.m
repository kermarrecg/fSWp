function [vP0,vLL]= fgrid_InitCondsBN(vy, dLambda , mX,cp  )
mm= length( dLambda  ) ;
m2=mm ;

vd0=linspace(.14,.64,4);
vd1=linspace(.14,.64,4);
vr=linspace(0.01,0.2,4);

if     m2==1; [  d1, vr1] = ndgrid( vd1,vr );  mD= [ d1(:),vr1(:) ]';
elseif m2==2; [  d1,d2 , vr1 ] = ndgrid( vd1, vd0 ,vr);  mD= [ d1(:),  d2(:),vr1(:)  ]';
elseif m2==3; [  d1,d2,d3 , vr1 ] = ndgrid( vd1, vd1,vd0, vr);  mD= [ d1(:),  d2(:),  d3(:),vr1(:)  ]';
elseif m2==4; [  d1,d2,d3,d4 , vr1 ] = ndgrid(vd1, vd1, vd1, vd0, vr);  mD= [ d1(:),  d2(:),  d3(:),  d4(:),vr1(:)  ]'; 
elseif m2==5; [  d1,d2,d3,d4,d5 , vr1 ] = ndgrid(vd1, vd1, vd1,vd1, vd0, vr);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:),vr1(:) ]';
elseif m2==6; [  d1,d2,d3,d4,d5,d6, vr1 ] = ndgrid(vd1, vd1, vd1, vd1,vd1, vd0,vr);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),vr1(:) ]';
elseif m2==7; [  d1,d2,d3,d4,d5,d6,d7, vr1  ] = ndgrid(vd1,vd1, vd1, vd1, vd1,vd1, vd0,vr);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:),vr1(:) ]'; 
elseif m2==8; [  d1,d2,d3,d4,d5,d6,d7,d8, vr1 ] = ndgrid(vd1,vd1,vd1, vd1, vd1, vd1,vd1, vd0,vr);  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:),d8(:),vr1(:) ]'; 
end 
NN=size(mD,2);
vLL=NaN(NN,1);
f_LLik2= @(theta)   fSsfLogLik_Harvey(theta, vy, dLambda, cp, mX)  ;  %%Log likelihood
mS=NaN(mm,NN);
%vthetaphi0=NaN(NN,1);
for i=1:NN


    for j=1:mm 
        if mD(j,i)<0.5; mS(j,i) = log( ((gamma(1-mD(j,i))^2)/(gamma(1-2*mD(j,i))))* (1-mD(end,i))*nanvar(vy))*0.5; else
            mS(j,i) = log(  (1-mD(end,i)) *nanvar(vy))*0.5; end
    end
 
 vLL(i)=-f_LLik2( [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
      mS(1:end,i);   log(  ( mD(end,i)) *nanvar(vy))*0.5 ;1.3]);%;   log( 0.01*nanvar(vy))*0.5 ]); %; log(  0.2*nanvar(vy))*0.5 ]);
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(end));
i=vi_star(1); 

vP0 = [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
      mS(1:end,i);   log(  ( mD(end,i)) *nanvar(vy))*0.5 ;1.3];%;   log( 0.01 *nanvar(vy))*0.5 ]; %; log(  0.2*nanvar(vy))*0.5 ]   ;
end