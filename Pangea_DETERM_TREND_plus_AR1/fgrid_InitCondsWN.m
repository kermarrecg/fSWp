function [vP0,vLL]= fgrid_InitCondsWN(vy, dLambda , mX,cp  )
mm= length( dLambda  )  ;
 

vd0=linspace(.14,.84,5);
vd1=linspace(.14,.84,5);
 
 
if mm==1; [  d1 ] = ndgrid( vd1  );  mD= [ d1(:)  ]';
elseif mm==2; [  d1,d2] = ndgrid( vd0, vd1  );  mD= [ d1(:),  d2(:) ]';
elseif mm==3; [  d1,d2,d3  ] = ndgrid( vd0, vd1, vd1 );  mD= [ d1(:),  d2(:),  d3(:) ]';
elseif mm==4; [  d1,d2,d3,d4  ] = ndgrid( vd0, vd1, vd1, vd1 );  mD= [ d1(:),  d2(:),  d3(:),  d4(:) ]'; 
elseif mm==5; [  d1,d2,d3,d4,d5  ] = ndgrid( vd0, vd1, vd1, vd1,vd1 );  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:) ]';
elseif mm==6; [  d1,d2,d3,d4,d5,d6  ] = ndgrid( vd0,vd1, vd1, vd1, vd1,vd1 );  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:) ]';
elseif mm==7; [  d1,d2,d3,d4,d5,d6,d7  ] = ndgrid( vd0,vd1,vd1, vd1, vd1, vd1,vd1 );  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:) ]'; 
elseif mm==8; [  d1,d2,d3,d4,d5,d6,d7,d8  ] = ndgrid( vd0,vd1,vd1,vd1, vd1, vd1, vd1,vd1 );  mD= [ d1(:),  d2(:),  d3(:),  d4(:), d5(:), d6(:),d7(:),d8(:) ]'; 
end 
NN=size(mD,2);
vLL=NaN(NN,1);
f_LLik2= @(theta)   fSsfLogLik_Harvey(theta, vy, dLambda, cp, mX)  ;  %%Log likelihood
mS=NaN(mm,NN);
%vthetaphi0=NaN(NN,1);
for i=1:NN
 
    for j=1:mm
        if mD(j,i)<0.5; mS(j,i) = log( ((gamma(1-mD(j,i))^2)/(gamma(1-2*mD(j,i))))  *nanvar(vy))*0.5; else
            mS(j,i) = log( 0.1*nanvar(vy))*0.5; end
    end
 
 
 vLL(i)=-f_LLik2(  [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i)    ]   );
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(end));
i=vi_star(1); 

vP0 = [ log(mD(1:mm,i)./(1-mD(1:mm,i) ));... 
       mS(:,i)   ]     ;
end