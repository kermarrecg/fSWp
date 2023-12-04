function [vP0,vLL]= fgrid_InitCondsPN(vy , mX,cp  )

vd0=linspace(-.3,.7,11);
 
NN=length(vd0);
vLL=NaN(NN,1);
f_LLik2= @(theta)     fSsfLogLik_Harvey(theta, vy, mX,cp);    %%Log likelihood

%vthetaphi0=NaN(NN,1);
for i=1:NN

% vLL(i)=-f_LLik2( [ log(nanvar(vy))*0.5;  log( vd0(i) ./(1- vd0(i)  )) ...
%       ]   );
  vLL(i)=-f_LLik2( [ log(nanvar(vy))*0.5;   vd0(i)   ...
       ]   );
    disp(['grid search optimal initial conditions ', num2str(i),' to ',num2str(NN)])
end

maxLL=maxk(vLL,11);
vi_star=find(vLL>=maxLL(1));
i=vi_star(1); 

% vP0 =[ log(nanvar(vy))*0.5;  log( vd0(i) ./(1- vd0(i)  ))  ]    ;
 vP0 =[ log(nanvar(vy))*0.5;  vd0(i) ]    ;

end