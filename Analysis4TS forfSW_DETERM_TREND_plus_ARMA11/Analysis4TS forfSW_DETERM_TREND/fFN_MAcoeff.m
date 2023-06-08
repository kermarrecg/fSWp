function   vMAcoeff = fFN_MAcoeff(dd, cq)
vMAcoeff = NaN(1,cq+1);
cm = min(50,cq);
vj = 1:cm;
vpsi = gamma(vj + dd) ./ (gamma(vj + 1) .* gamma(dd));
vMAcoeff(1:(cm+1))  = [1 vpsi];
if (cq>cm)
    dpsi = vpsi(end);
    for j = (cm+1):cq
        dpsi = dpsi * (j-1+dd)/j;
        vMAcoeff(j+1) = dpsi;
    end
end
end