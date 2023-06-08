function   vARcoeff = fFN_ARcoeff(dd, cq)
vARcoeff = NaN(1,cq+1);
cm = min(50,cq);
vj = 1:cm;
vpi = gamma(vj - dd) ./ (gamma(vj + 1) .* gamma(-dd));
vARcoeff(1:(cm+1))  = [1 vpi];
if (cq>cm)
    dpi = vpi(end);
    for j = (cm+1):cq;
        dpi = dpi * (j-1-dd)/j;
        vARcoeff(j+1) = dpi;
    end
end
end
 