function   vACF = fFN_ACF(dd, cq)
vACF = NaN(1,cq+1);
cm = min(50,cq);
vj = 1:cm;
vrho = ( gamma(1 - dd) .* gamma(vj + dd) ) ./ (gamma(vj + 1-dd) .* gamma(dd));
vACF(1:(cm+1))  = [1 vrho];
if (cq>cm)
    drho = vrho(end);
    for j = (cm+1):cq;
        drho = drho * (j-1+dd)/(j-dd);
        vACF(j+1) = drho;
    end
end
end


