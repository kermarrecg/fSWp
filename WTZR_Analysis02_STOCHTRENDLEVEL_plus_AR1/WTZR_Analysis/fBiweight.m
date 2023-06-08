function vpsi = fBiweight(vt, dc)
vpsi = (abs(vt) <= dc).* vt .* (1-(vt/dc).^2).^2; 
end

