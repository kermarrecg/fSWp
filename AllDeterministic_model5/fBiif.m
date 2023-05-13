function vpsi = fBiif(vt, dc)
%vpsi = (abs(vt) <= dc).* vt .* (1-(vt/dc).^2).^2; 
vpsi=vt;
for i=1:length(vt)
if abs(vt(i))>dc; vpsi(i)=dc*sign(vt(i));end
end
end

