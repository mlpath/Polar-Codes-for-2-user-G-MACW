function z = integrnd2(y,x1,x2) 
z = triplequad(@(y,x1,x2)((exp(-(x1.^2+x2.^2)/2).*exp(-((y - x1 - x2).^2)/2)/(2*pi)^1.5).*log2( (exp(-((y-x1-x2).^2)/2)/(sqrt(2*pi)))./integrnd(y))),-10,10,-5,5,-5,5); 
end
