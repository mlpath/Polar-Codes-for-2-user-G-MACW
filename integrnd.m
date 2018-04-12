function z = integrnd(y)
z = dblquad(@(x1,x2)(exp(-((y - x1 - x2).^2)/(2)).*exp(-(x1.^2 + x2.^2)/2)/(2*pi*sqrt(2*pi))),-5,5,-5,5);
end