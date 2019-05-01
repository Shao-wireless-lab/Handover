function f = gaussmf(x,sigma,c)

f = exp(-(x-c).^2/(2.*sigma.^2));

end