function vals = FFTgaussian(x,mu,sig)
%This takes the FFT of a Gaussian, which is, interestingly, a Gaussian.

vals=zeros(1,length(x));

for i=1:length(x)
 vals(i)=exp(-.5*(x(i)^2)*(sig^2));
end
