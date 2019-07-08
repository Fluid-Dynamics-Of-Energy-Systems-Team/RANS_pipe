

fA = 0.001;
fB = 0.001;

imax = 10;

for i = 1:imax
    fact = (i-0.)/(imax-0.);
    ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB));
    ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
end

delta(1) = 0.5*ru(1);
for i = 2:imax
    delta(i)=0.5*(ru(i)-ru(i-1))
end


figure(2); hold off
plot([1:10], ru, 'ro')