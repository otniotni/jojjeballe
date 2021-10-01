function f = tfunc(x)
global mu beta f0

f = mu*x-(1-beta*x).^sum(f0);