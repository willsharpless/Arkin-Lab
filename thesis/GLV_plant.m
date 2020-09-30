function [xdot] = GLV_plant(t,x,mu,A)
%GLV model for ophelias reduced community
%willsharpless@berkeley.edu

%population dead
x(x<0)=0;

xdot = x.*(mu + A*x);

end