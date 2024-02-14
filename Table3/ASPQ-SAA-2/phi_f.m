function phi=phi_f(theta,lambda,b)
temp=0.5*((theta*b).^2./(1-2*theta*lambda) - log(1-2*theta*lambda));
phi=sum(temp);