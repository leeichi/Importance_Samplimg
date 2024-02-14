function phiprime=phi_prime(the,lambda,b,x,a0)
temp=the*b.^2.*(1-the*lambda)./((1-2*the*lambda).^2)+lambda./(1-2*the*lambda);
phiprime=sum(temp)-(x-a0);