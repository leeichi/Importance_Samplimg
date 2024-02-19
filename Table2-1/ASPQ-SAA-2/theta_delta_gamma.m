%%Revised 2021/09/23
% sigma and S could be vector
function [theta,delta,gamma]=theta_delta_gamma(r,q,sigma,S,K,tau)
d1=(log(S/K)+(r-q+sigma.^2/2)*tau)./(sigma*sqrt(tau));
d2=d1-sigma*sqrt(tau);
%theta
theta.c=(-exp(-q*tau)*S.*normpdf(d1).*sigma/(2*sqrt(tau))-r*K*exp(-r*tau)*normcdf(d2)+q*S.*exp(-q*tau).*normcdf(d1));
theta.p=(-exp(-q*tau)*S.*normpdf(-d1).*sigma/(2*sqrt(tau))+r*K*exp(-r*tau)*normcdf(-d2)-q*S.*exp(-q*tau).*normcdf(-d1));
%delta
delta.c=exp(-q*tau)*normcdf(d1);
delta.p=-exp(-q*tau)*normcdf(-d1);
%gamma
gamma.c=exp(-q*tau)*normpdf(d1)./(S.*sigma*sqrt(tau));
gamma.p=gamma.c;