%q = 2;

u = 0 : .01 : ((q-1)/exp(1))^(q-1);
%0:.001:1/exp(1);

x = -(u.^(1/(q-1)))/(q-1);

W = lambertw(-1,x);

S = log(-x)-log(-log(-x));

K = 0; M = 0;

for k = 0:K
   for m = 1:M
      ckm = ((-1)^(k+m-1))*StirlingFirst(k+m,k+1)/gamma(m+1);
      S = S + ckm*((log(-log(x))).^m)./((log(x)).^(m+k));
   end
end

phi = u.*((-W).^(1-q)).*(W-(2-q)/(q-1))+u/(q-1);
%phi = -u.*((-W).^(2-q)).*(1-(2-q)*hypergeom(1,2,(1-q)*W))-(2-q)/(q-1)^(2-q);

subplot(3,1,1); plot(u,W,'r-',u,S,'b--');
subplot(3,1,2); plot(u,(1-S./W)*100,'r-');
subplot(3,1,3); plot(u,(q-1)*phi,'r-',u,u.*log(u),'b-');
%figure(2); plot(u,phi,'r-')