%q = 2;

u = [.001 .02:.02:((q-1)/exp(1))^(q-1)];
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

a1 = -1;
c1 = p * a1 * gamma(p+q-1); % annuler en u = 0;
b = - a1; % annuler au bord u = ( (q-1) / e )^(q-1)
%
% Vraie phi
phi = c1 + b*u + a1 * u.*(((1-q)*W).^p).*(1-p*hypergeom(1,p+q,(1-q)*W)/(p+q-1));

% Cas p = 2-q
phiq = a1 * ( u.*(((1-q)*W).^(1-q)).*((1-q)*W+2-q) +q-2 ) + b * u + c1;

% Approximation avec -log(u) cas p qcq
phiun = c1 + b * u + a1 * u.*((-log(u)).^p).*(1-hypergeom(1,p+1,-log(u)));

% Approximation avec -log(u) cas p=2-q
phil = c1 - a1 + (b + a1) * u - a1 * u.*log(u);

figure(1); plot(u,phi,'r-',u,phiun,'b--',u,phil,'k-.');

figure(2); plot(u,phiq,'r-',u,phil,'b--');

%subplot(3,1,1); plot(u,W,'r-',u,S,'b--');
%subplot(3,1,2); plot(u,(1-S./W)*100,'r-');
%subplot(3,1,3);
%figure(2); plot(u,phi,'r-')
%plot(u,phiT,'r-',u,phi,'b--')
%figure(2); plot(u,(1-q)*W,'r-',u,-log(u),'b--')