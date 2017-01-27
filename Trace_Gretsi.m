beta = 1; alpha = 5.5;%[2 5];
p = 2;
l0 = 0; lk1 = [1 ; -.1]; c = [0 ; 0];

Np = 500;
un = ones(1,Np);
phi = zeros(2,Np,length(alpha));
phimin = zeros(2,length(alpha));

for ia = 1:length(alpha)
  a = alpha(ia);
  %
  T = beta * ( (a-1)^(a-1) ) / gamma(a) / exp(a-1);
  tau = ( ( gamma(a) / beta )^( 1 / (a-1) ) ) / (a-1);
  %
  y = [ 0 : T / (Np-1) : T ];
  yl = - tau * ( y.^( 1 / (a-1) ) ); % argument de la Lambert
  lky = lk1 * y;% facteur l_{1,k} y sur deux ligne (k=0 et k=-1)
  W = -[ real(lambertw(0,yl)) ; real(lambertw(-1,yl))]; % -Lamberts
  %
  % pour avoir phi_k nul à l'origine
  c = [0 ; lk1(2)*p*beta*gamma(p+a-1)/gamma(a)/((a-1)^p)];
  %
  % branches de la fonctionnelle entropique
  phi(:,:,ia) = c*un + l0*[1;1]*y + (lk1*y).*(W.^p).*(1-p*hypergeom(1,p+a,(a-1)*W)/(p+a-1));
  %
  I0 = find(y==0); % indices y = 0
  %
  % équivalents à l'origine:
  phi(:,I0,ia) = 0;
  %
  % position des minima de phi
  r = (-l0 ./ lk1).^(1/p);
  phimin(:,ia) = ((r .* exp(-r) / tau).^(a-1)) .* (r>0);
  %
  % Tracés
  % ------
  figure(ia)
  h = plot(y,phi(1,:,ia),'k-',y,phi(2,:,ia),'k-');
  set(h,'linewidth',1);
  set(gca,'xlim',[y(1) y(end)]);
  set(gca,'ylim',[min(min(phi(:,:,ia))) max(max(phi(:,:,ia)))]);
  set(gca,'fontsize',12);
  set(gcf,'paperposition',[0 0 8 4.8]);
  tt = menu(['Impression de la figure a = ' num2str(a) ' ?'],'oui','non');
  if (tt==1)
    eval(['print Phi_Gamma_a' int2str(a) '_p' int2str(p) ' -deps -loose']);
  end
end