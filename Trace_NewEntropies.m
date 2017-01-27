LogE = 0; LogV = 0;
ArcE = 0;
Gam = 1;

ifig = 1;


% Logistique - ordre 1
if(LogE == 1)
  l0 = 0; l1 = -1; c = 0;
  %
  u = 0:1e-3:1;
  sq = sqrt(1-u); at = atanh(sq);
  %
  phi_lm = c + l0*u + l1*(u.*at - sq).*(u<=1);
  %
  figure(ifig)
  h = plot(u,phi_lm,'k-'); set(h,'linewidth',1);
  dt = (max(phi_lm) - min(phi_lm))*.01;
  set(gca,'xlim',[0 1],'ylim',[min(phi_lm) max(phi_lm)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu('tracé de la logistique concave - moyennes :','oui','non');
  if(ct==1); print Logistic_moy -loose -deps; end
  %
end

% -----

% Logistique - ordre 2
if(LogV == 1)
  l0 = 0; l1 = -1; c = 0;
  u = 0:1e-4:1;
  sq = sqrt(1-u); at = atanh(sq);
  %
  phi_lv = c + l0*u + l1*( u.*(at.^2) - 2*sq.*at - log(u) ).*(u<=1);
  %
  figure(ifig)
  h = plot(u,phi_lv,'k-'); set(h,'linewidth',1);
  dt = (max(phi_lv) - min(phi_lv))*.01;
  set(gca,'xlim',[0 1],'ylim',[min(phi_lv) max(phi_lv)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu('tracé de la logistique - variance :','oui','non');
  if(ct==1); print Logistic_var -loose -deps; end
end


% Arcsin centre - ordre 1
if(ArcE == 1)
  c = 0; l0 = 0; l1 = 1;
  u = 1:1e-4:2;
  %
  su = sqrt(u.^2-1);
  phi = c + l0*u + (su+atan(1./su)).*(u>=1);
  %
  figure(ifig)
  h = plot(u,phi,'k-'); set(h,'linewidth',1);
  dt = (max(phi) - min(phi))*.01;
  set(gca,'xlim',[u(1) u(end)],'ylim',[min(phi)-dt max(phi)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu(['tracé de l''arcsine centré dilaté - moyenne :'],'oui','non');
  if(ct==1); print Arcsine_moy -loose -deps; end
end

% -----

% Gamma
if(Gam == 1)
  al = [2 5]; pp = 2; k = [1 1.125];% coef multiplicateur pb échelle
  c1 = 0; l0 = 0; l01 = 1; l11 = -.1;
  %
  Np = 500; x = (1:Np)/Np;
  %phi0 = zeros(length(u),length(al),length(pp));
  %phi1 = phi0;
  %
  for ip = 1:length(pp)
    p = pp(ip);
    %
    for ia = 2:2;%1:length(al)
      a = al(ia);
      %
      u = x*exp(1-a);
      c0 = c1- p*gamma(p+a-1)*l11/((a-1)^(p+a-1));
      %
      W0 = -real(lambertw(0,-u.^(1/(a-1))));
      W1 = -real(lambertw(-1,-u.^(1/(a-1))));
      %
      phi0 = c0 + l0*u + l01*u.*(W0.^p).*(1-p*hypergeom(1,p+a,(a-1)*W0)/(p+a-1));
      phi1 = c1 + l0*u + l11*u.*(W1.^p).*(1-p*hypergeom(1,p+a,(a-1)*W1)/(p+a-1));
      %
      figure(ifig)
      h = plot(u,phi0,'k-',u,phi1,'k-'); set(h,'linewidth',1);
      dt = (max([phi0 phi1]) - min([phi0 phi1]))*.01;
      set(gca,'xlim',[0 u(end)],'ylim',[min([phi0 phi1])-dt max([phi0 phi1])+dt],'fontsize',6);
      %yt = get(gca,'ytick')'; set(gca,'yticklabel',num2str(yt));
      set(gcf,'paperposition',[0 0 5 3*k(ia)]);
      ifig = ifig+1;
      ct = menu(['tracé de la gamma(' int2str(a) ') - moment p=' int2str(p) ' : '],'oui','non');
      if(ct==1); eval(['print Gamma_a' int2str(a) '_p' int2str(p) ' -loose -deps;']); end
    end
  end
end