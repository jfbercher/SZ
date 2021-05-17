LogE = 0; LogV = 0;
ArcE = 0;
Gam = 1;

ifig = 1;


% Logistique - ordre 1
if(LogE == 1)
  a = 1; b = 0; c = 0; %l0 = 0; l1 = -1; c = 0;
  %
  u = 0:1e-3:1;
  sq = sqrt(1-u); at = atanh(sq);
  %
  phi_lm = -a*(u.*at - sq).*(u<=1) + b*u + c;
  %
  figure(ifig)
  h = plot([u 2],[phi_lm 0],'k-'); set(h,'linewidth',1.5);
  dt = max(phi_lm)*.01;
  set(gca,'xlim',[0 2],'ylim',[0 max(phi_lm)+dt],'fontsize',14);
  set(gcf,'paperposition',[0 0 5 3]); box on;
  ifig = ifig+1;
  ct = menu('trace de la logistique concave - moyennes :','oui','non');
  if(ct==1); print Logistic_moy -loose -deps; end
  %
end

% -----

% Logistique - ordre 2
if(LogV == 1)
  a = 1; b = 0; c = 0; %l0 = 0; l1 = -1; c = 0;
  u = 0:1e-4:1;
  sq = sqrt(1-u); at = atanh(sq);
  %
  phi_lv = -a*( u.*(at.^2) - 2*sq.*at - log(u) ).*(u<=1) + b*u + c;
  %
  figure(ifig)
  h = plot([u 2],[phi_lv 0],'k-'); set(h,'linewidth',1.5);
  dt = max(phi_lv)*.01;
  set(gca,'xlim',[0 2],'ylim',[0 max(phi_lv)+dt],'fontsize',14);
  set(gcf,'paperposition',[0 0 5 3]);box on;
  ifig = ifig+1;
  ct = menu('trace de la logistique - variance :','oui','non');
  if(ct==1); print Logistic_var -loose -deps; end
end


% Arcsin centre - ordre 1
if(ArcE == 1)
  a = 1; b = 0; c = -pi/2;%c = 0; l0 = 0; l1 = 1;
  u = 1:1e-4:2.25;
  %
  su = sqrt(u.^2-1);
  phi = a*(su+atan(1./su)).*(u>=1) + b*u + c;
  %
  figure(ifig)
  h = plot([0 u],[0 phi],'k-'); set(h,'linewidth',1.5);
  dt = max(phi)*.01 ;%- min(phi))*.01;
  set(gca,'xlim',[0 u(end)],'ylim',[0 max(phi)+dt],'fontsize',14);
  set(gcf,'paperposition',[0 0 5 3]); box on;
  ifig = ifig+1;
  ct = menu(['trace de l''arcsine centre dilate - moyenne :'],'oui','non');
  if(ct==1); print Arcsine_moy -loose -deps; end
end

% -----

% Gamma
if(Gam == 1)
  qq = [1.02 1.25 1.5 1.75 2 2.25 2.5]; pp = [1 2]; %k = [1 1.125];% coef multiplicateur pb échelle
  e = exp(1);
  %
  um = max(((qq-1)/e).^(qq-1)); % u maximum
  %
  a0 = 1; a1 = 1; b = a1; c0 = 0;
  %
  Np = 100; x = (1:Np)/Np/e;
  W0 = lambertw(0,-x);
  W1 = lambertw(-1,-x);
  %
  % A blanc pour loquer les figures
  if 0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ip = 1:length(pp)
    figure(2*(ip-1)+1); clf; hold on;
    figure(2*ip); clf; hold on;
  end
  %
  s = {'k-','k-.','k--','k:'};
  m0 = inf(size(pp)); m1=m0; % minimun axe y
  M0 = -inf(size(pp)); M1=M0; % maximum axe y
  for iq = 1:length(qq)
    %m0 = inf; m1 = inf;
    %M0 = -inf; M1 = -inf;% max axe y
    q = qq(iq); disp([' boucle pour q = ' num2str(q)]);
    %
    %u = ((q-1)*x).^(q-1);
    u = (1:Np)/Np*((q-1)/e)^(q-1);
    x = (u.^(1/(q-1)))/(q-1);
    mW0 = real((1-q)*lambertw(0,-x));
    mW1 = real((1-q)*lambertw(-1,-x));
    %'pause',pause
    %
    for ip = 1:length(pp)
      p = pp(ip);
      c1 = c0 - a1*p*gamma(p+q-1);
      %
      %
      if( (p==1) & (iq==1) )
        um = max(1,um);
        y = (1:Np)/Np*um; sha = y.*log(y);
	figure(2*(ip-1)+1)
        h = plot([0 um], [0 0] , 'k-'); set(h,'linewidth',.75);
	figure(2*ip)
        h = plot([0 y], [0 sha] , 'k-'); set(h,'linewidth',.75);
	%
	m1(ip) = min(m1(ip),min(sha)); M1(ip) = 0;
      end
      %
      p0 = u.*(mW0.^p).*(1-p*hypergeom(1,p+q,mW0)/(p+q-1));
      p1 = u.*(mW1.^p).*(1-p*hypergeom(1,p+q,mW1)/(p+q-1));
      %I=find(u==0); p1(I) = -p*gamma(p+q-1);% limite u -> 0
      p0 = [0 p0]; p1 = [-p*gamma(p+q-1) p1];
      %
      phi0 = c0  +  b*[0 u] + a0*p0;
      phi1 = c1  +  b*[0 u] - a1*p1;
      %
      figure(2*(ip-1)+1);
      h = plot([0 u],a0*p0,'k-'); set(h,'linewidth',1.5);
      %h = plot([0 u],a0*p0,char(s(iq))); set(h,'linewidth',1.5);
      %%h = plot([0 u],phi0,char(s(iq))); set(h,'linewidth',1.5);
      figure(2*ip);
      h = plot([0 u],phi1,'k-'); set(h,'linewidth',1.5);
      %h = plot([0 u],phi1,char(s(iq))); set(h,'linewidth',1.5);
      %
      m0(ip) = min(m0(ip),min(a0*p0)); m1(ip) = min(m1(ip),min(phi1));
      M0(ip) = max(M0(ip),max(a0*p0)); M1(ip) = max(M1(ip),max(phi1));
    end
  end
  %
  %
  end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pour scaler et deloquer les figures
  um = max(um,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ip = 1:length(pp)
    figure(2*(ip-1)+1);
    dt0=(M0(ip)-m0(ip))/20;
    set(gca,'xlim',[0 um],'ylim',[0 M0(ip)+dt0],'fontsize',14);
    set(gcf,'paperposition',[0 0 5 3]);box on;
    hold off;
    ct = menu(['trace p0 de la gamma - moment p = ' int2str(pp(ip)) ' : '],'oui','non');
    if(ct==1); eval(['print Gamma_phi0_p' int2str(pp(ip)) ' -loose -deps;']); end
    %
    figure(2*ip);
    dt1=(M1(ip)-m1(ip))/20;
    set(gca,'xlim',[0 um],'ylim',[m1(ip)-dt1 0],'fontsize',14);
    set(gcf,'paperposition',[0 0 5 3]);box on;
    ct = menu(['trace p1 de la gamma - moment p = ' int2str(pp(ip)) ' : '],'oui','non');
    if(ct==1); eval(['print Gamma_phi1_p' int2str(pp(ip)) ' -loose -deps;']); end
    hold off;
  end
  figure(2*length(pp)+1)
  q = 1:.01:4.5; D = ((q-1)/e).^(q-1);% domain for phi
  h = plot(q,D,'k-'); set(h,'linewidth',1.5);
  set(gca,'xlim',[1 q(end)],'ylim',[0 max(D)*1.01],'fontsize',14);
  set(gcf,'paperposition',[0 0 5 3]); box on;
  ct = menu(['trace D de la gamma : '],'oui','non');
  if(ct==1); print Gamma_Dphi -loose -deps; end
end

