LogE = 2; LogV = 2;
GamE = 2; GamV = 1;
ArcE = 2; ArcV = 2;
AcdE = 2; AcdV = 1;

ifig = 1;

% Logistique - ordre 1
if(LogE == 1)
  s = 1;
  y = 0:1e-3:1/4/s;
  sq = sqrt(1-4*s*y); at = atanh(sq);
  %
  phi_lm = sq/2 - 2*s*y.*at;
  %phi_lp = sq/2 + 2*s*y.*at;
  %
  figure(ifig)
  h = plot(y,phi_lm,'k-'); set(h,'linewidth',1);
  dt = (max(phi_lm) - min(phi_lm))*.01;
  set(gca,'xlim',[0 1/4/s],'ylim',[min(phi_lm)-dt max(phi_lm)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu('tracé de la logistique concave - moyennes :','oui','non');
  if(ct==1); print Logistic_moy -loose -deps; end
  %
  %figure(ifig)
  %h = plot(y,phi_lm,'k-',y,phi_lp,'k--'); set(h,'linewidth',1);
  %dt = (max([phi_lm phi_lp]) - min([phi_lm phi_lp]))*.01;
  %set(gca,'xlim',[0 1/4/s],'ylim',[min([phi_lm phi_lp])-dt max([phi_lm phi_lp]+dt)],'fontsize',6);
  %set(gcf,'paperposition',[0 0 5 3]);
  %ifig = ifig+1;
  %ct = menu('tracé de la logistique globale - moyenne :','oui','non');
  %if(ct==1); print Logistic_globale_moy -loose -deps; end
end

% -----

% Logistique - ordre 2
if(LogV == 1)
  s = 1;
  y = 0:1e-3:1/4/s;
  sq = sqrt(1-4*s*y); at = atanh(sq);
  %
  phi_lv = -4*s*y.*(at.^2)+2*sq.*at+log(4*s*y);
  %
  figure(ifig)
  h = plot(y,phi_lv,'k-'); set(h,'linewidth',1);
  dt = (max(phi_lv) - min(phi_lv))*.01;
  set(gca,'xlim',[0 1/4/s],'ylim',[min(phi_lv)-dt max(phi_lv)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu('tracé de la logistique - variance :','oui','non');
  if(ct==1); print Logistic_var -loose -deps; end
end

% -----

% Gamma alpha = 2 et 3 - ordre 2
if(GamE == 1)
  al = [2 5];
  beta = 3;
  for ia = 1:length(al)
    alpha = al(ia);
    fb = beta / gamma(alpha) * (((alpha-1)/exp(1))^(alpha-1));
    y = [0:1e-3:fb fb];
    %
    wy = -((y*gamma(alpha)/beta).^(1/(alpha-1))) / (alpha-1);
    W0 = real(lambertw(0,wy)); W1 = real(lambertw(-1,wy));
    %
    phi0 = zeros(size(wy)); phi1 = phi0;
    %
    c0 = 0;
    for im = 0:alpha-1
      mu = (-1)^(alpha-im)*gamma(alpha)/gamma(im+1)/(alpha-1)^(alpha-im);
      phi0 = phi0 + mu * W0.^(im+1-alpha);
      phi1 = phi1 + mu * W1.^(im+1-alpha);
      c0 = c0 + (-1)^(im+1-alpha) * mu;
    end
    c0 = (c0-1) * beta / gamma(alpha) * ((alpha-1)/exp(1))^(alpha-1);
    phi0 = -(phi0 + W0).*y + 2*c0; phi1 = (phi1 + W1).*y;
    %
    figure(ifig)
    h = plot(y,phi0,'k-',y,phi1,'k--'); set(h,'linewidth',1);
    dt = (max([phi0 phi1]) - min([phi0 phi1]))*.01;
    set(gca,'xlim',[0 y(end)],'ylim',[min([phi0 phi1])-dt max([phi0 phi1])+dt],'fontsize',6);
    set(gcf,'paperposition',[0 0 5 3]);
    ifig = ifig+1;
    ct = menu(['tracé de la gamma(' int2str(alpha) ') - moyenne :'],'oui','non');
    if(ct==1); eval(['print Gamma_' int2str(alpha) '_moy -loose -deps;']); end
    end
end

% -----

% Gamma alpha = 2 et 3 - ordre 2
if(GamV == 1)
  al = [5 10];
  beta = 3;
  for ia = 1:length(al)
    alpha = al(ia);
    fb = beta / gamma(alpha) * (((alpha-1)/exp(1))^(alpha-1));
    y = [0:1e-3:fb fb];
    %
    wy = -((y*gamma(alpha)/beta).^(1/(alpha-1))) / (alpha-1);
    W0 = real(lambertw(0,wy)); W1 = real(lambertw(-1,wy));
    %
    phi0 = zeros(size(wy)); phi1 = phi0;
    %
    c0 = 0;
    for im = 0:alpha
      mu = 2*(-1)^(alpha-im+1)*gamma(alpha+1)/gamma(im+1)/(alpha-1)^(alpha-im+1);
      phi0 = phi0 + mu * W0.^(im+1-alpha);
      phi1 = phi1 + mu * W1.^(im+1-alpha);
      c0 = c0 + (-1)^(im+1-alpha) * mu;
    end
    c0 = (c0-1) * beta / gamma(alpha) * ((alpha-1)/exp(1))^(alpha-1);
    c0 = 0;%%%%%%%%%
    phi0 = (phi0 + W0.^2).*y + 2*c0; phi1 = -(phi1 + W1.^2).*y;
    %
    figure(ifig)
    h = plot(y,phi0,'k-',y,phi1,'k--'); set(h,'linewidth',1);
    dt = (max([phi0 phi1]) - min([phi0 phi1]))*.01;
    set(gca,'xlim',[0 y(end)],'ylim',[min([phi0 phi1])-dt max([phi0 phi1])+dt],'fontsize',6);
    set(gcf,'paperposition',[0 0 5 3]);
    ifig = ifig+1;
    ct = menu(['tracé de la gamma(' int2str(alpha) ') - variance :'],'oui','non');
    if(ct==1); eval(['print Gamma_' int2str(alpha) '_var -loose -deps;']); end
    end
end

% -----

% Arcsin - ordre 1
if(ArcE == 1)
  y = 2/pi:1e-4:2;
  %
  sy = sqrt((pi*y).^2-4);
  phi = (sy+2*atan(2./sy))/pi-1;
  %
  figure(ifig)
  h = plot(y,phi-y+2/pi,'k-',y,phi+y-2/pi,'k--'); set(h,'linewidth',1);
  set(gca,'xlim',[y(1) y(end)],'ylim',[min([phi-y+2/pi]) max([phi+y-2/pi])],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu(['tracé de l''arcsine - moyenne :'],'oui','non');
  if(ct==1); print Arcsine_moy -loose -deps; end
end

% -----

% Arcsin - ordre 2
if(ArcV == 1)
  y = 2/pi:1e-4:1;
  %
  sy = sqrt((pi*y).^2-4);
  phi = (sy+2*atan(2./sy))/pi-1;
  dec = y+2/pi^2./y;
  %
  figure(ifig)
  h = plot(y,phi-dec+3/pi,'k-',y,phi+dec-3/pi,'k--'); set(h,'linewidth',1);
  set(gca,'xlim',[y(1) y(end)],'ylim',[min([phi-dec+3/pi]) max([phi+dec-3/pi])],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu(['tracé de l''arcsine - variance :'],'oui','non');
  if(ct==1); print Arcsine_var -loose -deps; end
end


% -----

% Arcsin centre dilate - ordre 1
if(AcdE == 1)
  s = 1;
  y = 1/(s*pi*sqrt(2)):1e-4:.5;
  %
  sy = sqrt(2*(pi*s*y).^2-1);
  phi = (sy+atan(1./sy))/pi;
  %
  figure(ifig)
  h = plot(y,phi,'k-'); set(h,'linewidth',1);
  dt = (max(phi) - min(phi))*.01;
  set(gca,'xlim',[y(1) y(end)],'ylim',[min(phi)-dt max(phi)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu(['tracé de l''arcsine centré dilaté - moyenne :'],'oui','non');
  if(ct==1); print ArcsineCD_moy -loose -deps; end
end

% -----

% Arcsin centre dilate - ordre 2
if(AcdV == 1)
  s = 1;
  y = 1/(s*pi*sqrt(2)):1e-4:.5;
  %
  phi = 2*s^2*y+1/pi^2./y;
  %
  figure(ifig)
  h = plot(y,phi,'k-'); set(h,'linewidth',1);
  dt = (max(phi) - min(phi))*.01;
  set(gca,'xlim',[y(1) y(end)],'ylim',[min(phi)-dt max(phi)+dt],'fontsize',6);
  set(gcf,'paperposition',[0 0 5 3]);
  ifig = ifig+1;
  ct = menu(['tracé de l''arcsine centré dilaté - variance :'],'oui','non');
  if(ct==1); print ArcsineCD_var -loose -deps; end
end