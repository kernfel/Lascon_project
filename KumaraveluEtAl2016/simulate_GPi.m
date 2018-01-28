function [vgi, vsngi, vgegi, vstrgi] = simulate_GPi(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  v4=-62+randn(n,1)*5; % line 105

  p=get_params(n);

  vgi=zeros(n,length(t)); % line 156
  vsngi = zeros(n,length(t));
  vgegi = zeros(n,length(t));
  vstrgi = zeros(n,length(t));

  vgi(:,1)=v4; % line 170
  
  % Lines 179ff
  N4=gpe_ninf(vgi(:,1));
  H4=gpe_hinf(vgi(:,1));
  R4=gpe_rinf(vgi(:,1));
  CA2=0.1; 
  CA4=CA2;
  
  for i=2:length(t)
    V4=vgi(:,i-1);
    
    % Lines 486ff
    m4=gpe_minf(V4);
    n4=gpe_ninf(V4);
    h4=gpe_hinf(V4);
    a4=gpe_ainf(V4);
    s4=gpe_sinf(V4);
    r4=gpe_rinf(V4);
    
    tn4=gpe_taun(V4);
    th4=gpe_tauh(V4);
    tr4=30;
    
    % Lines 563ff
    Il4=p.gl(3)*(V4-p.El(3));
    Ik4=p.gk(3)*(N4.^4).*(V4-p.Ek(3));
    Ina4=p.gna(3)*(m4.^3).*H4.*(V4-p.Ena(3));
    It4=p.gt(3)*(a4.^3).*R4.*(V4-p.Eca(3));
    Ica4=p.gca(3)*(s4.^2).*(V4-p.Eca(3));
    Iahp4=p.gahp(3)*(V4-p.Ek(3)).*(CA4./(CA4+p.k1(3)));
    Isngi=0; %(gsngi).*((V4-Esyn(4)).*(S2b+S21b));
    Igigi=0; %ggigi*((V4-Esyn(5)).*(S31b+S32b)); 
    Istrgpi=0; %gstrgpi*(V4-Esyn(6)).*(S9+S91+S92+S93+S94+S95+S96+S97+S98+S99);
    Iappgpi=0; %3;
    
    vsngi(:,i) = dt/p.Cm * .43 * .15 * V4;
    vgegi(:,i) = dt/p.Cm * .5 * .3 * (V4 + 85);
    vstrgi(:,i) = dt/p.Cm * .5 * .3 * (V4 + 85);
    
    % Lines 684ff
    vgi(:,i)=V4+dt*(1/p.Cm*(-Il4-Ik4-Ina4-It4-Ica4-Iahp4-Isngi-Igigi-Istrgpi+Iappgpi+Iinj(i)));
    N4=N4+dt*(0.1*(n4-N4)./tn4);
    H4=H4+dt*(0.05*(h4-H4)./th4);
    R4=R4+dt*(1*(r4-R4)./tr4);
    CA4=CA4+dt*(1*10^-4*(-Ica4-It4-p.kca(3)*CA4));
  end
end

function [ainf] = gpe_ainf(V)
    ainf=1./(1+exp(-(V+57)./2));
end

function [hinf] = gpe_hinf(V)
    hinf=1./(1+exp((V+58)./12));
end

function [minf] = gpe_minf(V)
    minf=1./(1+exp(-(V+37)./10));
end

function [ninf] = gpe_ninf(V)
    ninf=1./(1+exp(-(V+50)./14));
end

function [rinf] = gpe_rinf(V)
    rinf=1./(1+exp((V+70)./2));
end

function [sinf] = gpe_sinf(V)
    sinf=1./(1+exp(-(V+35)./2));
end

function [tau] = gpe_tauh(V)
    tau=0.05+0.27./(1+exp(-(V+40)./-12));
end

function [tau] = gpe_taun(V)
    tau=0.05+0.27./(1+exp(-(V+40)./-12));
end

