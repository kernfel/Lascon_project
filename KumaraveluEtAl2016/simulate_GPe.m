function [vge] = simulate_GPe(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  v3=-62+randn(n,1)*5; % line 104

  p=get_params(n);

  vge=zeros(n,length(t)); % line 155

  vge(:,1)=v3; % line 169
  
  % Lines 179ff
  N3=gpe_ninf(vge(:,1));
  H3=gpe_hinf(vge(:,1));
  R3=gpe_rinf(vge(:,1));
  CA2=0.1; 
  CA3=CA2;
  
  for i=2:length(t)
    V3=vge(:,i-1);
    
    % Lines 486ff
    m3=gpe_minf(V3);
    n3=gpe_ninf(V3);
    h3=gpe_hinf(V3);
    a3=gpe_ainf(V3);
    s3=gpe_sinf(V3);
    r3=gpe_rinf(V3);
    
    tn3=gpe_taun(V3);
    th3=gpe_tauh(V3);
    tr3=30;
    
    % Lines 550ff
    Il3=p.gl(3)*(V3-p.El(3));
    Ik3=p.gk(3)*(N3.^4).*(V3-p.Ek(3));
    Ina3=p.gna(3)*(m3.^3).*H3.*(V3-p.Ena(3));
    It3=p.gt(3)*(a3.^3).*R3.*(V3-p.Eca(3));
    Ica3=p.gca(3)*(s3.^2).*(V3-p.Eca(3));
    Iahp3=p.gahp(3)*(V3-p.Ek(3)).*(CA3./(CA3+p.k1(3)));
    Isngeampa=0; %(gsngea).*((V3-Esyn(2)).*(S2a+S21a)); 
    Isngenmda=0; %(gsngen).*((V3-Esyn(2)).*(S2an+S21an)); 
    Igege=0; %(0.25*(pd*3+1))*(ggege).*((V3-Esyn(3)).*(S31c+S32c)); 
    Istrgpe=0; %gstrgpe*(V3-Esyn(6)).*(S5+S51+S52+S53+S54+S55+S56+S57+S58+S59);
    Iappgpe=0; %3-2*corstim*~pd; %Modulation only during cortical stim to maintain mean firing rate
    
    % Lines 658ff
    vge(:,i)=V3+dt*(1/p.Cm*(-Il3-Ik3-Ina3-It3-Ica3-Iahp3-Isngeampa-Isngenmda-Igege-Istrgpe+Iappgpe+Iinj(i)));
    N3=N3+dt*(0.1*(n3-N3)./tn3);
    H3=H3+dt*(0.05*(h3-H3)./th3);
    R3=R3+dt*(1*(r3-R3)./tr3);
    CA3=CA3+dt*(1*10^-4*(-Ica3-It3-p.kca(3)*CA3));
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

