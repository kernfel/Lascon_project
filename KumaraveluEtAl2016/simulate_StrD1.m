function [vstr_dr] = simulate_StrD1(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  v6=-63.8+randn(n,1)*5; % line 107

  p=get_params(n);

  vstr_dr=zeros(n,length(t)); % line 158

  vstr_dr(:,1)=v6; % line 172
  
  % Lines 203ff
  m6=alpham(vstr_dr(:,1))./(alpham(vstr_dr(:,1))+betam(vstr_dr(:,1)));
  h6=alphah(vstr_dr(:,1))./(alphah(vstr_dr(:,1))+betah(vstr_dr(:,1)));
  n6=alphan(vstr_dr(:,1))./(alphan(vstr_dr(:,1))+betan(vstr_dr(:,1)));
  p6=alphap(vstr_dr(:,1))./(alphap(vstr_dr(:,1))+betap(vstr_dr(:,1))); 
  
  for i=2:length(t)
    V6=vstr_dr(:,i-1);
    
    % Lines 575ff
    Ina6=p.gna(4)*(m6.^3).*h6.*(V6-p.Ena(4));
    Ik6=p.gk(4)*(n6.^4).*(V6-p.Ek(4));
    Il6=p.gl(4)*(V6-p.El(4));
    %Im6=(2.6-1.1*pd)*gm*p6.*(V6-pEm);
    Im6=2.6*p.gm*p6.*(V6-p.Em);
    Igaba6=0; %(ggaba/3)*(V6-Esyn(7)).*(S81r+S82r+S83r);
    Icorstr6=0; %gcordrstr.*(V6-Esyn(2)).*(S6a);
    
    % Lines 734ff
   vstr_dr(:,i)=V6+(dt/p.Cm)*(-Ina6-Ik6-Il6-Im6-Igaba6-Icorstr6+Iinj(1));
   m6=m6+dt*(alpham(V6).*(1-m6)-betam(V6).*m6);
   h6=h6+dt*(alphah(V6).*(1-h6)-betah(V6).*h6);
   n6=n6+dt*(alphan(V6).*(1-n6)-betan(V6).*n6);
   p6=p6+dt*(alphap(V6).*(1-p6)-betap(V6).*p6);
   %S8=S8+dt*((Ggaba(V6).*(1-S8))-(S8/tau_i));
  end
end

function [ah] = alphah(V)
ah=0.128*exp((-50-V)/18);
end

function [am] = alpham(V)
am=(0.32*(54+V))./(1-exp((-54-V)/4));
end

function [an] = alphan(V)
an=(0.032*(52+V))./(1-exp((-52-V)./5));
end

function [ap] = alphap(V)
ap=(3.209*10^-4*(30+V))./(1-exp((-30-V)./9));
end

function [bh] = betah(V)
bh=4./(1+exp((-27-V)/5));
end

function [bn] = betan(V)
bn=0.5*exp((-57-V)./40);
end

function [bm] = betam(V)
bm=0.28*(V+27)./((exp((27+V)/5))-1);
end

function [bp] = betap(V)
bp=(-3.209*10^-4*(30+V))./(1-exp((30+V)./9));
end