function [vstr_indr] = simulate_StrD2(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  v5=-63.8+randn(n,1)*5; % line 106

  p=get_params(n);

  vstr_indr=zeros(n,length(t)); % line 157

  vstr_indr(:,1)=v5; % line 171
  
  % Lines 199ff
  m5=alpham(vstr_indr(:,1))./(alpham(vstr_indr(:,1))+betam(vstr_indr(:,1)));
  h5=alphah(vstr_indr(:,1))./(alphah(vstr_indr(:,1))+betah(vstr_indr(:,1)));
  n5=alphan(vstr_indr(:,1))./(alphan(vstr_indr(:,1))+betan(vstr_indr(:,1)));
  p5=alphap(vstr_indr(:,1))./(alphap(vstr_indr(:,1))+betap(vstr_indr(:,1)));
  
  for i=2:length(t)
    V5=vstr_indr(:,i-1);
    
    % Lines 575ff
    Ina5=p.gna(4)*(m5.^3).*h5.*(V5-p.Ena(4));
    Ik5=p.gk(4)*(n5.^4).*(V5-p.Ek(4));
    Il5=p.gl(4)*(V5-p.El(4));
    %Im5=(2.6-1.1*pd)*gm*p5.*(V5-Em);
    Im5=2.6*p.gm*p5.*(V5-p.Em); % No PD related downregulation of gm
    Igaba5=0; %(ggaba/4)*(V5-Esyn(7)).*(S11cr+S12cr+S13cr+S14cr);
    Icorstr5=0; %gcorindrstr*(V5-Esyn(2)).*(S6a);
    
    % Lines 709ff
   vstr_indr(:,i)=V5+(dt/p.Cm)*(-Ina5-Ik5-Il5-Im5-Igaba5-Icorstr5+Iinj(i));
   m5=m5+dt*(alpham(V5).*(1-m5)-betam(V5).*m5);
   h5=h5+dt*(alphah(V5).*(1-h5)-betah(V5).*h5);
   n5=n5+dt*(alphan(V5).*(1-n5)-betan(V5).*n5);
   p5=p5+dt*(alphap(V5).*(1-p5)-betap(V5).*p5);
   %S1c=S1c+dt*((Ggaba(V5).*(1-S1c))-(S1c/tau_i));
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