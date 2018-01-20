function [vth] = simulate_th(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  v1=-62+randn(n,1)*5; % line 102

  p=get_params(n);

  vth=zeros(n,length(t)); % line 153

  vth(:,1)=v1; % line 167
  
  % Lines 180f
  H1=th_hinf(vth(:,1)); 
  R1=th_rinf(vth(:,1));
  
  for i=2:length(t)
    V1=vth(:,i-1);
    
    % Lines 485-499
    m1=th_minf(V1);
    h1=th_hinf(V1);
    p1=th_pinf(V1);
    r1=th_rinf(V1);

    th1=th_tauh(V1);
    tr1=th_taur(V1);
    
    % Lines 528ff
    Il1=p.gl(1)*(V1-p.El(1));
    Ina1=p.gna(1)*(m1.^3).*H1.*(V1-p.Ena(1));
    Ik1=p.gk(1)*((0.75*(1-H1)).^4).*(V1-p.Ek(1));
    It1=p.gt(1)*(p1.^2).*R1.*(V1-p.Et);
    Igith=0; %igith=ggith*(V1-Esyn(6)).*(S4); % GPi -> Thalamus synapse
    Iappth=1.2;
    
    % Lines 600ff
    vth(:,i)= V1+dt*(1/p.Cm*(-Il1-Ik1-Ina1-It1-Igith+Iappth+Iinj(i)));
    H1=H1+dt*((h1-H1)./th1);
    R1=R1+dt*((r1-R1)./tr1);
  end
end

% Lines 866ff
function [hinf] = th_hinf(V)
    hinf=1./(1+exp((V+41)./4));
end

function [minf] = th_minf(V)
    minf=1./(1+exp(-(V+37)./7));
end

function [pinf] = th_pinf(V)
    pinf=1./(1+exp(-(V+60)./6.2));
end

function [rinf] = th_rinf(V)
    rinf=1./(1+exp((V+84)./4));
end

function [a] = ah(V)
    a=0.128*exp(-(V+46)./18); % part of th_tauh fxn
end

function [b] = bh(V)
    b=4./(1+exp(-(V+23)./5)); % part of th_tauh fxn
end

function [tau] = th_tauh(V)
    tau=1./(ah(V)+bh(V));
end

function [tau] = th_taur(V)
    tau=0.15*(28+exp(-(V+25)./10.5));
end