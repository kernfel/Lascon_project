function [vsn] = simulate_STN(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  v2=-62+randn(n,1)*5; % line 103

  p=get_params(n);
  CAsn2=0.005*ones(n,1);

  vsn=zeros(n,length(t)); % line 154

  vsn(:,1)=v2; % line 168
  
  % Lines 187ff
  N2=stn_ninf(vsn(:,1));
  H2=stn_hinf(vsn(:,1));
  M2=stn_minf(vsn(:,1));
  A2=stn_ainf(vsn(:,1));
  B2=stn_binf(vsn(:,1));
  C2=stn_cinf(vsn(:,1));
  D2=stn_d2inf(vsn(:,1));
  D1=stn_d1inf(vsn(:,1));
  P2=stn_pinf(vsn(:,1));
  Q2=stn_qinf(vsn(:,1));
  R2=stn_rinf(vsn(:,1));
  
  for i=2:length(t)
    V2=vsn(:,i-1);
    
    % Lines 501-523
    n2=stn_ninf(V2);
    m2=stn_minf(V2);
    h2=stn_hinf(V2);
    a2=stn_ainf(V2);
    b2=stn_binf(V2);
    c2=stn_cinf(V2);
    d2=stn_d2inf(V2);
    d1=stn_d1inf(V2);
    p2=stn_pinf(V2);
    q2=stn_qinf(V2);
    r2=stn_rinf(V2);
 
    td2=130;
    tr2=2;
    tn2=stn_taun(V2);
    tm2=stn_taum(V2);
    th2=stn_tauh(V2);
    ta2=stn_taua(V2);
    tb2=stn_taub(V2);
    tc2=stn_tauc(V2);
    td1=stn_taud1(V2);
    tp2=stn_taup(V2);
    tq2=stn_tauq(V2);

	Ecasn=p.con*log(p.Cao./CAsn2);
    
    % Lines 538ff
    Ina2=p.gna(2)*(M2.^3).*H2.*(V2-p.Ena(2));
    Ik2=p.gk(2)*(N2.^4).*(V2-p.Ek(2));
    Ia2=p.ga*(A2.^2).*(B2).*(V2-p.Ek(2));
    IL2=p.gL*(C2.^2).*(D1).*(D2).*(V2-Ecasn);
    It2=(p.gt(2)*(P2.^2).*(Q2).*(V2-Ecasn));
    Icak2=p.gcak*(R2.^2).*(V2-p.Ek(2));
    Il2=p.gl(2)*(V2-p.El(2));
    Igesn=0; %(ggesn*((V2-Esyn(1)).*(S3a+S31a))); 
    Icorsnampa=0; %gcorsna.*(V2-Esyn(2)).*(S6b+S61b);
    Icorsnnmda=0; %gcorsnn.*(V2-Esyn(2)).*(S6bn+S61bn);
    
    % Lines 623ff
    vsn(:,i)=V2+dt*(1/p.Cm*(-Ina2-Ik2-Ia2-IL2-It2-Icak2-Il2-Igesn-Icorsnampa-Icorsnnmda+Idbs(i)+Iinj(i))); %STN-DBS
    N2=N2+dt*((n2-N2)./tn2); 
    H2=H2+dt*((h2-H2)./th2);
    M2=M2+dt*((m2-M2)./tm2); 
    A2=A2+dt*((a2-A2)./ta2);
    B2=B2+dt*((b2-B2)./tb2); 
    C2=C2+dt*((c2-C2)./tc2);
    D2=D2+dt*((d2-D2)./td2); 
    D1=D1+dt*((d1-D1)./td1);
    P2=P2+dt*((p2-P2)./tp2); 
    Q2=Q2+dt*((q2-Q2)./tq2);
    R2=R2+dt*((r2-R2)./tr2); 
    
    CAsn2=CAsn2+dt*((-p.alp*(IL2+It2))-(p.Kca*CAsn2));
  end
end

function [ainf] = stn_ainf(V)
    ainf=1./(1+exp(-(V+45)./14.7));
end

function [binf] = stn_binf(V)
    binf=1./(1+exp((V+90)./7.5));
end

function [cinf] = stn_cinf(V)
    cinf=1./(1+exp(-(V+30.6)./5));
end

function [d1inf] = stn_d1inf(V)
    d1inf=1./(1+exp((V+60)./7.5));
end

function [d2inf] = stn_d2inf(V)
    d2inf=1./(1+exp((V-0.1)./0.02));
end

function [hinf] = stn_hinf(V)
    hinf=1./(1+exp((V+45.5)./6.4));
end

function [minf] = stn_minf(V)
    minf=1./(1+exp(-(V+40)./8));
end

function [ninf] = stn_ninf(V)
    ninf=1./(1+exp(-(V+41)./14));
end

function [pinf] = stn_pinf(V)
    pinf=1./(1+exp(-(V+56)./6.7));
end

function [qinf] = stn_qinf(V)
    qinf=1./(1+exp((V+85)./5.8));
end

function [rinf] = stn_rinf(V)
    rinf=1./(1+exp(-(V-0.17)./0.08));
end

function [tau] = stn_taua(V)
    tau=1+1./(1+exp(-(V+40)./-0.5));
end

function [tau] = stn_taub(V)
    tau=200./(exp(-(V+60)./-30)+exp(-(V+40)./10));
end

function [tau] = stn_tauc(V)
    tau=45+10./(exp(-(V+27)./-20)+exp(-(V+50)./15));
end

function [tau] = stn_taud1(V)
    tau=400+500./(exp(-(V+40)./-15)+exp(-(V+20)./20));
end

function [tau] = stn_tauh(V)
    tau=24.5./(exp(-(V+50)./-15)+exp(-(V+50)./16));
end

function [tau] = stn_taum(V)
    tau=0.2+3./(1+exp(-(V+53)./-0.7));
end

function [tau] = stn_taun(V)
    tau=11./(exp(-(V+40)./-40)+exp(-(V+40)./50));
end

function [tau] = stn_taup(V)
    tau=5+0.33./(exp(-(V+27)./-10)+exp(-(V+102)./15));
end

function [tau] = stn_tauq(V)
    tau=400./(exp(-(V+50)./-15)+exp(-(V+50)./16));
end