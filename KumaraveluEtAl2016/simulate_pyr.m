function [ve, vfsipyr, vthcor] = simulate_pyr(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  % 109ff
  ae=0.02;
  be=0.2;
  ce=-65;
  de=8;

  p=get_params(n);

  ve=zeros(n,length(t)); % line 159
  ue=zeros(n,length(t));
  vfsipyr=zeros(n,length(t));
  vthcor=zeros(n,length(t));

  ve(:,1)=ce; % line 172f
  ue(:,1)=be*ve(1);
  
  for i=2:length(t)
    V7=ve(:,i-1);
        
    % Lines 590ff
    Iie=0;%gie*(V7-Esyn(1)).*(S11br+S12br+S13br+S14br);
    Ithcor=0;%gthcor*(V7-Esyn(2)).*(S7);
    
    vfsipyr(:,i) = dt/p.Cm * 0.2 * 0.43 * (V7+85);
    vthcor(:,i) = dt/p.Cm * .15 * .43 * V7;
    
    % Lines 759ff
    ve(:,i)=V7+dt*((0.04*(V7.^2))+(5*V7)+140-ue(:,i-1)-Iie-Ithcor+Iinj(i));
    ue(:,i)=ue(:,i-1)+dt*(ae*((be*V7)-ue(:,i-1)));
    
   for j=1:n
        if ve(j,i-1)>=30
        ve(j,i)=ce;
        ue(j,i)=ue(j,i-1)+de;
        end
   end
  end
end
