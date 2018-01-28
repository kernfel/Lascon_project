function [vi, vpyrfsi] = simulate_fsi(dt, tmax, Iinj)
  n = 1;

  t=0:dt:tmax;
  
  if ( !isequal(size(Iinj), size(t)) )
    i = Iinj(1);
    Iinj = ones(size(t)) .* i;
  end

  Idbs = zeros(size(t));

  % 109ff
  ai=0.1;
  bi=0.2;
  ci=-65;
  di=2;

  p=get_params(n);

  vi=zeros(n,length(t)); % line 159
  ui=zeros(n,length(t));
  vpyrfsi=zeros(n,length(t));

  vi(:,1)=ci; % 174f
  ui(:,1)=bi*vi(1);
  
  for i=2:length(t)
    V8=vi(:,i-1);
        
    Iei=0;
    
    vpyrfsi(:,i) = dt/p.Cm * 0.1 * 0.43 * V8;
    
    % Lines 793ff
    vi(:,i)=V8+dt*((0.04*(V8.^2))+(5*V8)+140-ui(:,i-1)-Iei+Iinj(i));
    ui(:,i)=ui(:,i-1)+dt*(ai*((bi*V8)-ui(:,i-1)));
    
   for j=1:n
        if vi(j,i-1)>=30
        vi(j,i)=ci;
        ui(j,i)=ui(j,i-1)+di;
        end
   end
  end
end
