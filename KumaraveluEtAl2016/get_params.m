function [p] = get_params(n)
	% Returns a struct with parameters set as in lines 123ff
	p.Cm=1; %Membrane capacitance


	%Ionic conductance and Equilibrium potential values
	p.gl=[0.05 0.35 0.1 0.1]; p.El=[-70 -60 -65 -67];
	p.gna=[3 49 120 100]; p.Ena=[50 60 55 50]; 
	p.gk=[5 57 30 80]; p.Ek=[-75 -90 -80 -100];
	p.gt=[5 5 0.5]; p.Et=0;
	p.gca=[0 2 0.15]; p.Eca=[0 140 120];
	p.Em=-100; p.gm=1;					              				% gm from line 374
	p.gahp=[0 20 10]; p.k1=[0 15 10]; p.kca=[0 22.5 15];
	p.ga=5;
	p.gL=15;
	p.gcak=1;

	p.Kca=2*10^-3;
	p.Z=2;
	p.F=96485;
	p.Cao=2000;
	p.R=8314;
	p.T=298;

	p.alp=1/(p.Z*p.F);
	p.con=(p.R*p.T)/(p.Z*p.F);
