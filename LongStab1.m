function [results,stability] = LongStab1(results,JID,lattice,state,geo,ref)



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%						%%%
%%%		Derivatives		%%%
%%%						%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = 5e-02;

i=0;

%%	1. forward speed
i=i+1;
state.STW = state.STW + delta;

[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;

state.STW = state.STW - delta;

%%	2. vertical speed !
i=i+1;
state.W = state.W + delta;

[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;

state.W=state.W - delta;
%%	3. pitch vel
i=i+1;
state.Q = state.Q+delta;

[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;

state.Q = state.Q-delta;

%%	4. pitch	!
i=i+1;
state.theta = state.theta+delta;

[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;

state.theta = state.theta-delta;

%%	5. elavation	!
i=i+1;
state.ELA=state.ELA+delta;

[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;

state.ELA=state.ELA-delta;
%%	6. base case
i=i+1;


[lattice,ref,geo]=fLattice_setup2(geo,state,1);
results= solverloop6(results,1,JID,lattice,state,geo,ref);
stability.force(:,i)=results.FORCE;
stability.moment(:,i)=results.MOMENTS;


for j=1:i-1
	stability.forceder(:,j)=(stability.force(:,j)-stability.force(:,i))./delta;
	stability.momentder(:,j)=(stability.moment(:,j)-stability.moment(:,i))./delta;
end


%%% matrix derivative outlook
%		u	v	w	Q	P	phi		theta	h
%	X K
%	Y M
%	Z N
%%%



stability.elem = categorical({'\Delta u'; '\Delta z'; '\Delta \theta';'\Delta w'; '\Delta q'});
%	ALEC - this thing is not waterproof

%	possibly add a choice option on what stabilty calculations that need to
%	be performed: Bague, Masuyama, Drela
stability.StabAlec = [-1/geo.mass*stability.forceder(1,1) -1/geo.mass*stability.forceder(1,5) ...
	-1/geo.mass*stability.forceder(1,4) -1/geo.mass*stability.forceder(1,2) ...
	-1/geo.mass*stability.forceder(1,3); 0 0 0 -1 0; 0 0 0 0 1;...
	-1/geo.mass*stability.forceder(3,1) ...
	-1/geo.mass*stability.forceder(3,5) -1/geo.mass*stability.forceder(3,4) ...
	-1/geo.mass*stability.forceder(3,2) -1/geo.mass*stability.forceder(3,3); ...
	1/geo.I(2,2)*stability.momentder(2,1) 1/geo.I(2,2)*stability.momentder(2,5) ...
	1/geo.I(2,2)*stability.momentder(2,4) 1/geo.I(2,2)*stability.momentder(2,2) ...
	1/geo.I(2,2)*stability.momentder(2,3)];

%[V,W] = eig(stability.StabAlec);
save('last_calc');
% stab = [-0.27 -2.21 -4.99 -0.5 0;
% 	0 0 0 1 0;
% 	0 0 0 0 1;
% 	-2.24 -57.2 -156.21 -15.62 -6.51;
% 	0.01 -6.42 -78.59 -7.86 -14.75];
% 
% EigenVisual(stability);


end



