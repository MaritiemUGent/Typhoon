function [results,lattice]=solver10(results,state,geo,lattice,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 1999, 2007 Tomas Melin
%
% This file is part of Tornado
%
% Tornado is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2, or (at your option) any later version.
%
% Tornado is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public
% License along with Tornado; see the file GNU GENERAL 
% PUBLIC LICENSE.TXT.  If not, write to the Free Software 
% Foundation, 59 Temple Place -Suite 330, Boston, MA
% 02111-1307, USA.
%
% usage: [RESULTS] = solver8(results,state,geo,lattice,ref)
%
% This function computes forces and moments on each panel.
% Inputs are coordinades for old resluts, collocationpoints, 
%   vorticies and Normals, reference area and chord
%
% Example:
%
%   [results]=solver8(results,state,geo,lattice,ref);
%
% Calls:
%           Setboundary    
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado core function
%
% Revision History:
%   Bristol,  2007 06 27:  Addition of new header. TM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Running solver 8')
%f1=waitbar(0,'SOLVER10: starting calculation...');

latticem = mirror1(lattice,state);



[a vor_length void]=size(lattice.VORTEX);%extracting number of sections in 
										   %"horseshoes"

	%if vor_length < 8
    %   terror(1)
    %   return
    %end
    
%flops(0)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up regular vortices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u1 U1]=fastdw(lattice,0);      %void is the complete downwash, w2 is only the normal component

results.dwcond=cond(u1);		%ALEC - why?
%disp('dnwash... ok')
%count=flops         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up mirror vortices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[um1 Um1]=fastdw(latticem,1);
%waitbar(0.33,f1,'flow plus image complete, starting wave-making component...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up wavemaking contributions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state.message=0;
try
	state.wavemethod;
catch
	state.wavemethod = 0;
end


if (state.wavemethod)
	[uw1, Uw1, state] = wave(lattice,geo,state);
else
	[test1 test2 test3] = size(Um1);
	uw1 = zeros(test1,test2); Uw1 = zeros(test1,test2,test3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up right hand side %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhs=(setboundary5(lattice,state,geo))';
%disp('rhs... ok')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving for rhs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%	


%waitbar(0.5,f1,'wave-making component complete, start calculating gamma')

%gamma=(u1+uw1)\rhs';      %flow tangency condition test with no waves
gamma=(u1+um1+uw1)\rhs';      %flow tangency condition
%disp('gauss... ok')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporary test environment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tc=64;

%quiver3(lattice.COLLOC(:,1),lattice.COLLOC(:,2),lattice.COLLOC(:,3),lattice.N(:,1),lattice.N(:,2),lattice.N(:,3),'Color','r');
%hold on
results.uw1=uw1;
results.Uw1=Uw1;
results.u1=u1;
results.U1=U1;
results.um1=um1;
results.Um1=Um1;
results.tot = (u1+um1+uw1);
results.TOT = (U1+Um1+Uw1);

% for tc=25:120
% 	figure(1)
% 	set(gcf, 'Position', get(0, 'Screensize')-[0 0 200 100]);
% 	
% 	subplot(1,2,1)
% 	quiver3(lattice.COLLOC(:,1),lattice.COLLOC(:,2),lattice.COLLOC(:,3),Uw1(:,tc,1),Uw1(:,tc,2),Uw1(:,tc,3),'Color','r');
% 	set(gca,'DataAspectRatio',[1,1,1])
% 	subplot(1,2,2)
% 	colormap(hot);
% 	fill3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)',uw1(:,tc)');
% 	colorbar('vert')
% 	h = colorbar;
% 	set(h, 'ylim', [-0.2 0.2])
% 	axis equal
% 	pause(1)
% 	close all
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporary test environment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if state.pgcorr==1
        %tdisp('Trying PG correction')
        %Prandtl-Glauert compressibility correction
        [state.rho sos p_1]=ISAtmosphere(state.ALT);
        M=state.STW/sos;
        corr=1/(sqrt(1-M^2));  
        gamma=gamma*corr;
        %Yeah, this part is not validated yet... or even published, but it
        %seems to work. Do use it with caution though as the math rigour
        %isnt there yet.   
    end



b1=vor_length/2;

p1(:,:)=lattice.VORTEX(:,b1,:);		%Calculating panel vortex midpoint...	
p2(:,:)=lattice.VORTEX(:,b1+1,:);	%to use as a force locus
lattice.COLLOC(:,:)=(p1+p2)./2;	    % LOCAL control point, vortex midpoint.

q1(:,:)=lattice.VORTEX(:,b1,:);		%Calculating panel vortex midpoint...	
q2(:,:)=lattice.VORTEX(:,b1+1,:);	%to use as a force locus
lattice.COLLOC(:,:)=(q1+q2)./2;	    % LOCAL control point, vortex midpoint.



%%%
%	test
%%%

results.sumu1 = sum(sum(u1));
results.sumum = sum(sum(um1));
results.sumuw = sum(sum(uw1));
results.sumgamma = sum(gamma);
results.meanu1 = mean(mean(abs(u1)));
results.meanum = mean(mean(abs(um1)));
results.meanuw = mean(mean(abs(uw1)));
results.meangamma = mean(abs(gamma));


%%%
%	test
%%%




c3=lattice.COLLOC-ones(size(lattice.COLLOC,1),1)*geo.ref_point;



%waitbar(0.66,f1,'calculating flow for forces...')
[u2 U2]=fastdw(lattice,0);	        %Calculating downwash on vorticies, new AIC matrix
[um2 Um2]=fastdw(lattice,1);		%fastdw happens again as it now gets
%um2=-um2; Um2=-Um2;

if (state.wavemethod==1)
	[uw2, Uw2] = wave(lattice,geo,state);
else
	[test1 test2 test3] = size(Um1);
	uw2 = zeros(test1,test2); Uw2 = zeros(test1,test2,test3);
end
% %[uw2 Uw2]=wave(lattice,geo,state);	%calculated for new control point.
% uw2 = zeros(test1,test2); Uw2 = zeros(test1,test2,test3);
% %test no waves

[s1 s2]=size(p1);


%%

DWX=U2(:,:,1);
DWY=U2(:,:,2);
DWZ=U2(:,:,3);

DWXm=Um2(:,:,1);
DWYm=Um2(:,:,2);
DWZm=Um2(:,:,3);

VWX=Uw2(:,:,1);
VWY=Uw2(:,:,2);
VWZ=Uw2(:,:,3);


[void nofderiv]=size(gamma);        %nofderiv depends on rhs, thus on function setboundary5
le=(p2-p1);							%Vortex span vector


%for s=1:a
%	Lle(s)=norm(le(s,:));			%length of vortex span vector or panel span
%	lehat(s,:)=le(s,:)./Lle(s);	%
%end 


Lle=sqrt(sum(le.^2,2));
lehat(:,1)=le(:,1)./Lle;
lehat(:,2)=le(:,2)./Lle;
lehat(:,3)=le(:,3)./Lle;

%waitbar(0.9,f1,'calculating forces...');

for j=1:nofderiv
    IW(:,j,1)=DWX*gamma(:,j);
    IW(:,j,2)=DWY*gamma(:,j);
    IW(:,j,3)=DWZ*gamma(:,j);
	
	IWm(:,j,1)=DWXm*gamma(:,j);
    IWm(:,j,2)=DWYm*gamma(:,j);
    IWm(:,j,3)=DWZm*gamma(:,j);
	
	W(:,j,1)=VWX*gamma(:,j);
    W(:,j,2)=VWY*gamma(:,j);
    W(:,j,3)=VWZ*gamma(:,j);
	
	%%
	iw(:,j,1)=u2*gamma(:,j);
	
	iwm(:,j,1) = um2*gamma(:,j);
	
	wav(:,j,1) = uw2*gamma(:,j);
    
	%%
    

    
    G(:,1)=gamma(1:s1,j).*lehat(:,1);	%Aligning vorticity along panel vortex
    G(:,2)=gamma(1:s1,j).*lehat(:,2);
    G(:,3)=gamma(1:s1,j).*lehat(:,3);

    wind1=state.STW*([cos(state.psi) -sin(state.psi) 0])+[state.U state.V state.W]; %Aligning with wind

    for i=1:s1
		
		wind2(i,:)=-squeeze(W(i,j,:))';
		wind3(i,:)=-squeeze(IW(i,j,:))'-squeeze(IWm(i,j,:))';
		norm(i) = squeeze(iw(i,j,:))'+squeeze(iwm(i,j,:))'+ squeeze(wav(i,j,:));
		normal(i,:) = norm(i)'.*lattice.N(i,:);
		
		
		Wind(i,:)=wind1-squeeze(IW(i,j,:))'-squeeze(IWm(i,j,:))';
		Wind1(i,:)=Wind(i,:)-squeeze(W(i,j,:))';
        Rot(i,:)=cross((lattice.COLLOC(i,:)-geo.CG),[state.P state.Q state.R]); %Calculating rotations
    end                                   %^^^^^^^---new stuff in T131         %Thanks Luca for pointing out the error here

	Wind1=Wind1+Rot;								%Adding rotations
	
	
	%%
	%%%%%%%%
	% test %
	%%%%%%%%
	
	for kl=1:s1
		results.Windy(kl) = (Wind1(kl,1)^2+Wind1(kl,2)^2+Wind1(kl,3)^2)^(1/2);		%testpurpose - remove before release
		results.Windy2(:,kl) = [Wind1(kl,1) Wind1(kl,2) Wind1(kl,3) ];
		results.Wavy(kl) = (wind2(kl,1)^2+wind2(kl,2)^2+wind2(kl,3)^2)^(1/2);
		results.normaal(kl) = 0;
		results.induc(kl) = (wind3(kl,1)^2+wind3(kl,2)^2+wind3(kl,3)^2)^(1/2);
		results.induc2(:,kl) = [wind3(kl,1) wind3(kl,2) wind3(kl,3) ];
	end
	
	%%%%%%%%
	% test %
	%%%%%%%%
	%%

	Fprim(:,j,:)=state.rho*cross(Wind1,G);			    %Force per unit length, Kutta Joukowski
	
	
	F(:,j,1)=Fprim(:,j,1).*Lle;				%Force per panel
	F(:,j,2)=Fprim(:,j,2).*Lle;				%Force per panel
	F(:,j,3)=Fprim(:,j,3).*Lle;				%Force per panel
	
	
	C3(:,:,1)=c3(:,1)*ones(1,nofderiv);
	C3(:,:,2)=c3(:,2)*ones(1,nofderiv);
	C3(:,:,3)=c3(:,3)*ones(1,nofderiv);
	
end
results.F=F;
results.FORCE=sum(F,1);						%Total force
M=cross(C3,F,3);			                 %Moments per panel
results.M=M;
results.MOMENTS=sum(M,1);					%Summing up moments	
results.gamma=gamma;
%delete(f1);
end%FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% NEW DOWNWASH FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5five
function[dw,DW]=fastdw(lattice,mirror)
one_by_four_pi=1/(4*pi);

[psize vsize void]=size(lattice.VORTEX);


%disp('running right')
%psize=size(lattice.COLLOC,1);
lemma=ones(1,psize);

LDW=zeros(psize,psize,7,3);

mCOLLOC(:,:,1)=lattice.COLLOC(:,1)*lemma;
mCOLLOC(:,:,2)=lattice.COLLOC(:,2)*lemma;
mCOLLOC(:,:,3)=lattice.COLLOC(:,3)*lemma;

mN(:,:,1)=lattice.N(:,1)*lemma;
mN(:,:,2)=lattice.N(:,2)*lemma;
mN(:,:,3)=lattice.N(:,3)*lemma;

%for both m matrices, every column is the same

for j=1:(vsize-1)       %TP left and LP left; LP left and LP right; LP right and TP right
    
    lr1(:,:,1)=(lattice.VORTEX(:,j,1)*lemma)';
    lr1(:,:,2)=(lattice.VORTEX(:,j,2)*lemma)';
    lr1(:,:,3)=(lattice.VORTEX(:,j,3)*lemma)';
    
    lr2(:,:,1)=(lattice.VORTEX(:,j+1,1)*lemma)';
    lr2(:,:,2)=(lattice.VORTEX(:,j+1,2)*lemma)';
    lr2(:,:,3)=(lattice.VORTEX(:,j+1,3)*lemma)';
	
    %for every l matrix, every row is the same
    
    r1=lr1-mCOLLOC;
    r2=lr2-mCOLLOC;
	
	
    %by substracting every vector with every collocation point you create
    %the AIC
    %every column stands for a different vortex point, every row stands for
    %a different collocation point
	
	%every two consecutive vortex points define a vortex, by substracting 
	%lr1 and lr2 with mCOLLOC you create a matrix which contains the
	%vectors connecting the end and begin points of every vortex fillament with
	%every collocation poitn
    
    warning off
    LDW(:,:,j,:)=mega(r1,r2);
    warning on
	
	
end
LDW(find((isnan(LDW(:,:,:,:)))))=0;

DW=-squeeze(sum(LDW,3))*one_by_four_pi;     %summate the three different vortices
if mirror
	DW(:,:,3) = -DW(:,:,3);
end

dw=sum(DW.*mN,3);		%the normal component of the "induced velocity"
end


function[DW2]=mega(r1,r2)
%% First part
F1=cross(r1,r2,3);


LF1=(sum(F1.^2,3));

%	F2 is is F1 normalized
F2(:,:,1)=F1(:,:,1)./(LF1);
F2(:,:,2)=F1(:,:,2)./(LF1);
F2(:,:,3)=F1(:,:,3)./(LF1);
%clear('F1')


%% Next part

Lr1=sqrt(sum(r1.^2,3)); 
Lr2=sqrt(sum(r2.^2,3));

%normalizing r1 and r2

R1(:,:,1)=r1(:,:,1)./Lr1;
R1(:,:,2)=r1(:,:,2)./Lr1;
R1(:,:,3)=r1(:,:,3)./Lr1;

R2(:,:,1)=r2(:,:,1)./Lr2;
R2(:,:,2)=r2(:,:,2)./Lr2;
R2(:,:,3)=r2(:,:,3)./Lr2;



L1=(R2-R1);
%clear('R1','R2')



%% Third part
R0=(r2-r1);

radial_distance=sqrt((LF1./(sum(R0.^2,3))));
%because the overal function as defined by eg Drela results in a singular
%point when the perpendicular distance goes to zero a



%% combinging 2 and 3
L2=  R0(:,:,1).*L1(:,:,1)...
    +R0(:,:,2).*L1(:,:,2)...
    +R0(:,:,3).*L1(:,:,3);



%% Downwash
DW(:,:,1)=F2(:,:,1).*L2;
DW(:,:,2)=F2(:,:,2).*L2;
DW(:,:,3)=F2(:,:,3).*L2;

near=config('near');

DW2(:,:,1)=DW(:,:,1).*(1-(radial_distance<near));
DW2(:,:,2)=DW(:,:,2).*(1-(radial_distance<near));
DW2(:,:,3)=DW(:,:,3).*(1-(radial_distance<near));

%ALEC - matches with DRELA, non infinite trailing point, is this
%mathematically correct?

%add component for wavemaking disturbance here??





end



function[latticem]=mirror1(lattice,state)
latticem.COLLOC(:,1) = lattice.COLLOC(:,1);
latticem.COLLOC(:,2) = lattice.COLLOC(:,2);
latticem.COLLOC(:,3) = lattice.COLLOC(:,3);

latticem.VORTEX(:,:,1) = lattice.VORTEX(:,:,1);
latticem.VORTEX(:,:,2) = lattice.VORTEX(:,:,2);
latticem.VORTEX(:,:,3) = -lattice.VORTEX(:,:,3);

latticem.N(:,1) = lattice.N(:,1);
latticem.N(:,2) = lattice.N(:,2);
latticem.N(:,3) = lattice.N(:,3);

latticem.XYZ(:,:,1) = lattice.XYZ(:,:,1);
latticem.XYZ(:,:,2) = lattice.XYZ(:,:,2);
latticem.XYZ(:,:,3) = -lattice.XYZ(:,:,3);


end