 function [results,lattice,state]=equilibrium1(results,state,geo,lattice,ref)
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
% This function does a static calculation to determine where both force and
% moment balances are zero.
%
% Example:
%
%
%
% Calls:
%
%
% Author: Alec Bague <alec.bague@ugent.be>
% Keywords: Tornado core function
%
% Revision History:
%   Ghent: 11/06/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = input(' moment around y axis due to sail: ','s');
% if isinptok(data)==1
% 	state.Ms = zeros(3,1);
% 	state.Ms(2,1) = str2num(data);
% end


settings=config('startup');

state.integral1=32;
state.integral2=32;

%% initial set-up of calc

k = 1 ;

[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
results00 = solverloop6(results,1,'1',lattice,state,geo,ref);

FORCE(k,1,:)=results00.FORCE-[0; 0; geo.mass*9.81];
MOMENTS(k,1,:)=results00.MOMENTS+state.Ms;

cd(settings.odir)
	fileID = fopen('calc_equi_07072020.txt','a');
	fprintf(fileID, '%i \t %i \t %i \t %i \t \n', state.ELA, state.theta, ...
		FORCE(k,1,3),MOMENTS(k,1,2));
	fclose(fileID);
	cd(settings.hdir)

%derivative to z

state.ELA = state.ELA+0.1;

[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
results01 = solverloop6(results,1,'1',lattice,state,geo,ref);

FORCE(k,2,:)=results01.FORCE-[0; 0; geo.mass*9.81];
MOMENTS(k,2,:)=results01.MOMENTS+state.Ms;

cd(settings.odir)
	fileID = fopen('calc_equi_07072020.txt','a');
	fprintf(fileID, '%i \t %i \t %i \t %i \t \n', state.ELA, state.theta, ...
		FORCE(k,2,3),MOMENTS(k,2,2));
	fclose(fileID);
	cd(settings.hdir)

state.ELA = state.ELA-0.1;

%derivative to theta

state.theta = state.theta + 1*pi/180;

[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
results02 = solverloop6(results,1,'1',lattice,state,geo,ref);

FORCE(k,3,:)=results02.FORCE-[0; 0; geo.mass*9.81];
MOMENTS(k,3,:)=results02.MOMENTS+state.Ms;

cd(settings.odir)
	fileID = fopen('calc_equi_07072020.txt','a');
	fprintf(fileID, '%i \t %i \t %i \t %i \t \n', state.ELA, state.theta, ...
		FORCE(k,3,3),MOMENTS(k,3,2));
	fclose(fileID);
	cd(settings.hdir)

state.theta = state.theta-1*pi/180;



%%

J = [(FORCE(k,2,3)-FORCE(k,1,3))/(0.1) (FORCE(k,3,3)-FORCE(k,1,3))/(1*pi/180); ...
	(MOMENTS(k,2,2)-MOMENTS(k,1,2))/(0.1) (MOMENTS(k,3,2)-MOMENTS(k,1,2))/(1*pi/180)];


A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);


%this step limits the maximum  step size

if (abs(A(k,1))>0.1)
	A(k,1) = sign(A(k,1))*0.1;
end

if (abs(A(k,2))>0.0175)
	A(k,2) = sign(A(k,2))*0.0175;
end

inc=1;

cd(settings.odir)
fileID = fopen('calc_equi_07072020.txt','a');
fprintf(fileID, 'update : %i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
	FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
fclose(fileID);
cd(settings.hdir)



while (k<100 && inc>1.8e-07)
	angle_c=1e-02*pi/180;
	k=k+1;
	if (abs(FORCE(k-1,1,3))>(0.1*9.81*geo.mass) && abs(MOMENTS(k-1,1,2))>abs(state.Ms(2,1)*0.1))
		if ((abs(A(k-1,1))>1e-03) && (abs(A(k-1,2))>angle_c))
			
			%derivative to z
			
			state.ELA = state.ELA+A(k-1,1);
			
			[lattice,ref,geo]=fLattice_setup2(geo,state,0);
			stat=1;
			results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
			
			FORCE(k,2,:)=results01.FORCE-[0; 0; geo.mass*9.81];
			MOMENTS(k,2,:)=results01.MOMENTS+state.Ms;
			
			
			
			%derivative to theta
			
			state.theta = state.theta + A(k-1,2);
			
			[lattice,ref,geo]=fLattice_setup2(geo,state,0);
			stat=1;
			results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
			
			FORCE(k,3,:)=results02.FORCE-[0; 0; geo.mass*9.81];
			MOMENTS(k,3,:)=results02.MOMENTS+state.Ms;
			
			FORCE(k,1,:)=FORCE(k,3,:);
			MOMENTS(k,1,:)=MOMENTS(k,3,:);
			
			%%
			
			J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(A(k-1,1)) (FORCE(k,3,3)-FORCE(k-1,1,3))/(A(k-1,2)); ...
				(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(A(k-1,1)) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(A(k-1,2))];
			
			
			A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
			
			if (abs(A(k,1))>0.1)
				A(k,1) = sign(A(k,1))*0.1;
			end
			
			if (abs(A(k,2))>0.0175)
				A(k,2) = sign(A(k,2))*0.0175;
			end
			
			inc = abs(A(k,1))*abs(A(k,2));
			
			cd(settings.odir)
			fileID = fopen('calc_equi_07072020.txt','a');
			fprintf(fileID, '%i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
				FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
			fclose(fileID);
			cd(settings.hdir)
		else
			if (abs(A(k-1,1))>1e-03)
				
				
				%derivative to z, only movement
				
				state.ELA = state.ELA+A(k-1,1);
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				FORCE(k,2,:)=results02.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,2,:)=results02.MOMENTS+state.Ms;
				
				FORCE(k,1,:)=FORCE(k,2,:);
				MOMENTS(k,1,:)=MOMENTS(k,2,:);
				
				%derivative to theta
				
				state.theta = state.theta + 0.5*pi/180;
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				state.theta = state.theta - 0.5*pi/180;
				
				FORCE(k,3,:)=results01.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,3,:)=results01.MOMENTS+state.Ms;
				
				
				J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(A(k-1,1)) (FORCE(k,3,3)-FORCE(k-1,1,3))/(0.5*pi/180); ...
					(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(A(k-1,1)) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(0.5*pi/180)];
				
				A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
				
				if (abs(A(k,1))>0.1)
					A(k,1) = sign(A(k,1))*0.1;
				end
				
				if (abs(A(k,2))>0.0175)
					A(k,2) = sign(A(k,2))*0.0175;
				end
				
				inc = abs(A(k,1))*abs(A(k,2));
				
				cd(settings.odir)
				fileID = fopen('calc_equi_07072020.txt','a');
				fprintf(fileID, '%i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
					FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
				fclose(fileID);
				cd(settings.hdir)
			else
				if (abs(A(k-1,2))>angle_c)
					
					
					%derivative to thet, only movement
					
					state.theta = state.theta+A(k-1,2);
					
					[lattice,ref,geo]=fLattice_setup2(geo,state,0);
					stat=1;
					results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
					
					FORCE(k,3,:)=results02.FORCE-[0; 0; geo.mass*9.81];
					MOMENTS(k,3,:)=results02.MOMENTS+state.Ms;
					
					FORCE(k,1,:)=FORCE(k,3,:);
					MOMENTS(k,1,:)=MOMENTS(k,3,:);
					
					%derivative to z
					
					state.ELA = state.ELA + 0.05;
					
					[lattice,ref,geo]=fLattice_setup2(geo,state,0);
					stat=1;
					results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
					
					state.ELA = state.ELA - 0.05;
					
					FORCE(k,3,:)=results01.FORCE-[0; 0; geo.mass*9.81];
					MOMENTS(k,3,:)=results01.MOMENTS+state.Ms;
					
					J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(A(k-1,1)) (FORCE(k,3,3)-FORCE(k-1,1,3))/(0.5*pi/180); ...
						(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(A(k-1,1)) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(0.5*pi/180)];
					
					
					A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
					
					if (abs(A(k,1))>0.1)
						A(k,1) = sign(A(k,1))*0.1;
					end
					
					if (abs(A(k,2))>0.0175)
						A(k,2) = sign(A(k,2))*0.0175;
					end
					
					inc = abs(A(k,1))*abs(A(k,2));
					
					cd(settings.odir)
					fileID = fopen('calc_equi_07072020.txt','a');
					fprintf(fileID, '%i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
						FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
					fclose(fileID);
					cd(settings.hdir)
				else
					%% too small both ways, re-initialize
					%% initial set-up of calc
					
					
					
					%derivative to z
					
					state.ELA = state.ELA+0.05;
					
					[lattice,ref,geo]=fLattice_setup2(geo,state,0);
					stat=1;
					results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
					
					FORCE(k,2,:)=results01.FORCE-[0; 0; geo.mass*9.81];
					MOMENTS(k,2,:)=results01.MOMENTS+state.Ms;
					
					cd(settings.odir)
					fileID = fopen('calc_equi_07072020.txt','a');
					fprintf(fileID, '%i \t %i \t %i \t %i \t \n', state.ELA, state.theta, ...
						FORCE(k,2,3),MOMENTS(k,2,2));
					fclose(fileID);
					cd(settings.hdir)
					
					state.ELA = state.ELA-0.05;
					
					%derivative to theta
					
					state.theta = state.theta + 0.5*pi/180;
					
					[lattice,ref,geo]=fLattice_setup2(geo,state,0);
					stat=1;
					results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
					
					FORCE(k,3,:)=results02.FORCE-[0; 0; geo.mass*9.81];
					MOMENTS(k,3,:)=results02.MOMENTS+state.Ms;
					
					FORCE(k,1,:) = FORCE(k-1,1,:);
					MOMENTS(k,1,:) = FORCE(k-1,1,:);
					
					cd(settings.odir)
					fileID = fopen('calc_equi_07072020.txt','a');
					fprintf(fileID, '%i \t %i \t %i \t %i \t \n', state.ELA, state.theta, ...
						FORCE(k,3,3),MOMENTS(k,3,2));
					fclose(fileID);
					cd(settings.hdir)
					
					state.theta = state.theta-0.5*pi/180;
					
					
					
					%%
					
					J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(0.05) (FORCE(k,3,3)-FORCE(k-1,1,3))/(0.5*pi/180); ...
						(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(0.05) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(0.5*pi/180)];
					
					
					A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
					
					
					%this step limits the maximum  step size
					
					if (abs(A(k,1))>0.1)
						A(k,1) = sign(A(k,1))*0.1;
					end
					
					if (abs(A(k,2))>0.0175)
						A(k,2) = sign(A(k,2))*0.0175;
					end
					
					inc=1;
					
					cd(settings.odir)
					fileID = fopen('calc_equi_07072020.txt','a');
					fprintf(fileID, 'update : %i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
						FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
					fclose(fileID);
					cd(settings.hdir)
					
					
				end
			end
		end
	else
		%% if we are close enough to the solution, take alternating steps
		k2 = mod(k,2);
		switch k2
			case 0
				%% angle step
				
				
				%derivative to theta
				
				state.theta = state.theta + A(k-1,2);
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				FORCE(k,3,:)=results02.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,3,:)=results02.MOMENTS+state.Ms;
				
				
				
				%% side step for derivative to z
				
				state.ELA = state.ELA+0.05;
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				FORCE(k,2,:)=results01.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,2,:)=results01.MOMENTS+state.Ms;
				
				FORCE(k,1,:) = FORCE(k,3,:);
				MOMENTS(k,1,:) = MOMENTS(k,3,:);
				
				state.ELA = state.ELA -0.05;
				
				%%
				
				J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(0.05) (FORCE(k,3,3)-FORCE(k-1,1,3))/(A(k-1,2)); ...
					(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(0.05) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(A(k-1,2))];
				
				
				A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
				
				if (abs(A(k,1))>0.1)
					A(k,1) = sign(A(k,1))*0.1;
				end
				
				if (abs(A(k,2))>0.0175)
					A(k,2) = sign(A(k,2))*0.0175;
				end
				
				inc = abs(A(k,1))*abs(A(k,2));
				
				cd(settings.odir)
				fileID = fopen('calc_equi_07072020.txt','a');
				fprintf(fileID, '%i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
					FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
				fclose(fileID);
				cd(settings.hdir)
			case 1
				%% elavation step
				%derivative to z
				
				state.ELA = state.ELA+A(k-1,1);
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results02 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				FORCE(k,2,:)=results02.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,2,:)=results02.MOMENTS+state.Ms;
				
				
		
				
				%% side step for derivative to theta
				
				state.theta = state.theta+0.5*pi/180;
				
				[lattice,ref,geo]=fLattice_setup2(geo,state,0);
				stat=1;
				results01 = solverloop6(results,1,'1',lattice,state,geo,ref);
				
				FORCE(k,3,:)=results01.FORCE-[0; 0; geo.mass*9.81];
				MOMENTS(k,3,:)=results01.MOMENTS+state.Ms;
				
				FORCE(k,1,:) = FORCE(k,2,:);
				MOMENTS(k,1,:) = MOMENTS(k,2,:);
				
				state.theta = state.theta-0.5*pi/180;
				
				%%
				
				J = [(FORCE(k,2,3)-FORCE(k-1,1,3))/(A(k-1,1)) (FORCE(k,3,3)-FORCE(k-1,1,3))/(0.5*pi/180); ...
					(MOMENTS(k,2,2)-MOMENTS(k-1,1,2))/(A(k-1,1)) (MOMENTS(k,3,2)-MOMENTS(k-1,1,2))/(0.5*pi/180)];
				
				
				A(k,:) = transpose(J^(-1) * [-FORCE(k,1,3); -MOMENTS(k,1,2)]);
				
				if (abs(A(k,1))>0.1)
					A(k,1) = sign(A(k,1))*0.1;
				end
				
				if (abs(A(k,2))>0.0175)
					A(k,2) = sign(A(k,2))*0.0175;
				end
				
				inc = abs(A(k,1))*abs(A(k,2));
				
				cd(settings.odir)
				fileID = fopen('calc_equi_07072020.txt','a');
				fprintf(fileID, '%i \t %i \t %i \t %i \t %i \t %i \n', state.ELA, state.theta, ...
					FORCE(k,1,3),MOMENTS(k,1,2), A(k,1), A(k,2));
				fclose(fileID);
				cd(settings.hdir)
		end
	end
	if ((abs(FORCE(k,1,3))+abs(MOMENTS(k,1,2)))<(0.05*geo.mass*9.81))
		kk=k;
		k=100;
	end
	
	
end

results = results02;

state.theta = state.theta + A(kk-1,2);
state.ELA = state.ELA+A(kk-1,1);

results.stat_eq = [state.theta, state.ELA]





save('calc_equi');

end








