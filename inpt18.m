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
% usage: [GEO] = INPT (GEO)
%
% Invokes the geometry manipulation script of Tornado to allow 
% the user to load, save, edit or input a new geometry structure.
% Geo is the structure containing all necessary geometrical 
% definitions.
%
% Example:
%
%  [geo]=inpt18(geo);
%
% Calls:
%       questions                   Contain user interface queries in string format. 
%       terror                      Displays various Error messages.
%       isinptok                    Checks format of imput. 
%       tedit                       Edit geo interface.
%       f_version_geo_transform     converts old geometry files to new standard.
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado text based user interface
%
% Revision History:
%   Ghent,    2019 11 08:  Adjusted to cater for Hydrofoils
%                          Added the possibility for curved foils and
%                          removed blunt body type. Still some nomenclature
%   Bristol,  2007 06 27:  Addition of new header. TM.
%   Stockholm,2004 01 18:  Changed all geometry variables to new format. TM.
%                          Added a "convert geometry file" function. TM.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [geo]=inpt18(geo)

 settings=config('startup'); %setting directories


loop1=1;
while loop1==1;
q=questions(2);
stat=0;

if isempty(q)
   q=5; %will go to otherwise section
end


switch q
    
    %%
	case 1 %enter new geometry
		
		 ok1=0;
        geo=[];
        geo.flapped=0;
		
        disp(' ')
        disp(' ************ ')
        data = input(' what is the mass of the craft [kg]: ','s');
        if isinptok(data)==1
            geo.mass=str2num(data);
		end
		
		disp(' ')
        disp(' ************ ')
		disp(' moments of inertia')
        data = input(' Ixx: ','s');
        if isinptok(data)==1
            geo.I(1,1)=str2num(data);
		end
        data = input(' Ixy: ','s');
        if isinptok(data)==1
            geo.I(1,2)=str2num(data);
		end
		data = input(' Ixz: ','s');
        if isinptok(data)==1
            geo.I(1,3)=str2num(data);
        end
		data = input(' Iyx: ','s');
        if isinptok(data)==1
            geo.I(2,1)=str2num(data);
        end
		data = input(' Iyy: ','s');
        if isinptok(data)==1
            geo.I(2,2)=str2num(data);
        end
		data = input(' Iyz: ','s');
        if isinptok(data)==1
            geo.I(2,3)=str2num(data);
		end
		data = input(' Izx: ','s');
        if isinptok(data)==1
            geo.I(3,1)=str2num(data);
		end
		data = input(' Izy: ','s');
        if isinptok(data)==1
            geo.I(3,2)=str2num(data);
		end
		data = input(' Izz: ','s');
        if isinptok(data)==1
            geo.I(3,3)=str2num(data);
        end
		
		
       
        
        while ok1==0
    	  data=input('Number of Foils: ','s');
            if isinptok(data)==1;   %ALEC - checks if the input is within range
         	    geo.nwing=str2num(data);
                ok1=1;
            end
     	end
     
  		for s=1:geo.nwing			%Stepping over the wings.
   		%home	
  		 	disp('""""""""""""""""""""""""""""""""""""""""""""""""""""""')
			disp('"                                                     ')
            disp(strcat('" Data regarding wing number :', num2str(s),'  '))
   		    disp('"                                                     ')
			disp('""""""""""""""""""""""""""""""""""""""""""""""""""""""')
            disp(' ');
         
         ok2=0;
         while ok2==0
         data=input('Number of semispanwise partitions for this wing: ','s');
         	if isinptok(data)==1;
          		geo.nelem(s)=str2num(data);
        	 	ok2=1;
            end
         end 
         
   		t=0;
		for s2=1:geo.nelem(s)       %stepping over the partitions in the wings
         	t=t+1;
         
         	disp(' ')
  		 	   disp('""""""""""""""""""""""""""""""""""""""""""""""""""""""')
			   disp('"                                                    ')
               disp(strcat('" Data regarding partition number:  ',num2str(t),('                  ')))
            
         	disp(' ')   
            
            if s==1
               disp(' ')
               disp('b - Back one question. ')
               disp('q - Abort input sequence. ')
               disp(' ')
            end
            
           	
            
            
            stepper=-6;
            while stepper<99
            	
                 switch stepper
                 	
                     
                  case -7
                  	stepper=1;
                  case -6                   %CG data added TM20070329
                       if s==1
                            if t==1;
                                data=input('Center of gravity x-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.CG(1)=str2num(data);
                                end
                            end
                        end   
                         
                   case -5
                       if s==1
                            if t==1;
                                data=input('Center of gravity y-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.CG(2)=str2num(data);
                                end
                            end
                        end 
                   case -4
                       if s==1
                            if t==1;
                                data=input('Center of gravity z-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.CG(3)=str2num(data);
                                end
                            end
                        end
                    
                    
                    
                    case -3
                       if s==1
                            if t==1;
                                data=input('Reference point x-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.ref_point(1)=str2num(data);
                                end
                            end
                        end   
                         
                   case -2
                       if s==1
                            if t==1;
                                data=input('Reference point y-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.ref_point(2)=str2num(data);
                                end
                            end
                        end 
                   case -1
                       if s==1
                            if t==1;
                                data=input('Reference point z-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.ref_point(3)=str2num(data);
                                end
                            end
                        end                    
                    
                   case 0
                      
                    	if t==1;
                			data=input('Is the wing mirrored in the xz-plane [1 0]: ','s');
                         if isinptok(data)==1;
                            geo.symetric(s)=str2num(data);
                         end
                      end
                  case 1            %ALEC this is to determine the second wing's location, changed to also give location to first wing
                        %if s~=1
                            if t==1;
                                data=input('Apex x-coordinate: ','s');
                                if isinptok(data)==1;
                        	        geo.startx(s)=str2num(data);
                                end
                            end
                        %else
                        %    geo.startx(s)=0;
                        %end
                     
                  case 2
                        %if s~=1
                            if t==1;
                                data=input('Apex y-coordinate: ','s');
                                if isinptok(data)==1;
                                    geo.starty(s)=str2num(data);
                                end
                            end
                        %else
                        %    geo.starty(s)=0;
                        %end
                            
                  case 3
                        %if s~=1
                            if t==1;
                                data=input('Apex z-coordinate: ','s');
                                if isinptok(data)==1;
                                    geo.startz(s)=str2num(data);
                                end
                            end
                        %else
                            %geo.startz(s)=0;    
                        %end
                        

                            
                  case 4
                     if t==1;
                        data=input('Root chord: ','s');
                        if isinptok(data)==1;
                           geo.c(s)=str2num(data);
                        end
                     end
                  case 5
                     if t==1;      
                        ok=0;  
                        while ok==0;
                            try
                                cd(settings.afdir)
      

                                disp(' ')
                                disp('****************************** ')
                                disp(' AVAILABLE AIRFOILS: ')
                                ls
                                disp('****************************** ')
                                disp(' ')
                                disp('Enter profile filename from the list above (ex CLARKY.DAT) ')
                                disp('OR any NACA four digits series numer (ex: 2412)')
                                disp('0 (zero) for a flat plate. ')
                                disp(' ')
                                data=input('Base chord airfoil: ','s');
                                
                                try
                                    foo=str2num(data);
                                    cd(settings.hdir)
                                catch   
                                    load(data)  %Testload to see that the file exists
                                    cd(settings.hdir)
                                end
     
                                   geo.foil(s,t,1)={data};
                                   data=0;
                                   
                                ok=1;
                            catch
                                cd(settings.hdir)
                                terror(4)
                            end
                        end
                     end 
                  
                	 case 6   
                     if t==1;
                        data=input('Number of panels chord wise: ','s');
                        if isinptok(data)==1;
                           geo.nx(s,1)=str2num(data);
                        end
                        
                     end


                    
                     
                  case 7
                     if t==1;
                          data=input('Partition dihedral [deg]: ','s');
                          if isinptok(data)==1;
                             geo.dihed(s,t)=str2num(data)*pi/180;      %adjusted to make partitions tangent to previous panels
                          end
                     else
                         if geo.curved(s,t-1)==0
        		         data=input('Partition dihedral [deg]: ','s');
                          if isinptok(data)==1;
                             geo.dihed(s,t)=geo.dihed(s,t-1)+str2num(data)*pi/180;      %adjusted to make partitions tangent to previous panels
                          end
                         else
                             data=input('Partition dihedral [deg]: ','s');
                          if isinptok(data)==1;
                             geo.dihed(s,t)=geo.dihed(s,t-1)-geo.b(s,t-1)/geo.radius(s,t-1)+str2num(data)*pi/180;      %adjusted to make partitions tangent to previous panels
                          end
                         end
                     end
                     
                     
                     case 8        %currently symmetric pitch control
                     %different from sweep angle, as it affects the entire
                     %axis system of the foil as a whole, inserting a theta
                     %here affects also subsequent partitions.
                     
                      if t==1 %first partition only
                        
                            disp(' ')
                            disp(' *************** ')
                            disp(' Picth control of foil... ')
                            disp(' Rotated around normal direction  ')
                            data = input(' Pitch of this foil (AoA) [deg]: ','s');
                            if isinptok(data)==1
                                geo.foiltheta(s,t)=str2num(data)*pi/180;
                            end
                      end
                     
                     case 9
                     %if s==1	% First wing base chord is reference for twist 
              		%		geo.TW(s,1,1)=0;
                     %else
                            if t==1     %first partition only
                                data=input('Base chord twist [deg]: ','s');
                                if isinptok(data)==1;
                                    geo.TW(s,1,1)=str2num(data)*pi/180;
                                    geo.BTW(s)=str2num(data)*pi/180;
                                end                          
                            end
                            
                        
                     %end  % we take this away to also allow AoA for first
                     %wing
                     
                  case 10
         				data=input('Number of panels semi-span wise: ','s');
                     if isinptok(data)==1;
                        geo.ny(s,t)=str2num(data);
                     end
                     
                  case 11
                     data=input('Span of partition: ','s');
                     if isinptok(data)==1;
                        geo.b(s,t)=str2num(data);
                     end
                     
                  case 12
                     data=input('Taper ratio: ','s');
                     if isinptok(data)==1;
                        geo.T(s,t)=str2num(data);
                     end
                     
                  case 13
                             
                        ok=0;  
                        while ok==0;
                            try
                                cd(settings.afdir)
                                disp(' ')
                                disp('****************************** ')
                                disp(' AVAILABLE AIRFOILS: ')
                                ls
                                disp('****************************** ')
                                disp(' ')
                                disp('Enter profile filename from the list above (ex CLARKY.DAT) ')
                                disp('OR any NACA four digits series numer (ex: 2412)')
                                disp('0 (zero) for a flat plate. ')
                                disp(' ')
                                data=input('Tip chord airfoil: ','s');
                                
                                try
                                    foo=str2num(data);
                                    cd(settings.hdir)
                                catch   
                                    load(data)  %Testload to see that the file exists
                                    cd(settings.hdir)
                                end
     
                                   geo.foil(s,t,2)={data};
                                   data=0;
                                ok=1;
                            catch
                                cd(settings.hdir)
                                terror(4)
                            end
                        end
                      

                  case 14
					 disp('		SWEEP does not combine with curve...')
					 disp('		will be set to zero if (curve==1).')
					 disp('		')
					 %combining curve with sweep is work in progress
                     data=input('Quarter chord line sweep [deg]: ','s');
                     if isinptok(data)==1;
                        geo.SW(s,t)=str2num(data)*(pi/180);
                     end
                     
                  case 15   
 					  
					  disp('	please consider that twist and curve cannot be combined:')
					  disp('	twist will be set to zero if (curved==1)')
					  data=input('Outboard twist [deg]: ','s');
					  
                      if isinptok(data)==1;
                         geo.TW(s,t,2)=str2num(data)*(pi/180);
                      end

                  case 16
                      disp(' ')
                      disp(' *************** ')
                      disp(' Available mesh distribution types:')
                      disp('   [1] Linear')
                      disp('   [2] Spanwise half-cosine')
                      disp('   [3] Spanwise half-cosine, chordwise cosine')
                      disp('   [5] Spanwise cosine')
                      disp('   [6] Chordwise cosine')% (Added 22/08/2008 AT)
                      disp('   [7] 3:rd order centerpacking. (Not for wings)') 
                      disp(' ')
 					  data=input('Mesh type: ','s');
                      if isinptok(data)==1;
                         geo.meshtype(s,t)=str2num(data);
                      end
                      
                      
                      
                 
                     
                 case 17
                      disp(' ')
                      disp(' *************** ')
                      disp(' Curvature:')
 					  data=input(' Is the foil partition curved [1 0] ','s');
                      if isinptok(data)==1;
                         geo.curved(s,t)=str2num(data);
                      end
                      
                 case 18
                      if geo.curved(s,t)==1
                          disp('What is the radius: ')
                          disp('     Foil is assumed tangent to previous partition (or according present dihedral angle)')
                          disp('     ')
                          data = input(' Radius [m]: ','s');
                      end
                      if isinptok(data)==1
                          geo.radius(s,t)=str2num(data); 
                          
                      end
                     
                      
                         
                      
                    
                   case 19
                      %data=input('Is partition flapped [1 0]:','s');
                      %for now no flaps
                      data='0';
                      if isinptok(data)==1;
                         geo.flapped(s,t)=str2num(data);
                      end
                      if geo.flapped==0;
                         stepper=100;
                      end
                   case 20    
                  	if geo.flapped(s,t)==1;
            				data=input('Flap chord in fraction of local chord (0..1): ','s');
                        if isinptok(data)==1;
                           geo.fc(s,t)=str2num(data);
                        end 
                     end
                                       
                  case 21
                     if geo.flapped(s,t)==1;
                        data=input('Number of chord wise panels on flap: ','s'); 
                        if isinptok(data)==1;
                           geo.fnx(s,t)=str2num(data);
                        end
                        
                     end

                     case 22
                     if and(geo.symetric(s),geo.flapped(s,t));
              				data=input('Do control surfaces deflect symmetrically [1 0]:','s');
                          if isinptok(data)==1;
                             geo.fsym(s,t)=str2num(data);
                          end
                          
                     else
               			geo.fsym(s,t)=0;
                     end
                     
                   
                  end %caseblock
                  
                  if isinptok(data)==1
                     stepper=stepper+1;
                  elseif isinptok(data)==-1
                     stepper=stepper-1;
                  elseif isinptok(data)==-2
                     stepper=99;
                     return
				  else
					  
                  end
                                        
                  
               end% while partitionblock
               
               if geo.flapped(s,t)==0
                  stepper=100;
               	geo.fc(s,t)=0;
                geo.fnx(s,t)=0;
               	geo.fsym(s,t)=0;
         		end
         		geo.nx(s,t)=geo.nx(s,1);
            	%dihed(s,t)=dihed(s,1);
				
				if geo.curved(s,t)		%combination twist/sweep and curve not supported
					geo.TW(s,t,2)=0;
					geo.SW(s,t)=0;
				end
            	if t~=1
            		geo.TW(s,t,1)=geo.TW(s,t-1,2);	%continious twist
            		geo.foil(s,t,1)=geo.foil(s,t-1,2); %continious camber
               end
             
      	end%partitionblock
      end%wingblock
      
      geo.flap_vector=zeros(size(geo.flapped));
      if isinptok(data)~=-2
         loop1=0;				%loop if aborted
      end
%%
      
case 2
    geo=[];
    lattice=[];
  	try
        cd(settings.acdir)
        ls
     	fname1=input('Load file: ','s');
        load(fname1);  
        cd(settings.hdir)
      
        %if isempty(geo.flap_vector)                   %remove this passus in version current +2
        %    geo.flap_vector=zeros(size(geo.flapped)); %fix for old geometry files
        %end
%         if (tph_version)~=1
%             disp(' ')
%             disp('*** WARNING! *** ')
%             disp(' ')
%             disp('Old geometry standard, NO geometry loaded. ')
%             disp('Please use geometry conversion function and reload.')
%             disp(' ')
%             return
%         end
      
      loop1=0;
    catch
        cd(settings.hdir)
        terror(4)        
    end
      %%
case(3)   
   [geo]=tedit(geo);        %ALEC - not adjusted
   %%
case(4)
   cd(settings.acdir)
      ls
      tph_version=1;
      fname2=input('Save file as: ','s');
      if isempty(fname2)==1;
         disp('+++ no name, not saved +++')
         cd(settings.hdir)
      else
      		try
        	  save(fname2,'geo','tph_version') ;  
              disp('*** File Saved. ****');
              cd(settings.hdir)
      		catch
        			cd(settings.hdir)
        			terror(5)   
         	end
      end
   
 %{ 
case(5)
     %f_version_geo_transform;
     geo=bodyinput(geo);
     
     %}
   case (0)
      loop1=0;
   
	otherwise
   	terror(9);	
	end
end

    