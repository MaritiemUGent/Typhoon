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
% usage: [STATE]=statesetup(STATE)
%
% Enables a user to manipulate the STATE structure.
%
% Example:
%
%  [state]=statesetup(state);
%
% Calls:
%       questions       Contain user interface queries in string format. 
%       terror          Displays various Error messages.
%       isinptok        Checks format of imput. 
%       ISAtmosphere    Table inter polation of standard atmosphere.
%
% Author: Tomas Melin <melin@kth.se>
% Keywords: Tornado text based user interface
%
% Revision History:
%   Ghent, 2019-11-08:    Adjustments made to cater for hydrofoils.
%   Bristol, 2007-06-27:  Addition of new header. TM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state]=statesetup(state);
settings=config('startup');


q2=questions(4);
   if isempty(q2)
      q2=0;
     	terror(9)
   end
     
switch q2
%% Load file
    case 2

       ok=0;  
       while ok==0;  
          try
           	cd(settings.sdir)
   				disp(' ')   
            	ls
      	   	fname=input('Load file: ','s');
   				load(fname);
		   		cd(settings.hdir)
         	ok=1;
      	catch
         	cd(settings.hdir)
         	terror(4)
            ok=1;
      	end
   	end
   
      
%% Create State
	case 1
   		disp(' ')
         disp(' ')
         
         stepper=0;
         while stepper<99
            switch stepper
            case -1
               stepper=1;
           
            case 0
               data=input('Global Pitch (theta) [deg]: ','s');
               if isinptok(data)==1;
                  state.theta=str2num(data)*pi/180;
                  
               end
            case 1 
                data = input('Roll (phi) [deg]: ','s')
                if isinptok(data)==1
                    state.phi=str2num(data)*pi/180;
                end
            case 2
                data=input('Yaw (Psi) [deg]: ','s');
                if isinptok(data)==1;
                     state.betha=str2num(data)*pi/180;
                end
            
            case 3
				data=input('Roll angular velocity [deg/s]: ','s');
               if isinptok(data)==1;
                  state.P=str2num(data)*pi/180;
               end
   
      		case 4
               data=input('Pitch angular velocity [deg/s]: ','s');
               if isinptok(data)==1;
                  state.Q=str2num(data)*pi/180;
               end               
               
      		case 5
               data=input('Yaw angular velocity [deg/s]: ','s');
               if isinptok(data)==1;
                  state.R=str2num(data)*pi/180;
               end
               
      		case 6
               data=input('Angle of attack time derivative, (Alpha_dot), [deg/s]:','s');
               if isinptok(data)==1;
                  state.adot=str2num(data)*pi/180;
               end
                
            case 7
               data=input('Angle of sideslip time derivative, (Beta_dot), [deg/s]:','s');
               if isinptok(data)==1;
                  state.bdot=str2num(data)*pi/180;
               end
            case 8
                data = input('Speed Trough Water, (STW), [m/s]: ','s');
                if isinptok(data)==1
                    state.STW=str2num(data);
                end
            case 9
                data = input('Craft elevation (defined by distance of ref point above water surface), [m]:','s');
                if isinptok(data)==1
                    state.ELA=str2num(data);
				end
			case 10
			data = input('no of nodes in integral 1','s');
                if isinptok(data)==1
                    state.integral1=str2num(data);
				end	
			case 11
			data = input('no of nodes in integral 2','s');
                if isinptok(data)==1
                    state.integral2=str2num(data);
				end	
				
				case 12
			data = input('method used for free surface == 0','s');
                if isinptok(data)==1
                    state.fitt=0;
				end	
                       
%% Back to ordinary state            
             
            %{
                case 10
                    disp(' ')
                    disp('Caution, only use this option if you are sure, really sure, what you are doing.')
                    data=input('Apply Prandtl-Glauert Correction [0 1]: ','s');
                    if isinptok(data)==1;
                        state.pgcorr=str2num(data);
                        stepper=100;
                    end        
                    disp(' ')
                 %}
           
            
            	
               
               
            end         %caseblock
            
            if isinptok(data)==1
               stepper=stepper+1;
            elseif isinptok(data)==-1
                stepper=stepper-1;
            elseif isinptok(data)==-2
                stepper=99;
            else
            end
            
         end %whileblock
         
         %state.rho=config('rho'); %standard sealevel density

%% Save State        
  case 3 
         
     disp(' ');
     disp(' ');
     disp(' ');
     
     cd(settings.sdir)
     sfname=input('Save state as file: ','s');
     
     if isempty(sfname)==1;
         disp('+++ no name, not saved +++')
         cd(settings.hdir)
     else
    		save(sfname,'state');
          cd(settings.hdir)
          disp(' ');
          disp('*** File Saved. ****');

          
     end

     %% CHANGE ALPHA
     case(4)
     disp(' ');
     disp(' ');
     disp('Changing theta ');
     disp(strcat('Current theta is: ',num2str(state.theta*180/pi),' degrees.'))
     state.alpha=input('New theta [deg]: ')*pi/180;
     
     
     
     %ALEC - why not change other parameters here???
     
 %% exiting menu    
   case 0
      stat=1;
      return
   otherwise
      terror(9);
      stat=1;
      return    
   end
    stat=0;