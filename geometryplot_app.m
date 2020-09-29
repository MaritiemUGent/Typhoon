function[app_fig,h2,a3]=geometryplot_app(lattice,geo,ref,state,option_no);

[t q2 q3]=size(lattice.VORTEX);
% %%
% figure(1)
% %Changing variables to plot only partition outline
% g2=geo;
% g2.nx=double(g2.nx>0);
% g2.ny=double(g2.ny>0); 
% g2.fnx=double(g2.fnx>0); 
% s2.STW=1;
% s2.theta=0;
% s2.phi=0;
% s2.psi=0;
% s2.P=0;
% s2.Q=0;
% s2.R=0;
% s2.ELA=0;
% s2.rho=997;
% s2.pgcorr=0;
% 
% [l2,ref]=fLattice_setup2(g2,s2,1);
% 
% g=fill3(l2.XYZ(:,:,1)',l2.XYZ(:,:,2)',l2.XYZ(:,:,3)','w');
% set(g,'LineWidth',2);
% view([0,90]);
% axis equal,hold on
% xlabel('Aircraft body x-coordinate')
% ylabel('Aircraft body y-coordinate')
% zlabel('Aircraft body z-coordinate')
% title('3D wing and partition layout')
% grid on
% 
% h=line([ref.mac_pos(1) ref.mac_pos(1)+ref.C_mac],[ref.mac_pos(2) ref.mac_pos(2)],[ref.mac_pos(3) ref.mac_pos(3)]);
% set(h,'LineWidth',5);
% 
% a=plot3(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3)+state.ELA,'r+');
% set(a,'MarkerSize',15,'linewidth',3);
% a=plot3(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3)+state.ELA,'ro');
% set(a,'MarkerSize',15,'linewidth',3);
% 
% a=plot3(geo.CG(1),geo.CG(2),geo.CG(3),'ko');
% set(a,'MarkerSize',15,'linewidth',3);
% a=plot3(geo.CG(1),geo.CG(2),geo.CG(3),'kx');
% set(a,'MarkerSize',15,'linewidth',3);
% 
% %plotting legend
% L=gca;
% set(L,'Position',[0.1 0.1 0.6 0.8]);
% axes('position',[0.75 0.6 0.2 0.2]);
% axis([0 1 0 1])
% hold on
% h=line([0.1 0.4],[1 1]);
% set(h,'LineWidth',6);
% 
% a=plot(0.25,0.66,'r+');
% set(a,'MarkerSize',15,'linewidth',3);
% a=plot(0.25,0.66,'ro');
% set(a,'MarkerSize',15,'linewidth',3);
% 
% a=plot(0.25,0.33,'kx');
% set(a,'MarkerSize',15,'linewidth',3);
% a=plot(0.25,0.33,'ko');
% set(a,'MarkerSize',15,'linewidth',3);
% 
% text(0.5,1,'MAC');
% text(0.5,0.66,'ref point')
% text(0.5,0.33,'c.g.')
% 
% 
% axis off
% 
% 
% try
%     B=lattice.XYZ(:,:,3);           %Check if geometry is present.
% catch
%     terror(10);
%     return
% end

%%
h3=[]; h2=[]; a3=[];
if (option_no==1)
	clf
app_fig=figure(50);
set(gca,'Position',[0 0 1 1]);
set(app_fig,'Color','White');
set(app_fig,'Visible','off')



    h3=plot3(lattice.XYZ_o(:,:,1)',lattice.XYZ_o(:,:,2)',lattice.XYZ_o(:,:,3)','--r');
	hold on
	h2=plot3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)','k');
    set(gca,'Position',[0.05 0.05 0.9 0.90]);
    axis equal,hold on,axis off, view([25 15]);
    title('ISO')
    grid on
    h=line([ref.mac_pos(1) ref.mac_pos(1)+ref.C_mac],[ref.mac_pos(2) ref.mac_pos(2)],[ref.mac_pos(3) ref.mac_pos(3)]);
    set(h,'LineWidth',5);
    a=plot3(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3),'r+');
    set(a,'MarkerSize',15);
    a=plot3(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3),'ro');
    set(a,'MarkerSize',15);
	xLimits = get(gca,'XLim')+[-0.5 0.5];  yLimits = get(gca,'YLim')+[-0.5 0.5];  zLimits = get(gca,'ZLim')+[-0.5 0.5];  
	[X,Y]=meshgrid(xLimits(1):0.2:xLimits(2),yLimits(1):0.2:yLimits(2));
	Z=Y.*0-state.ELA;
	s=surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.8,'FaceColor','b');
	axis([xLimits, yLimits, zLimits]);
	quiver3(0,0,0,1,0,0);
	quiver3(0,0,0,0,1,0);
	quiver3(0,0,0,0,0,1);
	

	
	

   
    %a=text(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3)+state.ELA,'Reference point');
    
    %try
    %        a=plot3(geo.CG(1),geo.CG(2),geo.CG(3),'kx');
    %        set(a,'MarkerSize',15);
    %         a=plot3(geo.CG(1),geo.CG(2),geo.CG(3),'ko');
    %         set(a,'MarkerSize',15);
    %         a=text(geo.ref_point(1),geo.ref_point(2),geo.ref_point(3)+state.ELA,'Center of gravity');
    %         
    %end
%%    
end
if (option_no==2)
app_fig=figure(2)
plot3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)','m')
hold on
grid on  
for s=1:(t)	
   w=0;
   for u=1:q2-1
      w=w+1;
      VX(w,:)=[lattice.VORTEX(s,u,1) lattice.VORTEX(s,u+1,1)];
      VY(w,:)=[lattice.VORTEX(s,u,2) lattice.VORTEX(s,u+1,2)];
      VZ(w,:)=[lattice.VORTEX(s,u,3) lattice.VORTEX(s,u+1,3)];
	  
	  
	end   
   	rc=lattice.COLLOC(s,:);
   	A=rc+lattice.N(s,:);					%Check routine
      x=[rc(1) A(1)];				%Calculating normals
      y=[rc(2) A(2)];				         
      z=[rc(3) A(3)];
         
      NORMALS(s,:,1)=x;				%saving normals
      NORMALS(s,:,2)=y;
      NORMALS(s,:,3)=z;
      
      plot3(VX,VY,VZ,'r.-.')
 
   end
xlabel('Body x-coord')
ylabel('Body y-coord')
zlabel('Body z-coord')
title('3D wing configuration, vortex and wake layout.')
axis equal 
%%
a3=figure(3);
plot3(lattice.XYZ(:,:,1)',lattice.XYZ(:,:,2)',lattice.XYZ(:,:,3)','k')
hold on

for s=1:(t)
   plot3(lattice.COLLOC(s,1),lattice.COLLOC(s,2),lattice.COLLOC(s,3),'g*')
   plot3(NORMALS(s,:,1),NORMALS(s,:,2),NORMALS(s,:,3),'r:');
end

axis equal
xlabel('Body x-coordinate')
ylabel('Body y-coordinate')
zlabel('Body z-coordinate')
title('3D panels, collocation points and normals.')


grid on
end

end