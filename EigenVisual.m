function [] = EigenVisual(stability)

stab = stability.StabAlec;
%  stab = [-0.27 -2.21 -4.99 -0.5 0;
%  	0 0 0 1 0;
%  	0 0 0 0 1;
%  	-2.24 -57.2 -156.21 -15.62 -6.51;
%  	0.01 -6.42 -78.59 -7.86 -14.75];

[V,W] = eig(stab);


%%		Eigenvalue plot

 [~,idx]=sort(diag(real(W)),'descend');
 W=diag(W(idx,idx));
 M = max(imag(W));

figure(1)

[s1,s2]=size(W);
for i=1:s1
	plot(real(W(i)),imag(W(i)),'marker','o','MarkerFaceColor', 'r','MarkerEdgeColor','r') 
	hold on
	txt=strcat('\lambda',num2str(i));
	text(real(W(i))-0.75,imag(W(i))+0.75,txt,'Color', 'r')
end
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('Re')
ylabel('Im')
hi = real(W(1))+5;
lo = real(W(s1))-5;

xlim([-20 hi]);

lo = -M-2.5; hi = M+2.5;
ylim([lo hi])
set(gcf,'color','w');
set(gca,'box','off');




  		savefig('output\eigenvalue');


%%		Eigenvector plot

V(:,:) = V(:,idx);
last = 0;

for i=1:s1
	if (imag(W(i))==0)
		figure(i+1)
		x=categorical({'\Delta u'; '\Delta z'; '\Delta \theta';'\Delta w'; '\Delta q'});
		y=V(:,i)';
		bar(x,y)
		
		set(gcf,'color','w');
		set(gca,'box','off');
		
		txt=strcat('output\eigenvector',num2str(i));
		
			savefig(txt);
		
	else
		if (last~=real(W(i)))
			last = real(W(i));
			
			x=["\Delta u"; "\Delta z"; "\Delta \theta";"\Delta w"; "\Delta q"];
			
			Mr = max(real(V(:,i)));
			Mi = max(imag(V(:,i)));
			mr = min(real(V(:,i)));
			mi = min(imag(V(:,i)));
			
			figure(i+1)
			for j=1:s1
				quiver(diag(zeros(s1)),diag(zeros(s1)),real(V(:,i)),imag(V(:,i)))
				hold on
				
			end
			set(gcf,'color','w');
			set(gca,'box','off');
			xlabel('Re')
			ylabel('Im')
			ax = gca;
			ax.XAxisLocation = 'origin';
			ax.YAxisLocation = 'origin';
			xl = xlim;
			yl = ylim;
			xlim(1.5*xl)
			ylim(1.5*yl)
			for j=1:s1
				text(real(V(j,i))+sign(real(V(j,i)))*(xl(2)-xl(1))*0.1,...
					imag(V(j,i))+sign(real(V(j,i)))*(yl(2)-yl(1))*0.1,x(j))
				hold on
			end
			txt=strcat('output\eigenvector',num2str(i));
			
			savefig(txt);
			
		else
			last = real(W(i));
		end
	end
end

end

