function [vw, Vw,state] = wave(lattice, geo, state)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[void, m, void]=size(lattice.VORTEX);

[s1 s2] = size(lattice.COLLOC); 

lemma=ones(1,s1);
var.lemma = lemma;
%gamma_h = gamma(1:s1,:);

var.kappa0= 9.81/((state.STW)^2);


%% validation
var.integral1=30;
var.integral2=30;

var.fitt=0;	%fitt == zero seems to be the way to go, judging from surface plot
%% end validation

%% no wave region
cd = -0.05;


var.xi(:,:,1) = (lattice.VORTEX(:,m/2,1)*lemma)';
var.eta(:,:,1) = (lattice.VORTEX(:,m/2,2)*lemma)';
var.zeta(:,:,1) = (lattice.VORTEX(:,m/2,3)*lemma)';
var.xi(:,:,2) = (lattice.VORTEX(:,m/2+1,1)*lemma)';
var.eta(:,:,2) = (lattice.VORTEX(:,m/2+1,2)*lemma)';
var.zeta(:,:,2) = (lattice.VORTEX(:,m/2+1,3)*lemma)';
var.x = (lattice.COLLOC(:,1)*lemma);
var.y = (lattice.COLLOC(:,2)*lemma);
var.z = (lattice.COLLOC(:,3)*lemma);



% var.z = var.z + mask.*100;
% var.zeta(:,:,1) = var.zeta(:,:,1) + mask.*100;
% var.zeta(:,:,2) = var.zeta(:,:,2) + mask.*100;

mN(:,:,1)=(lattice.N(:,1)*lemma);
mN(:,:,2)=(lattice.N(:,2)*lemma);
mN(:,:,3)=(lattice.N(:,3)*lemma);

% var.message= state.message;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r1(:,:,1) = var.xi(:,:,1)-var.x(:,:);
% r1(:,:,2) = var.eta(:,:,1)-var.y(:,:);
% r1(:,:,3) = var.zeta(:,:,1)-var.z(:,:);
% r2(:,:,1) = var.xi(:,:,2)-var.x(:,:);
% r2(:,:,2) = var.eta(:,:,2)-var.y(:,:);
% r2(:,:,3) = var.zeta(:,:,2)-var.z(:,:);
% 
% R0(:,:,:) = r2(:,:,:)-r1(:,:,:);
% 
% F1=cross(r1,r2,3);
% 
% LF1=sum(F1.^2,3);
% 
% rd=sqrt((LF1./(sum(R0.^2,3))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

near =1e-6;

Vw = quadrature(-pi/2,pi/2,var);

% Vw(:,:,1) = Vw(:,:,1).*(1-(rd<near));
% Vw(:,:,2) = Vw(:,:,2).*(1-(rd<near));
% Vw(:,:,3) = Vw(:,:,3).*(1-(rd<near));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these lines were added to improve the behaviour of close to surface
% collocation points

mask = var.z>-0.1;

mask1 = Vw>0.5;
mask2 = Vw<-0.5;

Vw = Vw.*(1-mask1) + mask1.*0.5;

Vw = Vw.*(1-mask2) + mask2.*-0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vw=sum(Vw.*mN,3);
 
%Vw = zeros(s1,s1,3); vw = zeros(s1,s1);
% state.message=state.message+1;

% uw = integral(:,1);
% vw = integral(:,2);
% ww = integral(:,3);
% 
% 
% Vw = [uw vw ww];

end


function [integral] = quadrature(a, b, var)

[Nu, W] = lgwt(var.integral1,a,b);
nu=Nu';
w=W';

f = ftion(nu,var) ;

[a b c] = size(f);

for j=1:c
	ix(:,:,j) = imag(f(:,:,j)).*cos(nu(j));
	iy(:,:,j) = imag(f(:,:,j)).*sin(nu(j));
	iz(:,:,j) = real(f(:,:,j));
	
end



for j=1:c
Ix(:,:,j)=(ix(:,:,j).*w(j));
Iy(:,:,j)=(iy(:,:,j).*w(j));
Iz(:,:,j)=(iz(:,:,j).*w(j));
end

Ix = -1/(2*pi).*sum(Ix,3);
Iy = -1/(2*pi).*sum(Iy,3);
Iz = -1/(2*pi).*sum(Iz,3);
%is the minus sign correct?


integral(:,:,1) = Ix;
integral(:,:,2) = Iy;
integral(:,:,3) = Iz;
 





end

function [f] = ftion(nu,var)

[s1, s2] =size(nu);

for j=1:s2
	
	pre(:,:,j) = ((var.eta(:,:,2)-var.eta(:,:,1)).*sec(nu(j))+1i.*(var.zeta(:,:,2)-var.zeta(:,:,1)).*tan(nu(j)))./...
		((var.xi(:,:,2)-var.xi(:,:,1)).*cos(nu(j))+(var.eta(:,:,2)-var.eta(:,:,1)).*sin(nu(j))+...
		1i.*(var.zeta(:,:,2)-var.zeta(:,:,1)));
end

if var.fitt
	differ = integrand2(nu,var,2,s2)-integrand2(nu,var,1,s2);
else
	differ = integrand1(nu,var,2,s2)-integrand1(nu,var,1,s2);
end
f = pre.*differ;

% 	for a=10:10:160
% 		for b=10:10:160
% 			figure(4)
% 			plot(squeeze(nu),real(squeeze(f(a,b,:))),squeeze(nu),imag(squeeze(f(a,b,:))));
% 			grid on
% 			pause(0.5)
% 			close all
% 		end
% 	end

end

function [J] = integrand1(nu,var,rl,size)

% if var.message
% 	message1 = 'WAVE: starting... ';
% 	message2 = '(force calculation)';
% 	
% 	f2=waitbar(mod(rl,2)/2,strcat(message1,message2));
% else
% 	f2=waitbar(mod(rl,2)/2,'WAVE: starting... ');
% end

om = omega(nu,var,rl,size);

% t= linspace(0,5,100);

for j=1:size
	
	inte(:,:,j) = integrand4(nu(j),var,rl,om(:,:,j));
	
% 	part11(j) = ((-2*1i).*heaviside(om(:,:,j)).*Kappanu(nu(j),var).*...
% 		exp(Kappanu(nu(j),var).*((var.z+var.zeta(:,:,rl))+1i*om(:,:,j))));
% 	
% 	part12(:,:,j) = (1/pi.*(1./((var.z+var.zeta(:,:,rl))+1i*om(:,:,j)) - ...
% 		Kappanu(nu(j),var)*inte(:,:,j)));
% 	
% 	part13(:,:,j) = 0;
% 	
	J(:,:,j) = ((-2*1i).*heaviside(om(:,:,j)).*Kappanu(nu(j),var).*...
		exp(Kappanu(nu(j),var).*((var.z+var.zeta(:,:,rl))+1i*om(:,:,j)))+...
		1/pi.*(1./((var.z+var.zeta(:,:,rl))+1i*om(:,:,j)) - ...
		Kappanu(nu(j),var)*inte(:,:,j)));
	
	
%	f2=waitbar(mod(rl,2)/2+1/2*j/size,f2);
end

%delete(f2)
	
%% TEST - this code is for visualising the different parts of FS formula which lead to instabilities	
% 	for a=1:10:50
% 		for b=1:10:50
% 			figure(4)
% 			plot(squeeze(nu),real(squeeze(part11(:))),squeeze(nu),imag(squeeze(part11(:))));
% 			grid on
% 			pause(1)
% 			close all
% 			
% 			figure(5)
% 			plot(squeeze(nu),real(squeeze(part12(a,b,:))),squeeze(nu),imag(squeeze(part12(a,b,:))));
% 			grid on
% 			pause(1)
% 			
% 			
% 			figure(6)
% 			plot(squeeze(nu),real(squeeze(part13(a,b,:))),squeeze(nu),imag(squeeze(part13(a,b,:))));
% 			grid on
% 			pause(1)
% 			close all
% 			
% 			
% 			
% % 			figure(3)
% % 			plot(squeeze(nu),real(squeeze(inte(a,b,:))),squeeze(nu),imag(squeeze(inte(a,b,:))));
% % 			grid on
% % 			pause(0.5)
% % 			close all
% % 			%one singularity
% 		end
% 	end
	
	
	

end


function [I] = integrand2giesing(nu,var,rl,om)
% as described in Hess and Smith, 6 degree polynom, highly oscilatory

	m1=0.23721365; m2=0.020654300; m3=-0.00076329700; m4=0.0000097687007;
	n1=-1.49545886; d1=-0.76273617; n2= 0.041806426; d2=- 0.28388363;
	n3=-0.03000591; d3= -0.066786033; n4= 0.0019387339; d4= 0.0112982719;
	n5=-0.00051801555; d5=-0.00087008610; d6=0.00029892040;
	
	B = Kappanu(nu,var).*((var.z+var.zeta(:,:,rl))+1i.*om);
	gamma = 0.5772156649;
	
	[s1, s2, s3] = size(B);
	
	M = -(ones(s1,s2,s3)+m1.*B+m2.*(B.^2)+m3.*(B.^3)+m4.*((var.xi(:,rl)).^4)).*log(B);
	N = -gamma.*(0.99999207.*ones(s1,s2,s3)+n1.*B+n2.*(B.^2)+n3.*((var.z).^3)+n4.*(B.^4)+n5.*(B.^5));
	D = 1+d1.*B  + d2.*(B.^2) + d3.*(B.^4) + d4.*(B.^4) + d5.*(B.^5) + d6.*(B.^6);
	
	I = (M+N)./D;

end


function [I] = integrand3(nu,var,rl,om)
	B = Kappanu(nu,var).*((var.z+var.zeta(:,:,rl))+1i.*om);
	
	[Tu, W] = lgwt(10000,0,100);
	t=Tu';
	w=W';
	
	
	[s1, s2] =size(t);

	for j=1:s2
		II(:,:,j) = exp(-t(j))./(t(j)+B(:,:));
	end
	
	for j=1:s2
		I(:,:,j)=(II(:,:,j).*w(j));
	end
	I = sum(I,3);

end


function [II] = integrand4(nu,var,rl,om)
	B = Kappanu(nu,var).*((var.z+var.zeta(:,:,rl))+1i.*om);
	
	
	
	
	
	[Tu1, W1] = lgwt(round(var.integral2*0.9),0,5);
	[Tu2,W2] = lgwt(round(var.integral2*0.1),5,10000);
	tu1=Tu1';
	w1=W1';
	tu2=Tu2';
	w2=W2';
	
	tu = [tu1 tu2];
	
	w = [w1 w2];
	
	[s1, s2] =size(tu);

	for j=1:s2
		III(:,:,j) = exp(-tu(j))./(tu(j)+B(:,:));
	
		II(:,:,j)=(III(:,:,j).*w(j));
	end
	II = sum(II,3);
	
%% TEST - this code is for visualising the different parts of FS formula which lead to instabilities	

% 	for a=1:50
% 		for b=1:50
% 			figure(4)
% 			plot(squeeze(tu1),real(squeeze(III(a,b,1:80))),squeeze(tu1),imag(squeeze(III(a,b,1:80))));
% 			grid on
% 			pause(0.5)
% 			close all
% 		end
% 	end
	

end




%%%
% Method Alistair and Fitt
%%%



function [J] = integrand2(nu,var,rl,size)

% if var.message
% 	message1 = 'WAVE: starting... ';
% 	message2 = '(force calculation)';
% 	
% 	f2=waitbar(mod(rl,2)/2,strcat(message1,message2));
% else
% 	f2=waitbar(mod(rl,2)/2,'WAVE: starting... ');
% end

om = omega(nu,var,rl,size);

% t= linspace(0,5,100);

for j=1:size
	
	inte2(:,:,j) = integrand5(nu(j),var,rl,om(:,:,j));

	J(:,:,j) = ((-1*1i).*Kappanu(nu(j),var).*...
		exp(Kappanu(nu(j),var).*((var.z+var.zeta(:,:,rl))+1i*om(:,:,j)))+...
		1/pi.*(inte2(:,:,j)));
	
	
	f2=waitbar(mod(rl,2)/2+1/2*j/size,f2);
end

delete(f2)

end

function [II] = integrand5(nu,var,rl,om)
	B = Kappanu(nu,var).*((var.z+var.zeta(:,:,rl))+1i.*om);
	
	
	
	
	
	[Tu1, W1] = lgwt(100,0,5);
	[Tu2,W2] = lgwt(28,5,10000);
	tu1=Tu1';
	w1=W1';
	tu2=Tu2';
	w2=W2';
	
	tu = [tu1 tu2];
	w = [w1 w2];
	
	[s1, s2] =size(tu);

	for j=1:s2
		III(:,:,j) = (tu(j)./(Kappanu(nu,var)-tu(j))).*exp(tu(j).*((var.z+var.zeta(:,:,rl))+1i.*om));
	end
	
	for j=1:s2
		II(:,:,j)=(III(:,:,j).*w(j));
	end
	II = sum(II,3);
	
%% TEST - this code is for visualising the different parts of FS formula which lead to instabilities	

% 	for a=1:50
% 		for b=1:50
% 			figure(4)
% 			plot(squeeze(tu1),real(squeeze(III(a,b,1:80))),squeeze(tu1),imag(squeeze(III(a,b,1:80))));
% 			grid on
% 			pause(0.5)
% 			close all
% 		end
% 	end
	

end



function [om] = omega(nu,var,rl,size)


for j=1:size
	om(:,:,j) = (var.x-var.xi(:,:,rl)).*cos(nu(j))+(var.y-var.eta(:,:,rl)).*sin(nu(j));
end

end

function [kappanu] = Kappanu(nu,var)

	kappanu = var.kappa0.*(sec(nu)).^2;

end

function y = heaviside(x)
y = (x > 0) + 0.5*(x == 0);
end
