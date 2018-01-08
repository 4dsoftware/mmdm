%Fchain Run
clear
clc
figure(1);

B0 = 6;
Np = 5;
fac = 1.1;
a = 1.4e-6;
D = 2*a;
Nup = 4; %Total update number
q = fac*D;
th = (0:1:180)/180*pi;
Kai=1.6*0.6*3/3.96;  %M-270 based on instruction and Dynabeads word file
mu_0 = pi * 4e-7; %vacuum permeability
H = B0/1e4/mu_0; % in Gauss
Fmag1 = zeros(1,length(th));
Fmag2 = zeros(1,length(th));

Fmagmic1 = zeros(1,length(th));
Fmagmic2 = zeros(1,length(th));

for k = 1:length(q)
    rr = [1 0];
    cntr = [q(k)*(1:Np)' zeros(Np,1)]; %chain in x direction
    for j = 1:length(th) 
        for i = 1:Np
            m(i,:) = 4/3 * pi * a^3 * Kai * H * [cos(th(j)) sin(th(j))];  %dipole moment
            mmic(i,:) = 4/3 * pi * a^3 * Kai * H * [cos(th(j)) sin(th(j)) 0];  %micro dipole moment
        end
               
        iup = 1; %update number
        cntdu0 = [cntr zeros(Np,1)]; %initial centers
        cntdu = [cntr zeros(Np,1)]; %centers
        while iup <= Nup    
            m = mdm(m,cntr,H,th(j),Kai); %MDM
            mmic = mmdm_moment(mmic,cntdu,H,th(j),Kai); %MMDM, update moment
            cntdu = mmdm_position(mmic,cntdu0,H,th(j),6); %MMDM, update positions (centers)
            iup = iup+1;
        end
        for pn = 2:Np
                rrm = cntr(pn,:) - cntr(1,:);
                rrdu = cntdu(pn,:) - cntdu(1,:);
                F1 = (3*mu_0/(4*pi*(norm(rrm))^4))*(dot(m(1,:),rrm/norm(rrm))*m(pn,:) + dot(m(pn,:),rrm/norm(rrm))*m(1,:) + dot(m(1,:),m(pn,:))*rrm/norm(rrm) -5*dot(m(1,:),rrm/norm(rrm))*dot(m(pn,:),rrm/norm(rrm))*rrm/norm(rrm));
                Fmag1(j) = Fmag1(j) + F1(:,1); %sum forces from all other particles
                F2 = (3*mu_0/(4*pi*(norm(rrdu))^4))*(dot(mmic(1,:),(rrdu/norm(rrdu)))*mmic(pn,:) + dot(mmic(pn,:),(rrdu/norm(rrdu)))*mmic(1,:) + dot(mmic(1,:),mmic(pn,:))*(rrdu/norm(rrdu)) -5*dot(mmic(1,:),(rrdu/norm(rrdu)))*dot(mmic(pn,:),(rrdu/norm(rrdu)))*((rrdu/norm(rrdu))));
                Fmagmic1(j) = Fmagmic1(j) + F2(:,1); %sum forces from all other particles
        end
    end
end

plot(th/pi*180,Fmag1,'k:','Linewidth',2.5); %MDM
hold on
plot(th/pi*180,Fmagmic1,'b-','linewidth',2.5); %MMDM
hold off
legend({'MDM','MMDM'});
ylabel('F_1 (N)');
xlabel('\theta (deg)');
title('5-particle Chain');
xlim([0 180]);


