%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                                                                     %
% Programmed by Daniel Du 4/13/2014                                   %
% Status:            Working                                          %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This updater is more advanced and more accurate than the previous 
%microm updater due to the consideration of 3-D distribution.
%This one is so far the most meaningful simplfied model to approximate 
%magnetic force. 

function cnt_new = mmdm_position(m,cnt,H0,th,Ng)
%m is the old dipole moment matrix
%cnt is center position of all particles
%H is external magnetic field strength
%th is the current angle of H
%f is the new dipole moments
%Make sure the particle number Np is larger than 2
%This function has to be used along with micromomentupdateNmer_m.m as this
%one only updates the positions, and the other one updates dipole moments.

%Parameters
a = 1.4e-6;
w = 3;
tmp = size(cnt); 
Np = tmp(1);
h = a/Ng;
cnt_new = zeros(Np,3);
nn = 0;
for ix = 1:(2*Ng+1)
                for iy = 1:(2*Ng+1)
                    for iz = 1:(2*Ng+1)
                        xc = -a+(ix-1)*h;
                        yc = -a+(iy-1)*h;
                        zc = -a+(iz-1)*h;
                        if norm([xc yc zc])<=a
                            nn = nn + 1;
                        end
                    end
                end
end
Nind = nn;

for i = 1:Np
    Hind = zeros(Nind,3);
    for j = 1:Np
        if j ~= i
            nn = 0;
            for ix = 1:2*Ng+1
                for iy = 1:2*Ng+1
                    for iz = 1:2*Ng+1
                        xc = -a+(ix-1)*h;
                        yc = -a+(iy-1)*h;
                        zc = -a+(iz-1)*h;
                        if norm([xc yc zc])<=a
                            nn = nn + 1;
                            rr = (cnt(i,:) + [xc yc zc] - cnt(j,:));
                            cpos(nn,:) = [xc yc zc];
                            Hind(nn,:) = Hind(nn,:) + (1/4/pi)*(3*dot(m(j,:),rr)*rr/(norm(rr))^5-m(j,:)/(norm(rr))^3);                    %induced field
                        end
                    end
                end
            end
        end
    end
    Hind = Hind + [H0*cos(th)*ones(nn,1) H0*sin(th)*ones(nn,1) zeros(nn,1)];
    for ss = 1:nn
        Hc1(ss) = norm(Hind(ss,:));
    end 
    for ss = 1:nn
        Hc2(ss) = norm(Hind(ss,:));
    end   
    for ss = 1:nn
        Hc3(ss) = norm(Hind(ss,:));
    end   
    
    cnt_new(i,1) = cnt(i,1) + dot(Hc1.^w,cpos(:,1))/sum(Hc1.^w);
    cnt_new(i,2) = cnt(i,2) + dot(Hc2.^w,cpos(:,2))/sum(Hc2.^w);
    cnt_new(i,3) = cnt(i,3) + dot(Hc3.^w,cpos(:,3))/sum(Hc3.^w);
end
