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

function m_new = mmdm_moment(m,cnt,H0,th,Chi)
%m is the old dipole moment matrix
%cnt is center position of all particles
%H is external magnetic field strength
%th is the current angle of H
%f is the new dipole moments
%Make sure the particle number Np is larger than 2
%This function has to be used along with micromomentupdateNmer_cnt.m as 
%this one only updates dipole moments, and the other one updates positions.

%Parameters
a = 1.4e-6;
tmp = size(cnt); 
Np = tmp(1);
m_new = m;
Hext = [H0*cos(th) H0*sin(th) 0];   %External field
for i = 1:Np
    %i is the reference particle
    %j will be the current particle
        Hind = [0 0 0];
        for j = 1:Np
            if j~=i
                r = (cnt(i,:)- cnt(j,:));
                Hind = Hind + (1/4/pi)*(3*dot(m(j,:),r)*r/(norm(r))^5-m(j,:)/(norm(r))^3);                    %induced field
            end
        end
        m_new(i,:) = 4/3 * pi * a^3 * Chi * (Hext+Hind);
end
