%density-dependent dipole moment updater
function f = mdm(m,cnt,H,th,Kai)
%m is the old dipole moment matrix
%cnt is center position of all particles
%H is external magnetic field strength
%th is the current angle of H
%Fac is the excluded region factor, typically 5~10
%f is the new 
%margin is the cropped distance from the edge in pixel
%Basically H and th should be external variables

%Parameters
a = 1.4e-6;
H0 = [H*cos(th) H*sin(th)];   %External field
D = 2*a;
f = m;
for i = 1:length(cnt)
    %i is the reference particle
    %j will be the current particle
        Hind = [0 0];
        for j = 1:length(cnt)
            r = [cnt(i,1)-cnt(j,1) cnt(i,2)-cnt(j,2)];
            if norm(r)>=0.5*D && norm(r) <= 10*D
                Hind = Hind + (1/4/pi)*(3*dot(m(j,:),r)*r/(norm(r))^5-m(j,:)/(norm(r))^3);                    %induced field
            end
        end
        f(i,:) = 4/3 * pi * a^3 * Kai * (H0+Hind);
end
