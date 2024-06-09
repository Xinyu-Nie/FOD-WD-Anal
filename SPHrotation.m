function Sph = SPHrotation(v0,peakv,Sph0,SPHARM_Order,R2C1,C2R1,UY90p,UY90n)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
%  Rotate a FOD function to align its maximum peak to v0 direction
       v01=v0';
       v02=peakv';
       v01=vectornomrlization(v01);
       v02=vectornomrlization(v02);
       nv=cross(v02,v01);
       if norm(nv)>0.001
       nv=vectornomrlization(nv);
       alfa=acos(sum(v01.*v02));
       v21=rotationMatrix(nv,alfa)*v02;
       v22=rotationMatrix(nv,-alfa)*v02;
       inv1=abs(sum(v21.*v01));
       inv2=abs(sum(v22.*v01));
       if inv1 >inv2
           RM=rotationMatrix(nv,alfa);
       else
           RM=rotationMatrix(nv,-alfa);
       end
       Sph=SPHARMRotaion(Sph0,RM',SPHARM_Order,R2C1,C2R1,UY90p,UY90n);
       else
           Sph=Sph0;
       end


end