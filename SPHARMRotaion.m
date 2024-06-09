function YR = SPHARMRotaion(Y,Rt,SPHARM_Order,R2C1,C2R1,UY90p,UY90n)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

%R2C = Real2complex(SPHARM_Order);
%C2R = Complex2real(SPHARM_Order);
%R2C1=transpose(R2C);
%C2R1=transpose(C2R);
%thetay=pi/2;
%UY90p= SpYrotation (SPHARM_Order,thetay);
%thetay=-pi/2;
%UY90n= SpYrotation (SPHARM_Order,thetay);

    [alfa,beta,gamma] = DecomrotationPro(Rt);
    Ylr= Y;
    Ylc= C2R1*Ylr;
    Ualfa = SpZrotation(SPHARM_Order,alfa+pi/2);
    Ylc1=Ualfa*Ylc;
    UZbeta=SpZrotation(SPHARM_Order,beta);
    Ylc2=UY90n*UZbeta*UY90p*Ylc1;
    Ugamma = SpZrotation(SPHARM_Order,gamma-pi/2);
    Ylc3=Ugamma*Ylc2;
    Ylr3=R2C1*Ylc3;
    Ylr3= real(Ylr3);
    YR=Ylr3;

end