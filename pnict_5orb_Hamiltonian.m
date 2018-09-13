function [f] = pnict_5orb_Hamiltonian(n,k,mu,xs,sw,eta,soc,bc,g)
%sw=zeros(1,n);
% Tight-binding Hamiltonian (5*2*4*n by 5*2*4*n)
% n sites, 5 orbitals: xz, yz, x^2-y^2, xy, 3z^2-r^2
A = zeros(40*n,40*n); % init sparse matrix

z5=zeros(5);
% "unit" matrices
[M11,M22,M33,M44,M55,M12,M13,M14,M15,M23,M24,M25,M34,M35,M45]=deal(z5,z5,z5,z5,z5,z5,z5,z5,z5,z5,z5,z5,z5,z5,z5);
M11(1,1) = 1/2; % diags are 1/2 for adding with h.c.
M22(2,2) = 1/2;
M33(3,3) = 1/2;
M44(4,4) = 1/2;
M55(5,5) = 1/2;
M12(1,2) = 1;
M13(1,3) = 1;
M23(2,3) = 1;
M14(1,4) = 1;
M24(2,4) = 1;
M15(1,5) = 1;
M25(2,5) = 1;
M34(3,4) = 1;
M35(3,5) = 1;
M45(4,5) = 1;

O=eye(2);
O5 = eye(5);
%Pauli Matrices
X=[0 1; 1 0];
Y=[0 -1i; 1i 0];
Z=[1 0; 0 -1];
OP=(O+Z)/2;
OM=(O-Z)/2;


%%%

% On-site energies
E1 = .13;
E2 = .13;
E3 = -.22;
E4 = .3;
E5 = -.211;

% Intraorbital hopping parameters
t1x = -.14;
t1y = -.4;
t1xy = .28;
t1xx = .02;
t1xxy = -.035;
t1xyy = .005;
t1xxyy = .035;

t3x = .35;
t3xy = -.105;
t3xx = -.02;

t4x = .23;
t4xy = .15;
t4xx = -.03;
t4xxy = -.03;
t4xxyy = -.03;

t5x = -.1;
t5xx = -.04;
t5xxy = .02;
t5xxyy = -.01;

% Interorbital hopping parameters

t12xy = .05;
t12xxy = -.015;
t12xxyy = .035;

t13x = -.354;
t13xy = .099;
t13xxy = .021;

t14x = .339;
t14xy = .014;
t14xxy = .028;

t15x = -.198;
t15xy = -.085;
t15xxyy = -.014;

t34xxy = -.01;

t35x = -.3;
t35xxy = -.02;

t45xy = -.15;
t45xxyy = .01;

%%%

% Intraorbital kinetic energy terms

E1os = E1 + 2*t1y*cos(k) - 2*t1xx*cos(2*k) - mu;
E1c1 = 2*t1x + 4*t1xy*cos(k) + 4*t1xyy*cos(2*k);
E1c2 = 2*t1xx + 4*t1xxy*cos(k) + 4*t1xxyy*cos(2*k);
E1s1 = 0;
E1s2 = 0;

E2os = E2 + 2*t1x*cos(k) + 2*t1xx*cos(2*k) - mu;
E2c1 = 2*t1y + 4*t1xy*cos(k) + 4*t1xxy*cos(2*k);
E2c2 = -2*t1xx + 4*t1xyy*cos(k) + 4*t1xxyy*cos(2*k);
E2s1 = 0;
E2s2 = 0;

E3os = E3 + 2*t3x*cos(k) + 2*t3xx*cos(2*k) - mu;
E3c1 = 2*t3x + 4*t3xy*cos(k);
E3c2 = 2*t3xx;
E3s1 = 0;
E3s2 = 0;

E4os = E4 + 2*t4x*cos(k) + 2*t4xx*cos(2*k) - mu;
E4c1 = 2*t4x + 4*t4xy*cos(k) + 4*t4xxy*cos(2*k);
E4c2 = 2*t4xx + 4*t4xxy*cos(k) + 4*t4xxyy*cos(2*k);
E4s1 = 0;
E4s2 = 0;

E5os = E5 + 2*t5x*cos(k) + 2*t5xx*cos(2*k) - mu;
E5c1 = 2*t5x + 4*t5xxy*cos(2*k);
E5c2 = 2*t5xx + 4*t5xxy*cos(k) + 4*t5xxyy*cos(2*k);
E5s1 = 0;
E5s2 = 0;

% Interorbital kinetic energy terms

E12os = 0;
E12c1 = 0;
E12c2 = 0;
E12s1 = -4*t12xy*sin(k) - 4*t12xxy*sin(2*k);
E12s2 = -4*t12xxy*sin(k) - 4*t12xxyy*sin(2*k);

E13os = 2*1i*t13x*sin(k);
E13c1 = 4*1i*t13xy*sin(k) - 4*1i*t13xxy*sin(2*k);
E13c2 = 4*1i*t13xxy*sin(k);
E13s1 = 0;
E13s2 = 0;

E23os = 0;
E23c1 = 0;
E23c2 = 0;
E23s1 = -2*1i*t13x - 4*1i*t13xy*cos(k) - 4*1i*t13xxy*cos(2*k);
E23s2 = 4*1i*t13xxy*cos(k);

E14os = 0;
E14c1 = 0;
E14c2 = 0;
E14s1 = 2*1i*t14x + 4*1i*t14xy*cos(k);
E14s2 = 4*1i*t14xxy*cos(k);

E24os = 2*1i*t14x*sin(k);
E24c1 = 4*1i*t14xy*sin(k) + 4*1i*t14xxy*sin(2*k);
E24c2 = 0;
E24s1 = 0;
E24s2 = 0;

E15os = 2*1i*t15x*sin(k);
E15c1 = -4*1i*t15xy*sin(k);
E15c2 = -4*1i*t15xxyy*sin(2*k);
E15s1 = 0;
E15s2 = 0;

E25os = 0;
E25c1 = 0;
E25c2 = 0;
E25s1 = 2*1i*t15x - 4*1i*t15xy*cos(k);
E25s2 = -4*1i*t15xxyy*cos(2*k);

E34os = 0;
E34c1 = 0;
E34c2 = 0;
E34s1 = 4*t34xxy*sin(2*k);
E34s2 = -4*t34xxy*sin(k);

E35os = -2*t35x*cos(k);
E35c1 = 2*t35x - 4*t35xxy*cos(2*k);
E35c2 = 4*t35xxy*cos(k);
E35s1 = 0;
E35s2 = 0;

E45os = 0;
E45c1 = 0;
E45c2 = 0;
E45s1 = 4*t45xy*sin(k);
E45s2 = 4*t45xxyy*sin(2*k);

%%%

Aos_ur = E1os*M11 + E2os*M22 + E3os*M33 + E4os*M44 + E5os*M55 +...
    E13os*M13 + E24os*M24 + E15os*M15 + E35os*M35;

Ac1_ur = E1c1*M11 + E2c1*M22 + E3c1*M33 + E4c1*M44 + E5c1*M55 +...
    E12c1*M12 + E13c1*M13 + E14c1*M14 + E15c1*M15 +...
    E23c1*M23 + E24c1*M24 + E25c1*M25 +...
    E34c1*M34 + E35c1*M35 +...
    E45c1*M45;

Ac2_ur = E1c2*M11 + E2c2*M22 + E3c2*M33 + E4c2*M44 + E5c2*M55 +...
    E12c2*M12 + E13c2*M13 + E14c2*M14 + E15c2*M15 +...
    E23c2*M23 + E24c2*M24 + E25c2*M25 +...
    E34c2*M34 + E35c2*M35 +...
    E45c2*M45;

As1_ur = E1s1*M11 + E2s1*M22 + E3s1*M33 + E4s1*M44 + E5s1*M55 +...
    E12s1*M12 + E13s1*M13 + E14s1*M14 + E15s1*M15 +...
    E23s1*M23 + E24s1*M24 + E25s1*M25 +...
    E34s1*M34 + E35s1*M35 +...
    E45s1*M45;

As2_ur = E1s2*M11 + E2s2*M22 + E3s2*M33 + E4s2*M44 + E5s2*M55 +...
    E12s2*M12 + E13s2*M13 + E14s2*M14 + E15s2*M15 +...
    E23s2*M23 + E24s2*M24 + E25s2*M25 +...
    E34s2*M34 + E35s2*M35 +...
    E45s2*M45;

Aosk = Aos_ur + Aos_ur';
Ac1 = Ac1_ur + Ac1_ur';
Ac2 = Ac2_ur + Ac2_ur';
As1 = As1_ur + As1_ur';
As2 = As2_ur + As2_ur';


% Left hops
AL1k=.5*(Ac1 - 1i*As1);
AL2k=.5*(Ac2 - 1i*As2);

% Right hops
AR1k=.5*(Ac1 + 1i*As1);
AR2k=.5*(Ac2 + 1i*As2);


%%%%%%%%%%%%%%%%%%

% KE terms for ky+Q

% On-site energies
EQ1 = .13;
EQ2 = .13;
EQ3 = -.22;
EQ4 = .3;
EQ5 = -.211;

% Intraorbital kinetic energy terms

EQ1os = EQ1 + 2*t1y*cos(k+pi) - 2*t1xx*cos(2*k) - mu;
EQ1c1 = 2*t1x + 4*t1xy*cos(k+pi) + 4*t1xyy*cos(2*k);
EQ1c2 = 2*t1xx + 4*t1xxy*cos(k+pi) + 4*t1xxyy*cos(2*k);
EQ1s1 = 0;
EQ1s2 = 0;

EQ2os = EQ2 + 2*t1x*cos(k+pi) + 2*t1xx*cos(2*k) - mu;
EQ2c1 = 2*t1y + 4*t1xy*cos(k+pi) + 4*t1xxy*cos(2*k);
EQ2c2 = -2*t1xx + 4*t1xyy*cos(k+pi) + 4*t1xxyy*cos(2*k);
EQ2s1 = 0;
EQ2s2 = 0;

EQ3os = EQ3 + 2*t3x*cos(k+pi) + 2*t3xx*cos(2*k) - mu;
EQ3c1 = 2*t3x + 4*t3xy*cos(k+pi);
EQ3c2 = 2*t3xx;
EQ3s1 = 0;
EQ3s2 = 0;

EQ4os = EQ4 + 2*t4x*cos(k+pi) + 2*t4xx*cos(2*k) - mu;
EQ4c1 = 2*t4x + 4*t4xy*cos(k+pi) + 4*t4xxy*cos(2*k);
EQ4c2 = 2*t4xx + 4*t4xxy*cos(k+pi) + 4*t4xxyy*cos(2*k);
EQ4s1 = 0;
EQ4s2 = 0;

EQ5os = EQ5 + 2*t5x*cos(k+pi) + 2*t5xx*cos(2*k) - mu;
EQ5c1 = 2*t5x + 4*t5xxy*cos(2*k);
EQ5c2 = 2*t5xx + 4*t5xxy*cos(k+pi) + 4*t5xxyy*cos(2*k);
EQ5s1 = 0;
EQ5s2 = 0;

% Interorbital kinetic energy terms

EQ12os = 0;
EQ12c1 = 0;
EQ12c2 = 0;
EQ12s1 = -4*t12xy*sin(k+pi) - 4*t12xxy*sin(2*k);
EQ12s2 = -4*t12xxy*sin(k+pi) - 4*t12xxyy*sin(2*k);

EQ13os = 2*1i*t13x*sin(k+pi);
EQ13c1 = 4*1i*t13xy*sin(k+pi) - 4*1i*t13xxy*sin(2*k);
EQ13c2 = 4*1i*t13xxy*sin(k+pi);
EQ13s1 = 0;
EQ13s2 = 0;

EQ23os = 0;
EQ23c1 = 0;
EQ23c2 = 0;
EQ23s1 = -2*1i*t13x - 4*1i*t13xy*cos(k+pi) - 4*1i*t13xxy*cos(2*k);
EQ23s2 = 4*1i*t13xxy*cos(k+pi);

EQ14os = 0;
EQ14c1 = 0;
EQ14c2 = 0;
EQ14s1 = 2*1i*t14x + 4*1i*t14xy*cos(k+pi);
EQ14s2 = 4*1i*t14xxy*cos(k+pi);

EQ24os = 2*1i*t14x*sin(k+pi);
EQ24c1 = 4*1i*t14xy*sin(k+pi) + 4*1i*t14xxy*sin(2*k);
EQ24c2 = 0;
EQ24s1 = 0;
EQ24s2 = 0;

EQ15os = 2*1i*t15x*sin(k+pi);
EQ15c1 = -4*1i*t15xy*sin(k+pi);
EQ15c2 = -4*1i*t15xxyy*sin(2*k);
EQ15s1 = 0;
EQ15s2 = 0;

EQ25os = 0;
EQ25c1 = 0;
EQ25c2 = 0;
EQ25s1 = 2*1i*t15x - 4*1i*t15xy*cos(k+pi);
EQ25s2 = -4*1i*t15xxyy*cos(2*k);

EQ34os = 0;
EQ34c1 = 0;
EQ34c2 = 0;
EQ34s1 = 4*t34xxy*sin(2*k);
EQ34s2 = -4*t34xxy*sin(k+pi);

EQ35os = -2*t35x*cos(k+pi);
EQ35c1 = 2*t35x - 4*t35xxy*cos(2*k);
EQ35c2 = 4*t35xxy*cos(k+pi);
EQ35s1 = 0;
EQ35s2 = 0;

EQ45os = 0;
EQ45c1 = 0;
EQ45c2 = 0;
EQ45s1 = 4*t45xy*sin(k+pi);
EQ45s2 = 4*t45xxyy*sin(2*k);

%%%

AosQ_ur = EQ1os*M11 + EQ2os*M22 + EQ3os*M33 + EQ4os*M44 + EQ5os*M55 +...
    EQ13os*M13 + EQ24os*M24 + EQ15os*M15 + EQ35os*M35;

Ac1Q_ur = EQ1c1*M11 + EQ2c1*M22 + EQ3c1*M33 + EQ4c1*M44 + EQ5c1*M55 +...
    EQ12c1*M12 + EQ13c1*M13 + EQ14c1*M14 + EQ15c1*M15 +...
    EQ23c1*M23 + EQ24c1*M24 + EQ25c1*M25 +...
    EQ34c1*M34 + EQ35c1*M35 +...
    EQ45c1*M45;

Ac2Q_ur = EQ1c2*M11 + EQ2c2*M22 + EQ3c2*M33 + EQ4c2*M44 + EQ5c2*M55 +...
    EQ12c2*M12 + EQ13c2*M13 + EQ14c2*M14 + EQ15c2*M15 +...
    EQ23c2*M23 + EQ24c2*M24 + EQ25c2*M25 +...
    EQ34c2*M34 + EQ35c2*M35 +...
    EQ45c2*M45;

As1Q_ur = EQ1s1*M11 + EQ2s1*M22 + EQ3s1*M33 + EQ4s1*M44 + EQ5s1*M55 +...
    EQ12s1*M12 + EQ13s1*M13 + EQ14s1*M14 + EQ15s1*M15 +...
    EQ23s1*M23 + EQ24s1*M24 + EQ25s1*M25 +...
    EQ34s1*M34 + EQ35s1*M35 +...
    EQ45s1*M45;

As2Q_ur = EQ1s2*M11 + EQ2s2*M22 + EQ3s2*M33 + EQ4s2*M44 + EQ5s2*M55 +...
    EQ12s2*M12 + EQ13s2*M13 + EQ14s2*M14 + EQ15s2*M15 +...
    EQ23s2*M23 + EQ24s2*M24 + EQ25s2*M25 +...
    EQ34s2*M34 + EQ35s2*M35 +...
    EQ45s2*M45;

AosQ = AosQ_ur + AosQ_ur';
Ac1Q = Ac1Q_ur + Ac1Q_ur';
Ac2Q = Ac2Q_ur + Ac2Q_ur';
As1Q = As1Q_ur + As1Q_ur';
As2Q = As2Q_ur + As2Q_ur';


% Left hops
AL1Q=.5*(Ac1Q - 1i*As1Q);
AL2Q=.5*(Ac2Q - 1i*As2Q);

% Right hops
AR1Q=.5*(Ac1Q + 1i*As1Q);
AR2Q=.5*(Ac2Q + 1i*As2Q);

%%%%%%

Aos=[Aosk,zeros(5);zeros(5),AosQ];
AL1=[AL1k,zeros(5);zeros(5),AL1Q];
AL2=[AL2k,zeros(5);zeros(5),AL2Q];
AR1=[AR1k,zeros(5);zeros(5),AR1Q];
AR2=[AR2k,zeros(5);zeros(5),AR2Q];

%%%%%%%%%%%%%%%%%%%

% sc pairing
%
B=zeros(n,40,40);
Bs=zeros(n,40,40);
for j=1:n
%B(j,:,:)=xs(j)*kron(kron(X,O),kron(Z,O5));
B(j,:,:)=[zeros(20,20),xs(j)*kron(O,kron(Z,O5));
  conj(xs(j))*kron(O,kron(Z,O5)), zeros(20,20)];
Bs(j,:,:)=[zeros(20,20),sw(j)*kron(O,kron(Z,O5));
  conj(sw(j))*kron(O,kron(Z,O5)), zeros(20,20)];
end
%B(n,:,:)=xs(n)*kron(kron(X,O),kron(Z,O5));
HR1=kron(kron(Z,O),AR1);
HL1=kron(kron(Z,O),AL1);
Hos=kron(kron(Z,O),Aos);
HR2=kron(kron(Z,O),AR2);
HL2=kron(kron(Z,O),AL2);
%

%
% Antiferromagnetic SDW
C=eta*kron(kron(O,Z),eye(10)); % For even-sites. -C for odd-sites.
%

% Spin-orbit coupling
Atuu=[0 -1i 0 0 0; 1i 0 0 0 0; 0 0 0 -2*1i 0; 0 0 2*1i 0 0; 0 0 0 0 0]; 
Atdd=-Atuu;
Atdu=[0 0 1 1i -sqrt(3);...
      0 0 -1i 1 -1i*sqrt(3);...
      -1 1i 0 0 0;...
      -1i -1 0 0 0;...
      sqrt(3) 1i*sqrt(3) 0 0 0]; 
Atud=Atdu';

%Hsoc_d = soc*[kron(O,Atuu), kron(O,Atud);kron(O,Atdu),kron(O,Atdd)];
%Hsoc_od = soc*[kron(X,Atuu), kron(X,Atud);kron(X,Atdu),kron(X,Atdd)];

%
Hsoc_d_p = soc*[kron(O,Atuu), kron(O,Atud);kron(O,Atdu),kron(O,Atdd)];
Hsoc_od_p = soc*[kron(X,Atuu), kron(X,Atud);kron(X,Atdu),kron(X,Atdd)];
%{
Hsoc_d_p = soc*[kron(O,Atuu), kron(O,Atdu);kron(O,Atud),kron(O,Atdd)];
Hsoc_od_p = soc*[kron(X,Atuu), kron(X,Atdu);kron(X,Atud),kron(X,Atdd)];
%}

%Hsoc_d = [Hsoc_d_p,zeros(20);zeros(20),-conj(Hsoc_d_p)];
%Hsoc_od = [Hsoc_od_p,zeros(20);zeros(20),-conj(Hsoc_od_p)];

Hsoc_d = [Hsoc_d_p,zeros(20);zeros(20),-Hsoc_d_p];%-Hsoc_d_p
Hsoc_od = [Hsoc_od_p,zeros(20);zeros(20),-Hsoc_od_p];%-Hsoc_od_p

%Uxsbasis=kron([1,0,0,0;0,1,0,0;0,0,0,1;0,0,-1,0],eye(10));
%Hsoc_d = Uxsbasis*kron(Z,Hsoc_d_p)*Uxsbasis';
%Hsoc_od = Uxsbasis*kron(Z,Hsoc_od_p)*Uxsbasis';
%
% diagonal terms
for j=1:80:(40*n -79)
    Bsi(:,:)=g*Bs(1+(j-1)/40,:,:);
 %  A(j:(j+39),j:(j+39)) = Hos - C + kron(Z,Hsoc_d - Hsoc_od);
    A(j:(j+39),j:(j+39)) = Hos + Bsi - C + Hsoc_d - Hsoc_od;
 %  A(j:(j+39),j:(j+39)) = Hos - C + [Hsoc_d - Hsoc_od,zeros(20);zeros(20),-conj(Hsoc_d - Hsoc_od)];
end
for j=41:80:(40*n -39)
     Bsi(:,:)=g*Bs(1+(j-1)/40,:,:);
%    A(j:(j+39),j:(j+39)) = Hos + C + kron(Z,Hsoc_d + Hsoc_od);
A(j:(j+39),j:(j+39)) = Hos + Bsi + C + Hsoc_d + Hsoc_od;
%A(j:(j+39),j:(j+39)) = Hos + C + [Hsoc_d + Hsoc_od,zeros(20);zeros(20),-conj(Hsoc_d + Hsoc_od)];
end

% 
% right hopping terms (lower left triangle)
for j=41:40:(40*n -39)
    Bi(:,:)=B((j-1)/40,:,:);
    A(j:(j+39),(j-40):(j-1)) = HR1+g*Bi';
end
for j=81:40:(40*n -39) 
    A(j:(j+39),(j-80):(j-41)) = HR2;
end

% left hopping terms (upper right triangle)
for j=41:40:(40*n -39)
    Bi(:,:)=B((j-1)/40,:,:);
    A((j-40):(j-1),j:(j+39)) = HL1+g*Bi;
end
for j=81:40:(40*n -39)
    A((j-80):(j-41),j:(j+39)) = HL2;
end

if bc==1
% periodic bc
Bn(:,:)=B(n,:,:);
A((40*n -39):(40*n),1:40)=HL1+g*Bn;
A(1:40,(40*n -39):(40*n))=HR1+g*Bn';

A((40*n -79):(40*n -40),1:40)=HL2;
A(1:40,(40*n -79):(40*n -40))=HR2;

A((40*n -39):(40*n),41:80)=HL2;
A(41:80,(40*n -39):(40*n))=HR2;

elseif bc==0
    % open boundaries
end

f=A;
end