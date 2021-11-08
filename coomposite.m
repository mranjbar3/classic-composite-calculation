clear all;
clc;

Inpot=xlsread('source.xlsx',1);%Input Data from Excel document
%Input Data Valiable
E_x=Inpot(:,1);    E_y=Inpot(:,2);   E_s=Inpot(:,3);   No_xy=Inpot(:,4); theta=Inpot(:,5);
k=size(E_x);    ti=Inpot(:,6);  X_T=Inpot(:,7); X_C=Inpot(:,8); Y_T=Inpot(:,9);
Y_C=Inpot(:,10);    S=Inpot(:,11);

clear Inpot

k=k(1,1);

Inpot=xlsread('source.xlsx',2);
Si_on(:,:)=Inpot(1,1:3)';  Ep_on(:,:)=Inpot(2,1:3)';  Si_off(:,:)=Inpot(3,1:3)';  Ep_off(:,:)=Inpot(4,1:3)';
N(:,:)=Inpot(5,1:3)'; M(:,:)=Inpot(6,1:3)';

clear Inpot

No_yx=E_y/E_x*No_xy;%define minor p.r.

A=zeros(3);B=zeros(3);D=zeros(3);S_on=zeros(3,3,k);C=zeros(3,3,k);S_off=zeros(3,3,k);
Q=zeros(3,3,k);T_sig=zeros(3,3,k);T_eps=zeros(3,3,k);
E_1=zeros(k,1);E_2=zeros(k,1);E_6=zeros(k,1);No_12=zeros(k,1);No_21=zeros(k,1);
SC_61=zeros(k,1);SC_62=zeros(k,1);NC_16=zeros(k,1);NC_26=zeros(k,1);
e_wu=zeros(1,k);    e_fp=zeros(1,k);    e_fmp=zeros(1,k);   e_fm=zeros(1,k);
e_mp=zeros(1,k);    e_mm=zeros(1,k);

%find the laminate thickness
h=0;
for i=1:k
    h=h+ti(i);
end

%define C & S & Q & A & D matrices
h2=h/2;
for i=1:k
    h1=h2-ti(i);
    
    m=cosd(theta(i));
    n=sind(theta(i));
    
    S_on(:,:,i)=[1/E_x(i) -No_yx(i)/E_y(i) 0;-No_xy(i)/E_x(i) 1/E_y(i) 0;0 0 1/E_s(i)];%on-axis compliance matrix
    C(:,:,i)=inv(S_on(:,:,i));%on-axis stiffness matrix
    T_sig(:,:,i)=[m.^2 n.^2 2*m*n;n.^2 m.^2 -2*m*n;-m*n m*n m.^2-n.^2];%Transformation matrix of Stress
    T_eps(:,:,i)=[m.^2 n.^2 m*n;n.^2 m.^2 -m*n;-2*m*n 2*m*n m.^2-n.^2];%Transformation matrix of Strain
    Q(:,:,i)=T_sig(:,:,i)\C(:,:,i)*T_eps(:,:,i);%off-axis stiffness matrix
    S_off(:,:,i)=inv(Q(:,:,i));%off-axis compliance matrix
    B=B+Q(:,:,i).*(h2^2-h1^2)/2;
    A=A+Q(:,:,i).*ti(i);
    D=D+Q(:,:,i).*(h2^3-h1^3)/3;
    
    h2=h2-ti(i);
end
T=[A B;B D];
clear m n h1 h2 B D

% Mechanical Specification of symmetric Laminate
a=inv(A);
E_1o=1/(h*a(1,1));
E_2o=1/(h*a(2,2));
E_6o=1/(h*a(3,3));
No_12o=-a(1,2)/a(2,2);
No_21o=-a(2,1)/a(1,1);
NC_16o=a(1,3)/a(3,3);
NC_26o=a(2,3)/a(3,3);
SC_61o=a(3,1)/a(1,1);
SC_62o=a(3,2)/a(2,2);

% Mechanical Specification of Unidirectional Ply
for i=1:k
    E_1(i)=1/S_off(1,1,i);%Normal stiffness
    E_2(i)=1/S_off(2,2,i);%Normal stiffness
    E_6(i)=1/S_off(3,3,i);%Shear stiffness
    No_21(i)=-S_off(2,1,i)/S_off(1,1,i);%Longitudinal poisson's ratio
    No_12(i)=-S_off(1,2,i)/S_off(2,2,i);%Transverse poisson's ratio
    SC_61(i)=S_off(3,1,i)/S_off(1,1,i);%Shear Coupling Coefficent
    SC_62(i)=S_off(3,2,i)/S_off(2,2,i);%Shear Coupling Coefficent
    NC_16(i)=S_off(1,3,i)/S_off(3,3,i);%Normal Coupling Coefficent
    NC_26(i)=S_off(2,3,i)/S_off(3,3,i);%Normal Coupling Coefficent
end
n1=0;
% Stress and Strain of Laminate and U.D.ply
if abs(N(1))>0 || abs(N(2))>0 || abs(N(3))>0 || abs(M(1))>0 || abs(M(2))>0 || abs(M(3))>0 %find strain
    P=[N;M];
    St=T\P;
    n1=1;
end

if n1==1 %find Total strain of each ply
    h2=h/2;
    for i=1:k/2
        Ep_off(1,i)=St(1)+St(4)*h2;
        Ep_off(2,i)=St(2)+St(5)*h2;
        Ep_off(3,i)=St(3)+St(6)*h2;
        Ep_off(1,k-i+1)=St(1)-St(4)*h2;
        Ep_off(2,k-i+1)=St(2)-St(5)*h2;
        Ep_off(3,k-i+1)=St(3)-St(6)*h2;
        h2=h2-ti(i);
    end
else
    St=Ep_off;
    Ep_off(1,1:k)=St(1);
    Ep_off(2,1:k)=St(2);
    Ep_off(3,1:k)=St(3);
end

for i=1:k
    if abs(Ep_off(1))>0 || abs(Ep_off(2))>0 || abs(Ep_off(3))>0
        Si_off(:,i)=Q(:,:,i)*Ep_off(:,i);
        Ep_on(:,i)=T_eps(:,:,i)*Ep_off(:,i);
        Si_on(:,i)=C(:,:,i)*Ep_on(:,i);
        e_wu(i)=1/(X_T(i)*X_C(i))*Si_on(1,i).^2+1/(X_T(i)*X_C(i))*Si_on(2,i).^2+1/S(i).^2*Si_on(3,i).^2-0.5*(X_T(i)*X_C(i)*Y_T(i)*Y_C(i)).^-0.5*Si_on(1,i)*Si_on(2,i)+(1/X_T(i)-1/X_C(i))*Si_on(1,i)+(1/Y_T(i)-1/Y_C(i))*Si_on(2,i);%Tsai_Wu
        %modified hashin theory
            e_fp(i)=(abs((Si_on(1,i)/X_T(i))^2+Si_on(3,i)/S(i)))^0.5;%fiber breaking
            e_fmp(i)=(abs(Si_on(1,i)/X_C(i)+Si_on(3,i)/S(i)))^0.5;%fiber_matrix shearing
            if e_fmp>1
                e_fm(i)=Si_on(1,i)/X_C(i);%fiber buckling
            end
            e_mp(i)=(abs((Si_on(2,i)/Y_T(i))^2+(Si_on(3,i)/S(i))^2))^0.5;%matrix in tension
            e_mm(i)=(abs((Si_on(2,i)/Y_C(i))^2+(Si_on(3,i)/S(i))^2))^0.5;%matrix in compression
    elseif abs(Ep_on(1))>0 || abs(Ep_on(2))>0 || abs(Ep_on(3))>0
        Si_on=C(:,:,i)*Ep_on;
    elseif abs(Si_on(1))>0 || abs(Si_on(2))>0 || abs(Si_on(3))>0
        Ep_on=S_on(:,:,i)*Si_on;
    elseif abs(Si_off(1))>0 || abs(Si_off(2))>0 || abs(Si_off(3))>0
        Si_on=T_sig(:,:,i)*Si_off;
        Ep_off=A\Si_off*h;
        Ep_on=S_on(:,:,i)*Si_on;
    end
end

clear E_x E_y E_s No_xy No_yx theta i h k ti M N h2 n1 a A