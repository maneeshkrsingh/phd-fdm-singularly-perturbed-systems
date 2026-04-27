clc
clear all
m1=1;
m2=2;
for N_counter=1
  N=2^(5+N_counter);

e=10^(-2);
mu=10^(-1);
sig2=min(1/2,2*mu*log(N));
sig1=1/4;
c=(1-exp(-1/e));
d=(1-exp(-1/mu));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(N+1,2);
U_num=zeros(2*N-2,1);
U1=zeros(N+1,1);
U_1=zeros(N-1,1);
V1=zeros(N+1,1);
V_1=zeros(N-1,1);
f1=zeros(N-1,1);
f2=zeros(N-1,1);
A=zeros(N-1,N-1);
B=zeros(N-1,N-1);
C=zeros(N-1,N-1);
D=zeros(N-1,N-1);
M=zeros(2*N-2,2*N-2);
x=zeros(N+1,1);
h=zeros(N,1);
%%%%%%%%%%%%%%%%%%%%%%%
% construction of mesh points along  spatial variable
for i=1:(N/4)
      h(i)=(4*sig1)/N;
end

for i=(N/4)+1:N/2
    h(i)=(4*((sig2)-(sig1)))/N;
 end
for i=(N/2)+1:N
    h(i)=(2*(1-(sig2)))/N;
end
for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
 U1(1)=0;
 U1(N+1)=0;
 
 V1(1)=0;
 V1(N+1)=0;
 
 for i=1:N-1
      f1(i)=((exp(-(x(i))/e)));
%     %      f1(i)=(-2*c)*(exp(-(x(i))/e))+(exp(-(x(i))/mu))*((d*e/(mu^2))-(d/mu)-d)+sin(pi*(x(i)/2))*((-e*(pi^2)/4)-4)-pi*cos(pi*(x(i)/2))+(x(i)*exp(x(i)-1));
% f2(i)=(exp(-(x(i))/mu))*(-(d/mu)-3*d)+c*(exp(-(x(i))/e))+exp(x(i)-1)*(mu*(x(i))+2*mu-6*x(i)-2)+2*sin(pi*(x(i)/2));
  f2(i)=2*sin(pi*(x(i)/2));
 end
 %%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for i=2:N-2
     A(i,i-1)=-2*e/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e/(h(i+1)*h(i)))+2+(1/h(i+1));
     A(i,i+1)=(-2*e/(h(i+1)+h(i)))-(1/h(i+1));
 end
 i=1;
      A(i,i)=(2*e/(h(i+1)*h(i)))+2+(1/h(i+1));
     A(i,i+1)=(-2*e/(h(i+1)+h(i)))-(1/h(i+1));
     
 i=N-1;   
     A(i,i-1)=-2*e/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e/(h(i+1)*h(i)))+2+(1/h(i+1));
     
   %D  
 for i=2:N-2
     D(i,i-1)=-2*mu/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*mu/(h(i+1)*h(i)))+2+(2/h(i+1));
     D(i,i+1)=(-2*mu/(h(i+1)+h(i)))-(2/h(i+1));
 end
 i=1;
      D(i,i)=(2*mu/(h(i+1)*h(i)))+2+(1/h(i+1));
     D(i,i+1)=(-2*mu/(h(i+1)+h(i)))-(2/h(i+1));
     
 i=N-1;   
     D(i,i-1)=-2*mu/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*mu/(h(i+1)*h(i)))+2+(2/h(i+1));
% B
for i=1:N-1
   b(i)=-1;
   B=diag(b);
   C=diag(b);
   
end


M=[A,B;C,D];
P=inv(M);
f=[f1;f2];
U_num=inv(M)*f;
for i=1:N-1
    U_1(i)=U_num(i);
    U_1;
end
 for i=N:2*N-2
     V_1((i-(N-1)))=U_num(i);
     V_1;
 end

 U1=[0;U_1;0];
  V1=[0;V_1;0];
  U=[U1,V1]';
end
