clc
clear all
format short e
m=1;
for  ep_counter=1:30
    e=2^(-2*(ep_counter-1));
 for N_counter =1:8
    N=8*2^(N_counter-1);
     sig=min(1/4,2*(sqrt(e/m))*log(N));
    
T1=zeros(N-1);
T2=zeros(N-1);
D1=zeros(N-1,1);
D2=zeros(N-1,1);
U=zeros(2*N+2,1);
U_exact=zeros(2*N+2,1);
W=zeros(N+1,1);
Z=zeros(N+1,1);
f1=zeros(N+1,1);
f2=zeros(N+1,1);
f=zeros(2*N+2,1);
a11=zeros(N+1,1);
a12=zeros(N+1,1);
a21=zeros(N+1,1);
a22=zeros(N+1,1);
x=zeros(N+1,1);
h=zeros(N,1);
% construction of mesh points along  spatial variable
for i=1:(N/4)
      h(i)=(4*sig)/N;
end

for i=(N/4)+1:3*N/4
    h(i)=(2*(1-2*(sig)))/N;
 end

for i=(3*N/4)+1:N
    h(i)=(4*sig)/N;
 end

x(1)=0;
 for i=1:N
    x(i+1)=x(i)+h(i);
 end
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
for i=1:N+1
    a11(i)=2*((x(i)+1)^2);
    a12(i)=-(1+(x(i)^3));
    a21(i)=-2*cos(pi*x(i)/4);
    a22(i)=2.2*exp(1-x(i));
end

for i=1:N+1
    f1(i)=2*exp(x(i));
    f2(i)=(10*x(i))+1;
end
W(1)=0;
W(N+1)=0;
Z(1)=0;
Z(N+1)=0;
%%%%%%%%%%%%construction of finite diffrence matrix T1  T2 D1 D2
%T1
for i=2:N-2
    T1(i,i-1)=-((2*e)/((h(i)+h(i+1))*h(i)));
    T1(i,i)=((2*e)/(h(i)*h(i+1)))+a11(i);
    T1(i,i+1)=-((2*e)/((h(i)+h(i+1))*h(i+1)));
end

i=1
T1(i,i)=((2*e)/(h(i)*h(i+1)))+a11(i);
T1(i,i+1)=-((2*e)/((h(i)+h(i+1))*h(i+1)));

i=N-1
T1(i,i-1)=-((2*e)/((h(i)+h(i+1))*h(i)));
T1(i,i)=((2*e)/(h(i)*h(i+1)))+a11(i);
%T2
for i=2:N-2
    T2(i,i-1)=-((2)/((h(i)+h(i+1))*h(i)));
    T2(i,i)=((2)/(h(i)*h(i+1)))+a22(i);
    T2(i,i+1)=-((2)/((h(i)+h(i+1))*h(i+1)));
end

i=1
T2(i,i)=((2)/(h(i)*h(i+1)))+a11(i);
T2(i,i+1)=-((2)/((h(i)+h(i+1))*h(i+1)));

i=N-1
T2(i,i-1)=-((2)/((h(i)+h(i+1))*h(i)));
T2(i,i)=((2)/(h(i)*h(i+1)))+a22(i);
%D1
for i=2:N-1
    D1(i,i-1)=0;
    D1(i,i)=a12(i);
    D1(i,i+1)=0;
end
D1(1,1)=a12(1);
D1(1,2)=0;
D1(N-1,N-2)=0;
D1(N-1,N-1)=a12(N-1);
%D2
for i=2:N-1
    D2(i,i-1)=0;
    D2(i,i)=a21(i);
    D2(i,i+1)=0;
end
D2(1,1)=a21(1);
D2(1,2)=0;
D2(N-1,N-2)=0;
D2(N-1,N-1)=a21(N-1);
%%%%%%%%%%%%%%%%%%%% Matrix using for iteration
P=inv(T2)*(D2)*inv(T1)*(D1);
Q=inv(T1)*(D1)*inv(T2)*(D2);
 end
end




    
    
