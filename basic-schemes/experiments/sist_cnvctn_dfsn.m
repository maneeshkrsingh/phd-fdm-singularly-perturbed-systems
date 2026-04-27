clc
clear all
format short e
m1=1;
m2=1;
for  ep_counter=1:2
   e1=2^(-2*(ep_counter-1));
   e2=2^(-2*(ep_counter-1));

  
for N_counter =1:3
    N=2^(6+(N_counter-1));
    
sig2=min(4*max(e1/m1,e2/m2)*log(N),1/2);
sig1=min(2*min(e1/m1,e2/m2)*log(N),(sig2)/2);

T1=zeros(N-1);
T2=zeros(N-1);

U1=zeros(N+1,1);
U2=zeros(N+1,1);
f1=zeros(N+1,1);
f2=zeros(N+1,1);

a11=zeros(N+1,1);
a21=zeros(N+1,1);
a22=zeros(N+1,1);
x=zeros(N+1,1);
h=zeros(N,1);

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
h_1=h(2:N,1);
 for i=1:N
     x(1)=0;
    x(i+1)=x(i)+h(i);
 end
 % numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
for i=1:N+1
    a11(i)=3;
    a21(i)=2.75;
    a22(i)=2;
end

for i=1:N+1
    f1(i)=15*(x(i)^4);
    f2(i)=0.6*exp(x(i)); 
     f_1=f1(2:N,1);
      f_2=f2(2:N,1);
end
% %%%%%%%%%%%construction of finite diffrence matrix T1  T2 D1 D2
% T1
for i=2:N-2
    T1(i,i-1)=-((2*e1)/((h(i)+h(i+1))*h(i)));
    T1(i,i)=((2*e1)/(h(i)*h(i+1)))+3/h(i+1);
    T1(i,i+1)=-((2*e1)/((h(i)+h(i+1))*h(i+1)))-3/h(i+1);
end
i=1;
    T1(i,i)=((2*e1)/(h(i)*h(i+1)))+3/h(i+1);
    T1(i,i+1)=-((2*e1)/((h(i)+h(i+1))*h(i+1)))-3/h(i+1);
i=N-1;
    T1(i,i-1)=-((2*e1)/((h(i)+h(i+1))*h(i)));
    T1(i,i)=((2*e1)/(h(i)*h(i+1)))+3/h(i+1);
    
 
      U1(1)=0;
     U1(N+1)=0;

     U1(2:N)=(inv(T1))*(f_1);
 %T2    
for i=2:N-2
    T2(i,i-1)=-((2*e2)/((h(i)+h(i+1))*h(i)));
    T2(i,i)=((2*e2)/(h(i)*h(i+1)))+2/h(i+1);
    T2(i,i+1)=-((2*e2)/((h(i)+h(i+1))*h(i+1)))-2/h(i+1);
end
i=1;
    T2(i,i)=((2*e2)/(h(i)*h(i+1)))+3/h(i+1);
    T2(i,i+1)=-((2*e2)/((h(i)+h(i+1))*h(i+1)))-2/h(i+1);
i=N-1;
    T2(i,i-1)=-((2*e2)/((h(i)+h(i+1))*h(i)));
    T2(i,i)=((2*e2)/(h(i)*h(i+1)))+2/h(i+1);
    
    U_1=abs(U1(2:N)-U1(1:N-1));
 
 D(1:N-1,1)= ((f_2)-(2.75/(h_1))*(U_1));
 
  U2(2:N)=(inv(T2))*D(1:N-1,1);
 
 

end
end
    
    
    
    
    