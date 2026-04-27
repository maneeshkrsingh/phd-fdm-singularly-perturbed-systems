clc
clear all
format short e
m=1;
for  ep_counter=1:30
    e(ep_counter)=2^(-2*(ep_counter-1));
 for N_counter =1:6
    N=8*2^(N_counter-1);
    
    sig=min(1/4,2*(sqrt(e(ep_counter)/m))*log(N));
 U1=zeros(2,N+1);  
 W_new=zeros(N+1,1);
Z_new=zeros(N+1,1);
for flag=1:2     
        
T1=zeros(N-1);
T2=zeros(N-1);
D1=zeros(N-1,1);
D2=zeros(N-1,1);
U=zeros(2,N+1);
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
W1=zeros(N-1,1);
Z1=zeros(N-1,1);
W2=zeros(N-1,1);
Z2=zeros(N-1,1);
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


 for i=1:N
     x(1)=0;
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
    f_1=f1(2:N,1);
    f_2=f2(2:N,1);
end


%     W1(1)=0;
%     W1(N-1)=0;
for i=1:N-1
    W1(i,1)=0.1;
end

%     Z1(1)=0;
%     Z1(N-1)=0;
 for i=1:N-1
    Z1(i,1)=0.1;
end   
W(1)=0;
W(N+1)=0;
Z(1)=0;
Z(N+1)=0;
%%%%%%%%%%%%construction of finite diffrence matrix T1  T2 D1 D2
%T1
for i=2:N-2
    T1(i,i-1)=-((2*e(ep_counter))/((h(i)+h(i+1))*h(i)));
    T1(i,i)=((2*e(ep_counter))/(h(i)*h(i+1)))+a11(i);
    T1(i,i+1)=-((2*e(ep_counter))/((h(i)+h(i+1))*h(i+1)));
end

i=1;
T1(i,i)=((2*e(ep_counter))/(h(i)*h(i+1)))+a11(i);
T1(i,i+1)=-((2*e(ep_counter))/((h(i)+h(i+1))*h(i+1)));

i=N-1;
T1(i,i-1)=-((2*e(ep_counter))/((h(i)+h(i+1))*h(i)));
T1(i,i)=((2*e(ep_counter))/(h(i)*h(i+1)))+a11(i);
%T2
for i=2:N-2
    T2(i,i-1)=-((2)/((h(i)+h(i+1))*h(i)));
    T2(i,i)=((2)/(h(i)*h(i+1)))+a22(i);
    T2(i,i+1)=-((2)/((h(i)+h(i+1))*h(i+1)));
end

i=1;
T2(i,i)=((2)/(h(i)*h(i+1)))+a11(i);
T2(i,i+1)=-((2)/((h(i)+h(i+1))*h(i+1)));

i=N-1;
T2(i,i-1)=((2)/((h(i)+h(i+1))*h(i)));
T2(i,i)=((2)/(h(i)*h(i+1)))+a22(i);
%D1
for i=2:N-2
    D1(i,i-1)=0;
    D1(i,i)=a12(i);
    D1(i,i+1)=0;
end
D1(1,1)=a12(1);
D1(1,2)=0;
D1(N-1,N-2)=0;
D1(N-1,N-1)=a12(N-1);
%D2
for i=2:N-2
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
R=inv(T2)*((f_2)-((D2)*inv(T1)*(f_1)));
S=inv(T1)*((f_1)-((D1)*inv(T2)*(f_2)));
%%%%%%%%%%%%%%%%%%%%% Iterative Scheme
if(max(abs(W2-W1))>10^(-5) || max(abs(Z2-Z1))>10^(-5))
W2=P*W1+R;
W1=W2;
Z2=Q*Z1+S;
Z1=Z2;
end
W(2:N)=W2;
Z(2:N)=Z2;

U=[W';Z'];
% 
if (flag==1)
   U1=U;
   W_new=W;
   Z_new=Z;
   x1=x;
    N=8*N;
end 

if (flag==2)
    W_interp=interp1(x1,W_new',x,'linear');
    Z_interp=interp1(x1,Z_new',x,'linear');
    U_interp=[W_interp';Z_interp'];
    error=max(abs(U_interp-U));
    maxerr(ep_counter,N_counter)=max(error);
end

end
end
end




    
    
