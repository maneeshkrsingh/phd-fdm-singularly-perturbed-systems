clc
clear all
format short e
C=1.05;
for ep_counter=8:8
    clear A
    clear B
    %clear C
    clear D
    e2=2^(-(ep_counter-1));
    e1=e2*(2^(-2*(ep_counter-1)));
for N_counter =1:2
    Nx1=16*2^(N_counter-1);
    del_t=0.1*2^(-(N_counter-1));
    Nt1=1/del_t;
   % Nt1=2*4^(N_counter-1);
 Nx1;
 Nt1;
% sig2
% sig1

for j=1:Nt1+1
    t(j)=(1*(j-1))/Nt1;
end
del_t=1/Nt1;

 for t_counter=1:Nt1+1
 for l=1:2
     if(l==1)
     Nx=Nx1;
    Nt=Nt1; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,Nt+1);
U_num=zeros(2*Nx-2,Nt+1);
U1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
%U_1=zeros(Nx-1,Nt);
U2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
%U_2=zeros(Nx-1,Nt);
f1=zeros(Nx-1,1);
f2=zeros(Nx-1,1);
A=zeros(Nx-1);
B=zeros(Nx-1);
%C=zeros(Nx-1);
D=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
d1=zeros(Nx-1,1);
d2=zeros(Nx-1,1);
b=zeros(2*Nx-2,1);
d=zeros(2*Nx-2,1);
flag=1;
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);
c=zeros(Nx-1,1);
%%%%%%%%%%%%%%%%%%%%%%%
% construction of mesh points along  spatial variable and time variable
x;

t;
del_t;
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
 for j=1:Nt+1
 for i=1:Nx+1
    U1(i,1)=0;
    U2(i,1)=0;    
end

    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;

% for j=1:Nt+1
    for i=1:Nx+1
        f1(i,j)=1;
        f2(i,j)=1;
    end
% end
 end
 x(1)=0;    
 for i=1:Nx
     x(i+1)=x(i)+(1/Nx);
 end 
        
 for i=2:Nx+1
     j=i-1;
     h(j)=x(i)-x(i-1);
 end 

while (flag==1)
    
        [U1,U2]=solution(e1,e2,x,Nx,Nt,del_t);
        [V1,V2]=smooth(e1,e2,x,Nx,Nt,del_t);
    

W_sin1=U1'-V1';
W_sin2=U2'-V2';
W1=W_sin1(round(1/del_t),:);
W2=W_sin2(round(1/del_t),:);
%           U1=U(1:Nx+1,:);
%           U2=U(Nx+2:2*Nx+2,:);
  for i=2:Nx+1
     j=i-1;
     h(j)=x(i)-x(i-1);
  end 
         a_dis=0;
       
        a_dis=h(1)*((abs(2*(((W1(3)-W1(2))/h(2))-((W1(2)-W1(1))/h(1)))/(h(2)+h(1))))^(1/2)+(abs(2*(((W2(3)-W2(2))/h(2))-((W2(2)-W2(1))/h(1)))/(h(2)+h(1))))^(1/2)); %h(i) has been changed for adjusting index
       
        a_dis=a_dis+h(Nx)*((abs(2*(((W1(Nx+1)-W1(Nx))/h(Nx))-((W1(Nx)-W1(Nx-1))/h(Nx-1)))/(h(Nx)+h(Nx-1))))^(1/2)+(abs(2*(((W2(Nx+1)-W2(Nx))/h(Nx))-((W2(Nx)-W2(Nx-1))/h(Nx-1)))/(h(Nx)+h(Nx-1))))^(1/2)); %h(i) has been changed for adjusting index
       
        for i=3:Nx
            a_dis=a_dis+h(i-1)*(((abs(2*(((W1(i)-W1(i-1))/h(i-1))-((W1(i-1)-W1(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2+((abs(2*(((W2(i)-W2(i-1))/h(i-1))-((W2(i-1)-W2(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2); %h(i) has been changed for adjusting index
        end
        M=zeros(Nx+1,1);
        M(1)=M(2);
        M(Nx+1)=M(Nx);
        a_dis;
        for i=2:Nx
            M(i)=a_dis+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2);
        end
        H=zeros(Nx,1);
        for i=2:Nx+1
            j=i-1;
            H(j)=((M(i-1)+M(i))/2)*h(j);           
        end
        L=zeros(Nx+1,1);
        L(1)=0;
        for i=1:Nx
            sum=0;
            for j=1:i
                sum=sum+H(j);
            end
            L(i+1)=sum;
        end
         C_tol=Nx*max(H)/L(Nx+1);
        if(C_tol<=C)
            flag=0;
            x_old=x;
        else
            for i=1:Nx+1
                Y(i)=(i-1)*L(Nx+1)/Nx;
            end
            new_x=interp1(L,x,Y);
            x=new_x;
        end
        %%%%%%%%%%%%%%%%%%%%%
          F1=U1;
          F2=U2;
 end
end
       if(l==2)
          for i=1:Nx+1
              x_double(2*i-1)=x(i);
          end
          for i=1:Nx
              x_double(2*i)=x(i)+0.5*(x(i+1)-x(i));
          end
          Nx=2*Nx;
                  
        for i=2:Nx+1
            j=i-1;
            h_double(j)=x_double(i)-x_double(i-1);
        end 
          [U1,U2]=solution(e1,e2,x_double,Nx,2*Nt,del_t/2); 
          G1=U1;
          G2=U2;
       end
     end
 end 
  err1=abs(F1-G1(1:2:Nx+1,1:2:2*Nt+1));
    err_sup1=max(err1);
         max_err1=max(err_sup1);
       error1(ep_counter,N_counter)=max_err1;
          
          err2=abs(F2-G2(1:2:Nx+1,1:2:2*Nt+1));
        err_sup2=max(err2);
         max_err2=max(err_sup2);
       error2(ep_counter,N_counter)=max_err2;
          
end
end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rt1(j,i)=log2(error1(j,i)/error1(j,i+1));
    end
end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rt2(j,i)=log2(error2(j,i)/error2(j,i+1));
    end
end
   error1
   max(error1)
    error2
    max(error2)
    cnv_rt1
    max(cnv_rt1)
    cnv_rt2
    max(cnv_rt2)