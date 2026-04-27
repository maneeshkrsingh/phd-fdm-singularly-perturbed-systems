clc
clear all
format long
C=1.1;
for ep_counter=1:1
    e=10^(-(ep_counter));
    m1=exp(-1/e);
    m2=1-m1;
    c=1/(1-exp(-1/e));
    %e1=e*(2^(-2*(ep_counter-1)));
%     figure
%     hold on
for N_counter =3:3
    disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_counter)]);
    Nx=2*2^(N_counter+1);
    Nt=Nx/2;
%     del_t=1/Nt;
   % k=0.1^(-(N_counter));
    k=1/Nt;
   %Nt=1/k;
   % Nt=2*4^(N_counter-1);
%sig1=min(1/2,2*e*log(Nx));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,Nt+1);
%U_num=zeros(2*Nx-2,Nt+1);
U1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
%U_1=zeros(Nx-1,Nt);
U2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
U1_exct=zeros(Nx+1,Nt+1);
U2_exct=zeros(Nx+1,Nt+1);
U1_final=zeros(Nx+1,Nt+1);
U2_final=zeros(Nx+1,Nt+1);
U1_intl=zeros(Nx+1,1);
U2_intl=zeros(Nx+1,1);

V1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
%U_1=zeros(Nx-1,Nt);
V2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
V1_exct=zeros(Nx+1,Nt+1);
V2_exct=zeros(Nx+1,Nt+1);
V1_final=zeros(Nx+1,Nt+1);
V2_final=zeros(Nx+1,Nt+1);
V1_intl=zeros(Nx+1,1);
V2_intl=zeros(Nx+1,1);

W1_sin=zeros(Nx+1,Nt+1);
W2_sin=zeros(Nx+1,Nt+1);
Y=zeros(Nx+1,1);
%M=zeros(Nx+1,1);
H=zeros(Nx,1);
L=zeros(Nx+1,1);
%U_2=zeros(Nx-1,Nt);
% f1=zeros(Nx-1,1);
% f2=zeros(Nx-1,1);
% A=zeros(Nx-1);
% B=zeros(Nx-1);
% C=zeros(Nx-1);
% D=zeros(Nx-1);
% b1=zeros(Nx-1,1);
% b2=zeros(Nx-1,1);
% b=zeros(2*Nx-2,1);

h=zeros(Nx,1);
t=zeros(Nt+1,1);

for j=1:Nt+1
    t(j)=(j-1)/Nt;
end


%x_old=zeros(Nt+1,Nx+1);

     flag=1;

counter=0;    

x=zeros(1,Nx+1);
 x(1)=0;    
 for i=1:Nx
     x(i+1)=x(i)+(1/Nx);
 end 
 

        
 for i=2:Nx+1
     j=i-1;
     h(j)=x(i)-x(i-1);
 end 
 %P(:,1)=x;
while (flag==1)

%        P(:,counter+1)=x;
%        counter=counter+1;
    for p=1:Nt    
        if(p==1)
            for i=1:Nx+1
                U1_intl(i)=m1+m2*(1-x(i))-exp(-x(i)/e);
                U2_intl(i)=c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1);
                V1_intl(i)=m1+m2*(1-x(i))-exp(-x(i)/e);
                V2_intl(i)=c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1);   
                
           U1_final(:,1)=U1_intl;
           U2_final(:,1)=U2_intl;
           
           V1_final(:,1)=V1_intl;
           V2_final(:,1)=V2_intl;
    %            P(:,1)=x;
            end

    %x_old(1,:)=x;
    %P(:,p)=x;
        else
            U1_intl=U1;
            U2_intl=U2;
            V1_intl=V1;
            V2_intl=V2;            
%     p
%     disp(['Hi ']);
    %x_old(p+1,:)=x;
        end

        [U1,U2]=solution(e,x,Nx,Nt,U1_intl,U2_intl);
        [V1,V2]=smooth(e,x,Nx,Nt,V1_intl,V2_intl);
        
  

            U1_final(:,p+1)=U1;
            U2_final(:,p+1)=U2;
            

            V1_final(:,p+1)=V1;
            V2_final(:,p+1)=V2;     
        

    end
                    W1_sin=U1_final-V1_final;
                    W2_sin=U2_final-V2_final;
   
%                     w1_sin=u1_sol-v1_smo;
%                     w2_sin=u2_sol-v2_smo;

                    W1=W1_sin(:,round((Nt)));
                    W2=W2_sin(:,round((Nt)));
                    N=Nx;
                    for i=2:Nx+1
                         j=i-1;
                         h(j)=x(i)-x(i-1);
                    end 
                    a_dis=0;
       
                     a_dis=h(1)*((abs(2*(((W1(3)-W1(2))/h(2))-((W1(2)-W1(1))/h(1)))/(h(2)+h(1))))^(1/2)+(abs(2*(((W2(3)-W2(2))/h(2))-((W2(2)-W2(1))/h(1)))/(h(2)+h(1))))^(1/2)); %h(i) has been changed for adjusting index
       
                    a_dis=a_dis+h(N)*((abs(2*(((W1(N+1)-W1(N))/h(N))-((W1(N)-W1(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/2)+(abs(2*(((W2(N+1)-W2(N))/h(N))-((W2(N)-W2(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/2)); %h(i) has been changed for adjusting index
       
                    for i=3:N
                        a_dis=a_dis+h(i-1)*(((abs(2*(((W1(i)-W1(i-1))/h(i-1))-((W1(i-1)-W1(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2+((abs(2*(((W2(i)-W2(i-1))/h(i-1))-((W2(i-1)-W2(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2); %h(i) has been changed for adjusting index
                    end
                    M=zeros(N+1,1);

                    a_dis;
                    for i=2:N
                      M(i)=a_dis+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2);
                    end
                    M(1)=M(2);
                    M(N+1)=M(N);
                    H=zeros(N,1);
                    for i=2:N+1
                     j=i-1;
                     H(j)=((M(i-1)+M(i))/2)*h(j);           
                    end
                     L=zeros(N+1,1);
                     L(1)=0;
                     for i=1:N
                       sum=0;
                          for j=1:i
                          sum=sum+H(j);
                          end
                          L(i+1)=sum;
                     end
                     C_tol=N*max(H)/L(N+1);
                     if(C_tol<=C)
                        flag=0;
                        x_old=x;
                     else
                        for i=1:N+1
                            Y(i)=(i-1)*L(N+1)/N;
                        end
                        new_x=interp1(L,x,Y);
                        x=new_x;
                     end
end
 

 for j=1:Nt+1
    for i=1:Nx+1
       U1_exct(i,j)=exp(-t(j))*(m1+m2*(1-x(i))-exp(-x(i)/e));
       U2_exct(i,j)=exp(-t(j))*(c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1));
       
    end 
 end

 err1=abs(U1_final-U1_exct);
 err2=abs(U2_final-U2_exct);



err1_sup=max(err1);
max_err1=max(err1_sup);
error1(ep_counter,N_counter)=max_err1;

err2_sup=max(err2);
max_err2=max(err2_sup);
error2(ep_counter,N_counter)=max_err2;   


% for i=1:Nx+1
%     plot(P(i,:),0:counter-1);
% end

end
% hold off
end
% for j=1:ep_counter
%     for i=1:N_counter-1
%          conv1(j,i)=log2(error1(j,i)/error1(j,i+1));
%          conv2(j,i)=log2(error2(j,i)/error2(j,i+1));
%     end
% end
% error1
% error2
% conv1
% conv2
% figure1
 surf(x,t,U1_final')
%figure2
%surf(x,t,U1_exct')