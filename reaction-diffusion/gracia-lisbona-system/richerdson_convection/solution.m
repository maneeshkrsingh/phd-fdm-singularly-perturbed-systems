
function [U1,U2]= solution(e1,e2,x,Nx,k,u1_d,u2_d,time_counter)   
A=zeros(Nx-1);
B=zeros(Nx-1);
C=zeros(Nx-1);
D=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
d=zeros(2*Nx-2,1);
c=zeros(Nx-1,1);
U1=zeros(1,Nx+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
U2=zeros(1,Nx+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
for i=2:Nx+1
     j=i-1;
     h(j)=x(i)-x(i-1);
 end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
% for j=1:Nt+1

    U1(1,1)=0;
    U1(1,Nx+1)=0;
    

% for j=1:Nt+1
    for i=1:Nx+1
        f1(i)=1;
        f2(i)=1;
    end
% end
%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for i=2:Nx-2
     A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
      A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
%D
 for i=2:Nx-2
     D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
      D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     
  for i=1:Nx-1
   c(i,1)=-k;
  end
     B=diag(c);
   C=diag(c);
   M=[A,B;C,D];
   
   
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    for j=2:Nt 
     for i=2:Nx-2
             b1(i)=(k+u1_d(1,i+1));
             b2(i)=(k+u2_d(1,i+1));
     end
     
      b1(1)=k+u1_d(1,2);
      b2(1)=k+u2_d(1,2);
      b1(Nx-1)=k+u1_d(1,Nx);
      b2(Nx-1)=k+u2_d(1,Nx);
      b=[b1;b2];
     
      d=M\b;   
      %d=inv(M)*b;
      for i=2:(Nx)
          U1(1,i)=d(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          U2(1,i-Nx+1)=d(i-1);
      end
 %   end
    %U1(2:Nx,:)=U_num(1:Nx-1,:);
    %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:);
   U=[U1;U2];
end