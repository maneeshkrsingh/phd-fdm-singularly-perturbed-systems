function [V1,V2]= smooth(e1,e2,x,Nx,k,v1_d,v2_d,time_counter)

A=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
c=zeros(Nx-1,1);
V1=zeros(1,Nx+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
V2=zeros(1,Nx+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 


    V1(1,1)=0;
  %  V1(Nx+1,j)=0;
    
    V2(1,1)=0;
%    V2(Nx+1,j)=0;

% for j=1:Nt+1
    for i=1:Nx+1
        f1(i)=1;
        f2(i)=1;
    end
% end
%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for i=2:Nx-2
     A(i,i-1)=0;
     A(i,i)=k+1;
     A(i,i+1)=0;
 end
 i=1;
      A(i,i)=k+1;
     A(i,i+1)=0;
     
 i=Nx-1;   
     A(i,i-1)=0;
     A(i,i)=k+1;
%D
 for i=2:Nx-2
     D(i,i-1)=0;
     D(i,i)=k+1;
     D(i,i+1)=0;
 end
 i=1;
      D(i,i)=k+1;
     D(i,i+1)=0;
     
 i=Nx-1;   
     D(i,i-1)=0;
     D(i,i)=k+1;
     
  for i=1:Nx-1
   c(i,1)=-k;
  end
     B=diag(c);
   C=diag(c);
   M=[A,B;C,D];
   
   
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    for j=2:Nt 
     for i=2:Nx-2
             b1(i)=(k+v1_d(1,i));
              b2(i)=(k+v2_d(1,i));
     end
     
      b1(1)=k;
      b2(1)=k;
      b1(Nx-1)=k+v1_d(1,Nx-1);
      b2(Nx-1)=k+v2_d(1,Nx-1);
      b=[b1;b2];
     
         
      d=inv(M)*b;
      for i=2:(Nx)
          V1(1,i)=d(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          V2(1,i-Nx+1)=d(i-1);
      end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
