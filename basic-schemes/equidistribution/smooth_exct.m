function [V1,V2]= smooth(e,x,Nx,k,v1_d,v2_d,time_counter)
A=zeros(Nx-1);
B=zeros(Nx-1);
C=zeros(Nx-1);
D=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
c=zeros(Nx-1,1);
V1=zeros(1,Nx+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
V2=zeros(1,Nx+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
  
   V1(Nx+1,1)=0;
  
    
    V2(Nx+1,1)=0;
for i=1:Nx+1
     
        f1(i)=(exp(-x(i)/e)*(-(x(i)/e)-2*x(i)+c)+x(i)*exp(x(i)-1)+(2*m1*x(i)+m2*(3*x(i)-2*(x(i)^2)+1)-c));
        f2(i)=(exp(-x(i)/e)*(-2*c*(x(i)/e)+x(i)-c*x(i))+exp(x(i)-1)*((1+2*x(i))*(1+x(i))+e*(x(i)+2)-(x(i)^2))+(x(i)*c-m1*x(i)-m2*x(i)*(1-x(i))));
end

%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 
 for i=2:Nx-2
     A(i,i-1)=0;
     A(i,i)=k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;
     A(i,i+1)=-k*(1+x(i))/h(i+1);
 end
 i=1;
      A(i,i)=k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;
     A(i,i+1)=-k*(1+x(i))/h(i+1);
     
 i=Nx-1;   
    A(i,i-1)=0;
     A(i,i)=k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;

 for i=2:Nx-2
     D(i,i-1)=0;
     D(i,i)=k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     D(i,i+1)=-k*(1+2*x(i))/h(i+1);
 end
 i=1;
      D(i,i)=k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     D(i,i+1)=-k*(1+2*x(i))/h(i+1);
     
 i=Nx-1;   
     D(i,i-1)=0;
     D(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     

  for i=2:Nx-2
     B(i,i-1)=0;
     B(i,i)=-k;
     B(i,i+1)=0;
 end
 i=1;
     B(i,i)=-k;
     B(i,i+1)=0;
     
 i=Nx-1;   
     B(i,i-1)=0;
     B(i,i)=-k;
     
for i=2:Nx-2
     C(i,i-1)=0;
     C(i,i)=-k*(x(i));
     C(i,i+1)=0;
 end
 i=1;
     C(i,i)=-k*(x(i));
     C(i,i+1)=0;
     
 i=Nx-1;   
     C(i,i-1)=0;
     C(i,i)=-k*(x(i));  
     

   M_mat=[A,B;C,D];
  
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     for i=2:Nx-2
             b1(i)=(k*f1(i)+V1_intl(i+1));
              b2(i)=(k*f2(i)+V2_intl(i+1));
     end
     
      b1(1)=(k*f1(1)+V1_intl(2));
      b2(1)=(k*f2(1)+V2_intl(2));
      b1(Nx-1)=k*f1(Nx-1)+V1_intl(Nx);
      b2(Nx-1)=k*f2(Nx-1)+V2_intl(Nx);
      b=[b1;b2];
     
      P=M_mat\b;
  %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for i=2:(Nx)
          V1(i)=P(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          V2(i-Nx+1)=P(i-1);
      end

end