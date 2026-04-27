  clc
  clear all
  format short e 

  e=2^(-4);
 % bt1=1;
 % bt2=1;
  m1=-exp(-1/e);
m2=-1-m1;
 % c1=1;
% c2=exp(-1/sqrt(e))-1;

for l1=1:3
    Nx=16*2^(l1-1);
    Ny=16*2^(l1-1);
    Nt=4*2^(l1-1);
  sigx=min(1/4,4*e*log(Nx));
  sigy=min(1/4,4*e*log(Ny));

 % for l2=1:1
     % Nx=Nx1*2^(l2-1);
     % Ny=Ny1*2^(l2-1);
     % Nt=Nt1*2^(l2-1);
  
 x=zeros(Nx+1,1);
 y=zeros(Ny+1,1);
 t=zeros(Nt+1,1);
 t1=zeros((2*Nt)+1,1);
 hy=zeros(Ny,1);
 hx=zeros(Nx,1);
 u1_exct=zeros(Nx+1,Ny+1,Nt+1);
 u2_exct=zeros(Nx+1,Ny+1,Nt+1);
 u1=zeros(Nx+1,Ny+1,Nt+1);
 u2=zeros(Nx+1,Ny+1,Nt+1);
 u1_temp=zeros(Nx+1,Ny+1,(2*Nt)+1);
 u2_temp=zeros(Nx+1,Ny+1,(2*Nt)+1);
 b11=zeros(Nx-1,1);
 b12=zeros(Nx-1,1);
 b21=zeros(Nx-1,1);
 b22=zeros(Nx-1,1);
 d11=zeros(Nx-1,1);
 d12=zeros(Nx-1,1);
 d21=zeros(Nx-1,1);
 d22=zeros(Nx-1,1);
 A1=zeros(Nx-1,Nx-1);
 B1=zeros(Nx-1,Nx-1);
 C1=zeros(Nx-1,Nx-1);
 D1=zeros(Nx-1,Nx-1);
 A2=zeros(Ny-1,Ny-1);
 B2=zeros(Ny-1,Ny-1);
 C2=zeros(Ny-1,Ny-1);
 D2=zeros(Ny-1,Ny-1);
 M1=zeros(2*Nx-2,2*Nx-2);
 M2=zeros(2*Ny-2,2*Ny-2);
 P1=zeros(2*Nx-2,1);
 P2=zeros(2*Ny-2,1);
 
 
  %%%%%%%%%%%%%%%%%mesh construction%%%%%%%%%%%%
for k=1:Nt+1
    t(k)=((k-1))/Nt;
end
  ht=1/(Nt);
  ht1=ht/2;
  for k=1:(2*Nt)+1 
      t1(k)=(k-1)/(2*Nt); % mesh t1 and t1/2
  end
%%%%%%%%%%%%%%%%%%%%%% x direction %%%%%%%%
for i=1:(Nx/2)
    hx(i)=2*(1-sigx)/Nx;
end

for i=(Nx/2)+1:Nx
    hx(i)=2*sigx/Nx;
end

x(1)=0;
for i=1:Nx-1
     x(i+1)=x(i)+hx(i);
     
 end
x(Nx+1)=1;
%%%%%%%%%%%%%%%%%%%%%% y direction %%%%%%%%

for j=1:(Ny/2)
    hy(j)=2*(1-sigy)/Ny;
end

for j=(Ny/2)+1:Ny
    hy(j)=2*sigy/Ny;
end

y(1)=0;
for j=1:Ny-1
     y(j+1)=y(j)+hy(j);
end
 y(Ny+1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%intial conditions and boundary conditon for both steps
%%%%%%%%%%%%%%%%%%%%%

for i=1:Nx+1
    for j=1:Ny+1
        
        %u1(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));
        %u2(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));
        u1(i,j,1)=0;
        u2(i,j,1)=0;
        u1_temp(i,j,1)=0;
        u2_temp(i,j,1)=0;
        %u1_temp(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));
        %u2_temp(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));

    end
end

%%%%%%%%%%%%%%%%%%%%%boundary conditon for both steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:Ny+1
    for k=1:Nt+1

           u1(1,j,k)=0;
           u2(1,j,k)=0;
           u1(Nx+1,j,k)=0;
           u2(Nx+1,j,k)=0;
    end
end

for i=1:Nx+1
    for k=1:Nt+1

         u1(i,1,k)=0;
         u2(i,1,k)=0;
         u1(i,Ny+1,k)=0;
         u2(i,Ny+1,k)=0;
    end
end

%%%%%%%%%%%%%%
for j=1:Ny+1
    for k=1:(2*Nt)+1
           u1_temp(1,j,k)=0;
           u2_temp(1,j,k)=0;
           u1_temp(Nx+1,j,k)=0;
           u2_temp(Nx+1,j,k)=0;
    end
    
end

for i=1:Nx+1
    for k=1:(2*Nt)+1

           u1_temp(i,1,k)=0;
           u2_temp(i,1,k)=0;
           u1_temp(i,Ny+1,k)=0;
           u2_temp(i,Ny+1,k)=0;
    end
end

%%%%%%%%%%%%%%%%%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for k=1:Nt+1
  for j=1:Ny+1
     for i=1:Nx+1
        %U2_exct(i,j)=exp(-t(j))*(m1+m2*(1-x(i))-exp(-x(i)/e2));
        u1_exct(i,j,k)=(1+t(k)-exp(-t(k)))*(m1+m2*x(i)+exp(-(1-x(i))/e))*(m1+m2*y(j)+exp(-(1-y(j))/e));
        u2_exct(i,j,k)=(1+((t(k))^2)-exp(-2*t(k)))*(m1+m2*x(i)+exp(-(1-x(i))/e))*(m1+m2*y(j)+exp(-(1-y(j))/e));
    end
 end
 end


%%%%%%%%%%%%%%%%%% matrix construction %%%%%%%%%%%%%%%%%%%%


it=0;
% 
for k=2:(Nt)+1
% for k=2:6
it;

    for j=1:Ny-1
        
      for i=2:Nx-2

        A1(i,i-1)=-2*e*ht/((hx(i)+hx(i+1))*hx(i))-ht*(1+x(i)*(1-x(i)))/hx(i);
        A1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+x(i)*(1-x(i)))/hx(i);  % b11=b11/2+b11/2 as defined for fractional step
        A1(i,i+1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i+1));
       
     end
     i=1;
        A1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+x(i)*(1-x(i)))/hx(i); 
        A1(i,i+1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i+1));

     i=Nx-1; 
         A1(i,i-1)=-2*e*ht/((hx(i)+hx(i+1))*hx(i))-ht*(1+x(i)*(1-x(i)))/hx(i);
        A1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+x(i)*(1-x(i)))/hx(i); 
        
  %%%%%%%%%%%%%%%%%%%%%%   for second component equation %%%%%%%%%%%%%%%%%
    for i=2:Nx-2

        D1(i,i-1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i))-ht*(1+(x(i)^2)*(1-x(i)))/hx(i);
        D1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+(x(i)^2)*(1-x(i)));  % b22=b22/2+b22/2 as defined for fractional step
        D1(i,i+1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i+1));
       
     end
     i=1;
        D1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+(x(i)^2)*(1-x(i)));  
        D1(i,i+1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i+1));

     i=Nx-1; 
          D1(i,i-1)=-(2*e*ht)/((hx(i)+hx(i+1))*hx(i))-ht*(1+(x(i)^2)*(1-x(i)));
        D1(i,i)=1+(2*e*ht)/(hx(i)*hx(i+1))+ht*(1+(x(i)^2)*(1-x(i)));
%%%%%%%%%%%%%%%%%%for second component in first equation %%%%%%%%%%%%%%%%%%%%
 for i=2:Nx-2
     B1(i,i-1)=0;
     B1(i,i)=0;%a12=b12/2+b12/2 as defined for fractional step
     B1(i,i+1)=0;
 end
 i=1;
     B1(i,i)=0;
     B1(i,i+1)=0;
     
 i=Nx-1;   
     B1(i,i-1)=0;
     B1(i,i)=0;
%%%%%%%%%%%%%%%%%%for first component in second equation %%%%%%%%%%%%%%%%%%%%
 for i=2:Nx-2
     C1(i,i-1)=0;
     C1(i,i)=0;%b21=b21/2+a21/2 as defined for fractional step
     C1(i,i+1)=0;
 end
 i=1;
     C1(i,i)=0;
     C1(i,i+1)=0;
     
 i=Nx-1;   
     C1(i,i-1)=0;
     C1(i,i)=0;
     
      M1=[A1,B1;C1,D1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%left hand side%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=2:Nx-2
           
         b11(i)=ht*func11(x(i),y(j),t(k),e)+u1(i+1,j,k-1);
         b12(i)=ht*func12(x(i),y(j),t(k),e)+u2(i+1,j,k-1);
    end

    
         b11(1)=ht*func11(x(1),y(j),t(k),e)+u1(2,j,k-1);
         b12(1)=ht*func12(x(i),y(j),t(k),e)+u2(2,j,k-1);
        
         b11(Nx-1)=ht*func11(x(Nx-1),y(j),t(k),e)+u1(Nx,j,k-1);
         b12(Nx-1)=ht*func12(x(Nx-1),y(j),t(k),e)+u2(Nx,j,k-1);
         
         b1=[b11;b12];
       
          
         P1=M1\b1;

       for i=2:Nx
           u1_temp(i,j+1,2*(k-1))=P1(i-1); %storing the value of P1(each time level) in column in u1.
      end
    
      for i=Nx+1:2*Nx-1
          u2_temp(i+1-Nx,j+1,2*(k-1))=P1(i-1);     
      end
         
         
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    for i=1:Nx-1
      for j=2:Ny-2

        A2(j,j-1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j))-ht*(1+y(j)*(1-y(j)))/hy(j);
        A2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+y(j)*(1-y(j)))/hy(j);  % b11=b11/2+a11/2 as defined for fractional step
        A2(j,j+1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j+1));
       
     end
     j=1;
         A2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+y(j)*(1-y(j)))/hy(j);
        A2(j,j+1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j+1));
        
     j=Ny-1; 
        A2(j,j-1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j))-ht*(1+y(j)*(1-y(j)))/hy(j);
        A2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+y(j)*(1-y(j)))/hy(j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%for second component%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=2:Ny-2

        D2(j,j-1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j))-ht*(1+(y(j)^3)*(1-y(j)))/hy(j);
        D2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+(y(j)^3)*(1-y(j)))/hy(j);  % b22=b22/2+b22/2 as defined for fractional step
        D2(j,j+1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j+1));
       
 end
     j=1;
        D2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+(y(j)^3)*(1-y(j)))/hy(j);
        D2(j,j+1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j+1));
        
     j=Ny-1; 
         D2(j,j-1)=-(2*e*ht)/((hy(j)+hy(j+1))*hy(j))-ht*(1+(y(j)^3)*(1-y(j)))/hy(j);
        D2(j,j)=1+(2*e*ht)/(hy(j)*hy(j+1))+ht*(1+(y(j)^3)*(1-y(j)))/hy(j);
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%second component in first equation %%%%%%%%%%%%%%%%%%
 
 for j=2:Ny-2
     B2(j,j-1)=0;
     B2(j,j)=0;%a12=a12/2+a12/2 as defined for fractional step
     B2(j,j+1)=0;
     
 end
   j=1;
     B2(j,j)=0;
     B2(j,j+1)=0;
     
   j=Ny-1; 
     B2(j,j-1)=0;
     B2(j,j)=0;
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%first component in  second equation %%%%%%%%%%%%%%%%%% 
  for j=2:Ny-2
     C2(j,j-1)=0;
     C2(j,j)=0;%a21=a21/2+a21/2 as defined for fractional step
     C2(j,j+1)=0;
     
 end
   j=1;
     C2(j,j)=0;
     C2(j,j+1)=0;
     
   j=Ny-1; 
     C2(j,j-1)=0;
     C2(j,j)=0;
  
    M2=[A2,B2;C2,D2];  
%%%%%%%%%%%%%%%%%%%%%%%left hand side %%%%%%%%%%%%%%%%%%%%%%%%     
     
    for j=2:Ny-2
           
         b21(j)=ht*func21(x(i),y(j),t(k),e)+u1_temp(i,j+1,2*(k-1));
         b22(j)=ht*func22(x(i),y(j),t(k),e)+u2_temp(i,j+1,2*(k-1));
    end

         b21(1)=ht*func21(x(i),y(1),t(k),e)+u1_temp(i,2,2*(k-1));
         b22(1)=ht*func22(x(i),y(1),t(k),e)+u2_temp(i,2,2*(k-1));
         
         b21(Ny-1)=ht*func21(x(i),y(Ny-1),t(k),e)+u1_temp(i,Ny,2*(k-1));
         b22(Ny-1)=ht*func22(x(i),y(Ny-1),t(k),e)+u2_temp(i,Ny,2*(k-1));
   
           b2=[b21;b22];
 
           P2=M2\b2;

          
         
      for j=2:Ny
           u1(i+1,j,k)=P2(j-1); %storing the value of P1(each time level) in column in u1.
      end
  
      for j=Ny+1:2*Ny-1
          u2(i+1,j+1-Ny,k)=P2(j-1);     
      end
    end
%      
  it=it+1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end
%      if(l2==1)
%          F1=u1;
%          F2=u2;
%      end
%      if(l2==2)
%          G1=u1;
%          G2=u2;
%      end



 
 % end
    max_err1(l1)=max(max(max(abs(u1- u1_exct))))
    max_err2(l1)=max(max(max(abs(u2- u2_exct))))

    if(l1>1)
        cnv_rt1(l1-1)=log2(max_err1(l1-1)/max_err1(l1))
        cnv_rt2(l1-1)=log2(max_err2(l1-1)/max_err2(l1))
    end
end
% surf(x,y,u(:,:,Nt+1)')
