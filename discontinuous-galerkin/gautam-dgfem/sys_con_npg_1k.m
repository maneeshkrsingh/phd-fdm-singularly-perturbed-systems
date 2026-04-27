%%%%%%%%new convection problem 
clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input valuess by user.......................
t=input('enter the value of epsilon in btwn 0 and 1:');
mu=input('enter the value of mu in btwn 0 and 1 and less than eps:');
n=input('enter the no of interval:');
ss=1;
penal=1;
%%%%%some defined values
% dimension of local matrices
locdim = 2;
 N=zeros(5,1);
 ord=zeros(4,1);
 err=zeros(5,1);
 w=zeros(5,1);
     for p=1:5
         nel=(2^(p-1))*n;
        N(p)=nel;
% dimension of global matrix
glodim = nel * locdim;
% initialize to zero matrix and right-hand side vector
Aglobal = zeros(glodim,glodim);
Dglobal =zeros(glodim,glodim);
Bglobal=zeros(glodim,glodim);
Cglobal=zeros(glodim,glodim);
rhsglobal = zeros(glodim,1);
rhsgloball= zeros(glodim,1);
% Gauss quadrature weights and points
wg(1) = 1.0;
wg(2) = 1.0;
sg(1) = -0.577350269189;
sg(2) = 0.577350269189;

%%%Shishkin mesh.......
sigma1=min(1/2,2*mu*log(nel));
sigma2=min(sigma1/2,2*t*log(nel));
x(1)=0;
 for i=1:(nel/4)
      h(i)=(4*(sigma2))/nel;
     for j=1:(locdim-1)
     x(locdim*(i-1)+j+1)=x(locdim*(i-1)+j)+h(i);
     end
     if(i~=nel)
    x(locdim*(i-1)+locdim+1)=x(locdim*(i-1)+locdim); 
 end
 end
 
 for i=(nel/4)+1:(nel/2)
       h(i)=(4*(sigma1-sigma2))/nel;
     for j=1:(locdim-1)
     x(locdim*(i-1)+j+1)=x(locdim*(i-1)+j)+h(i);
     end
     if(i~=nel)
    x(locdim*(i-1)+locdim+1)=x(locdim*(i-1)+locdim); 
 end
 end
 
 for i=(nel/2)+1:(nel)
       h(i)=(2*(1-sigma1))/nel;
     for j=1:(locdim-1)
     x(locdim*(i-1)+j+1)=x(locdim*(i-1)+j)+h(i);
     end
     if(i~=nel)
    x(locdim*(i-1)+locdim+1)=x(locdim*(i-1)+locdim); 
 end
 end
 
 
 
 
 for i=1:nel
     y(i)=x(2*(i-1)+1);
 end
 y(nel+1)=x(glodim);
 
 
 c=2;
 b=1;
 l=1;

Amat(1,1,l) =c*h(l);
Amat(1,2,l)=-2*b;
Amat(2,1,l)=0;
Amat(2,2,l)=c*(h(l)/3)+((4*t)/h(l));

Bmat(1,1,l)=nel*penal;
Bmat(1,2,l)=(t/h(l))-((nel)*penal);
Bmat(2,1,l)=((-t*ss)/h(l))-(nel*penal);
Bmat(2,2,l)=((-t+t*ss)/h(l))+(nel*penal);


for l=2:(nel/2)
  
Amat(1,1,l) =c*h(l);
Amat(1,2,l)=-2*b;
Amat(2,1,l)=0;
Amat(2,2,l)=c*(h(l)/3)+((4*t)/h(l));

Bmat(1,1,l)=nel*penal;
Bmat(1,2,l)=(t/h(l))-(nel)*penal;
Bmat(2,1,l)=((-t*ss)/h(l))-(nel*penal);
Bmat(2,2,l)=((-t+t*ss)/h(l))+(nel*penal);

Cmat(1,1,l)=b+(nel*penal);
Cmat(1,2,l)=b+(-t/h(l-1))+(nel*penal);
Cmat(2,1,l)=b+((t*ss)/h(l-1))+(nel*penal);
Cmat(2,2,l)=b+((-t+t*ss)/h(l-1))+(nel*penal);

Dmat(1,1,l)=-b-(nel*penal);
Dmat(1,2,l)=b+(-t/h(l))+(nel*penal);
Dmat(2,1,l)=-b+((-t*ss)/h(l-1))-(nel*penal);
Dmat(2,2,l)=b+(-t/h(l))+((t*ss)/h(l-1))+(nel*penal);

Emat(1,1,l)=(nel)*(-penal);
Emat(1,2,l)=(t/h(l-1))-(nel*penal);
Emat(2,1,l)=(t*ss/h(l))+(nel*penal);
Emat(2,2,l)=(-t/h(l-1))+(t*ss/h(l))+(nel*penal);
                
end
 

for l=(nel/2)+1:nel
    
Amat(1,1,l) =c*h(l);
Amat(1,2,l)=-2*b;
Amat(2,1,l)=0;
Amat(2,2,l)=c*(h(l)/3)+((4*t)/h(l));

Bmat(1,1,l)=1*penal;
Bmat(1,2,l)=(t/h(l))-(1)*penal;
Bmat(2,1,l)=((-t*ss)/h(l))-(1*penal);
Bmat(2,2,l)=((-t+t*ss)/h(l))+(1*penal);

Cmat(1,1,l)=b+(1*penal);
Cmat(1,2,l)=b+(-t/h(l-1))+(1*penal);
Cmat(2,1,l)=b+((t*ss)/h(l-1))+(1*penal);
Cmat(2,2,l)=b+((-t+t*ss)/h(l-1))+(1*penal);

Dmat(1,1,l)=-b-(1*penal);
Dmat(1,2,l)=b+(-t/h(l))+(1*penal);
Dmat(2,1,l)=-b+((-t*ss)/h(l-1))-(1*penal);
Dmat(2,2,l)=b+(-t/h(l))+((t*ss)/h(l-1))+(1*penal);

Emat(1,1,l)=(1)*(-penal);
Emat(1,2,l)=(t/h(l-1))-(1*penal);
Emat(2,1,l)=(t*ss/h(l))+(1*penal);
Emat(2,2,l)=(-t/h(l-1))+(t*ss/h(l))+(1*penal);  
end

F0mat =[-b+(nel)*penal b+((2*t)/h(1))-(nel)*penal; b-((2*t*ss)/h(1))-(nel)*penal -b-((2*t)/h(1))+((2*t*ss)/h(1))+(nel)*penal ];
FNmat =[1*penal (-2*t/h(nel))+1*penal; ((2*t*ss)/h(nel))+1*penal ((-2*t+2*t*ss)/h(nel))+1*penal ];


% assemble global matrix and right-hand side......................
% first block row


for ii=1:locdim
for jj=1:locdim
Aglobal(ii,jj) = Aglobal(ii,jj)+Amat(ii,jj,1)+F0mat(ii,jj)+Cmat(ii,jj,2);
je = locdim+jj;
Aglobal(ii,je) =Aglobal(ii,je)+Dmat(ii,jj,2);
end; %jj
end; %ii
% compute right-hand side.........
for ig=1:2
rhsglobal(1)=rhsglobal(1)+wg(ig)*(h(1)/2)*sourceff(((sg(ig)+1)*(h(1)/2)),t,mu);
rhsglobal(2)=rhsglobal(2)+wg(ig)*(h(1)/2)*sg(ig)*sourceff(((sg(ig)+1)*(h(1)/2)),t,mu);
end; %ig

% intermediate block rows
% loop over elements
for i=2:(nel-1)
for ii=1:locdim
ie = ii+(i-1)*locdim;
for jj=1:locdim
je = jj+(i-1)*locdim;
Aglobal(ie,je) =Aglobal(ie,je)+Amat(ii,jj,i)+Bmat(ii,jj,i)+Cmat(ii,jj,i+1);
je = jj+(i-2)*locdim;
Aglobal(ie,je) = Aglobal(ie,je)+Emat(ii,jj,i);
je = jj+(i)*locdim;
Aglobal(ie,je) =Aglobal(ie,je)+Dmat(ii,jj,i+1);
end; %jj

% compute right-hand side
for ig=1:2
rhsglobal(ie)=rhsglobal(ie)+wg(ig)*(sg(ig)^(ii-1))*(h(i)/2)*sourceff(((h(i)/2)*sg(ig))+y(i)+(h(i)/2),t,mu);
end; %ig
end; %ii
end; %i

% last block row
for ii=1:locdim
ie = ii+(nel-1)*locdim;
for jj=1:locdim
je = jj+(nel-1)*locdim;
Aglobal(ie,je) = Aglobal(ie,je)+Amat(ii,jj,nel)+FNmat(ii,jj)+Bmat(ii,jj,nel);
je = jj+(nel-2)*locdim;
Aglobal(ie,je) = Aglobal(ie,je)+Emat(ii,jj,nel);
end; %jj
% compute right-hand side
for ig=1:2
rhsglobal(ie)=rhsglobal(ie)+wg(ig)*(sg(ig)^(ii-1))*(h(nel)/2)*sourceff(((h(nel)/2)*sg(ig))+y(nel)+(h(nel)/2),t,mu);
end; %ig
end; %ii

%%%%%only for system......
 c1=4;
 b1=2;
l=1;

Amatt(1,1,l) =c1*h(l);
Amatt(1,2,l)=-2*b1;
Amatt(2,1,l)=0;
Amatt(2,2,l)=c1*(h(l)/3)+((4*mu)/h(l));

Bmatt(1,1,l)=nel*penal;
Bmatt(1,2,l)=(mu/h(l))-((nel)*penal);
Bmatt(2,1,l)=((-mu*ss)/h(l))-(nel*penal);
Bmatt(2,2,l)=((-mu+mu*ss)/h(l))+(nel*penal);


for l=2:(nel/2)
  
Amatt(1,1,l) =c1*h(l);
Amatt(1,2,l)=-2*b1;
Amatt(2,1,l)=0;
Amatt(2,2,l)=c1*(h(l)/3)+((4*mu)/h(l));

Bmatt(1,1,l)=nel*penal;
Bmatt(1,2,l)=(mu/h(l))-(nel)*penal;
Bmatt(2,1,l)=((-mu*ss)/h(l))-(nel*penal);
Bmatt(2,2,l)=((-mu+mu*ss)/h(l))+(nel*penal);

Cmatt(1,1,l)=b1+(nel*penal);
Cmatt(1,2,l)=b1+(-mu/h(l-1))+(nel*penal);
Cmatt(2,1,l)=b1+((mu*ss)/h(l-1))+(nel*penal);
Cmatt(2,2,l)=b1+((-mu+mu*ss)/h(l-1))+(nel*penal);

Dmatt(1,1,l)=-b1-(nel*penal);
Dmatt(1,2,l)=b1+(-mu/h(l))+(nel*penal);
Dmatt(2,1,l)=-b1+((-mu*ss)/h(l-1))-(nel*penal);
Dmatt(2,2,l)=b1+(-mu/h(l))+((mu*ss)/h(l-1))+(nel*penal);

Ematt(1,1,l)=(nel)*(-penal);
Ematt(1,2,l)=(mu/h(l-1))-(nel*penal);
Ematt(2,1,l)=(mu*ss/h(l))+(nel*penal);
Ematt(2,2,l)=(-mu/h(l-1))+(mu*ss/h(l))+(nel*penal);
                
end
 

for l=(nel/2)+1:nel
    
Amatt(1,1,l) =c1*h(l);
Amatt(1,2,l)=-2*b1;
Amatt(2,1,l)=0;
Amatt(2,2,l)=c1*(h(l)/3)+((4*mu)/h(l));

Bmatt(1,1,l)=1*penal;
Bmatt(1,2,l)=(mu/h(l))-(1)*penal;
Bmatt(2,1,l)=((-mu*ss)/h(l))-(1*penal);
Bmatt(2,2,l)=((-mu+mu*ss)/h(l))+(1*penal);

Cmatt(1,1,l)=b1+(1*penal);
Cmatt(1,2,l)=b1+(-mu/h(l-1))+(1*penal);
Cmatt(2,1,l)=b1+((mu*ss)/h(l-1))+(1*penal);
Cmatt(2,2,l)=b1+((-mu+mu*ss)/h(l-1))+(1*penal);

Dmatt(1,1,l)=-b1-(1*penal);
Dmatt(1,2,l)=b1+(-mu/h(l))+(1*penal);
Dmatt(2,1,l)=-b1+((-mu*ss)/h(l-1))-(1*penal);
Dmatt(2,2,l)=b1+(-mu/h(l))+((mu*ss)/h(l-1))+(1*penal);

Ematt(1,1,l)=(1)*(-penal);
Ematt(1,2,l)=(mu/h(l-1))-(1*penal);
Ematt(2,1,l)=(mu*ss/h(l))+(1*penal);
Ematt(2,2,l)=(-mu/h(l-1))+(mu*ss/h(l))+(1*penal);  
end





F0matt =[-b1+(nel)*penal b1+((2*mu)/h(1))-(nel)*penal; b1-((2*mu*ss)/h(1))-(nel)*penal -b1-((2*mu)/h(1))+((2*mu*ss)/h(1))+(nel)*penal ];
FNmatt =[1*penal (-2*mu/h(nel))+1*penal; ((2*mu*ss)/h(nel))+1*penal ((-2*mu+2*mu*ss)/h(nel))+1*penal ];
%%aseembling matrixxxxxx.....
%%first block.........................
for ii=1:locdim
for jj=1:locdim
Dglobal(ii,jj) = Dglobal(ii,jj)+Amatt(ii,jj,1)+F0matt(ii,jj)+Cmatt(ii,jj,2);
je = locdim+jj;
Dglobal(ii,je) =Dglobal(ii,je)+Dmatt(ii,jj,2);
end; %jj
end; %ii
% compute right-hand side.........
for ig=1:2
rhsgloball(1)=rhsgloball(1)+wg(ig)*(h(1)/2)*sourceff1(((sg(ig)+1)*(h(1)/2)),t,mu);
rhsgloball(2)=rhsgloball(2)+wg(ig)*(h(1)/2)*sg(ig)*sourceff1(((sg(ig)+1)*(h(1)/2)),t,mu);
end; %ig

% intermediate block rows
% loop over elements
for i=2:(nel-1)
for ii=1:locdim
ie = ii+(i-1)*locdim;
for jj=1:locdim
je = jj+(i-1)*locdim;
Dglobal(ie,je) =Dglobal(ie,je)+Amatt(ii,jj,i)+Bmatt(ii,jj,i)+Cmatt(ii,jj,i+1);
je = jj+(i-2)*locdim;
Dglobal(ie,je) = Dglobal(ie,je)+Ematt(ii,jj,i);
je = jj+(i)*locdim;
Dglobal(ie,je) =Dglobal(ie,je)+Dmatt(ii,jj,i+1);
end; %jj

% compute right-hand side
for ig=1:2
rhsgloball(ie)=rhsgloball(ie)+wg(ig)*(sg(ig)^(ii-1))*(h(i)/2)*sourceff1(((h(i)/2)*sg(ig))+y(i)+(h(i)/2),t,mu);
end; %ig
end; %ii
end; %i

% last block row
for ii=1:locdim
ie = ii+(nel-1)*locdim;
for jj=1:locdim
je = jj+(nel-1)*locdim;
Dglobal(ie,je) = Dglobal(ie,je)+Amatt(ii,jj,nel)+FNmatt(ii,jj)+Bmatt(ii,jj,nel);
je = jj+(nel-2)*locdim;
Dglobal(ie,je) = Dglobal(ie,je)+Ematt(ii,jj,nel);
end; %jj
% compute right-hand side
for ig=1:2
rhsgloball(ie)=rhsgloball(ie)+wg(ig)*(sg(ig)^(ii-1))*(h(nel)/2)*sourceff1(((h(nel)/2)*sg(ig))+y(nel)+(h(nel)/2),t,mu);
end; %ig
end; %ii

for m=1:nel
    Bglobal(2*m-1,2*m-1)=-h(m);
end
for m=1:nel
    Bglobal(2*m,2*m)=-h(m)/3;
end
Cglobal=Bglobal;
Lglobal=[Aglobal Bglobal; Cglobal Dglobal];
for m=1:glodim
    bglobal(m)=rhsglobal(m);
end
for m=1:glodim
    bglobal(glodim+m)=rhsgloball(m);
end


ysol = Lglobal\bglobal';



ysol1=zeros(glodim,1);
ysol2=zeros(glodim,1);

for m=1:glodim
    ysol1(m)=ysol(m);
end
for m=1:glodim
    ysol2(m)=ysol(glodim+m);
end

base=zeros(1,glodim);
base1=zeros(1,glodim);
base2=zeros(1,glodim);
base3=zeros(1,glodim);
base4=zeros(1,glodim);
sol1=zeros(glodim,1);
sol2=zeros(glodim,1);
sol3=zeros(nel+1,1);
sol4=zeros(nel+1,1);
sol5=zeros(nel+1,1);
sol6=zeros(nel+1,1);
%%%numerical solutionnn........u_1
for i=1:nel
    base=zeros(1,glodim);
    for k=(locdim*(i-1)+1):locdim*i
        base=zeros(1,glodim);
        for j=1:locdim
           base(locdim*(i-1)+j)=(-1)^(k*(j-1));
        end
     sol1(k)=base*ysol1; 
    end       
end

%%%Numerical solution.....u2
for i=1:nel
    base2=zeros(1,glodim);
    for k=(locdim*(i-1)+1):locdim*i
        base2=zeros(1,glodim);
        for j=1:locdim
           base2(locdim*(i-1)+j)=(-1)^(k*(j-1));
        end
     sol2(k)=base2*ysol2; 
    end       
end
v=0;
w=0;
vv=0;
ww=0;
%%%Exact Solution............
for m=1:glodim
 vv(m)=(((1-exp(-x(m)/t))/(1*(1-exp(-1/t))))+((1-exp(-x(m)/mu))/(1*(1-exp(-1/mu))))-2*sin((pi/2)*x(m)));
end
for m=1:glodim
  ww(m)=((1-exp(-x(m)/mu))/(1-exp(-1/mu)))-(x(m)*exp(x(m)-1));
end

%%%%exact solution...............
for m=1:(nel+1)
 v(m)=(((1-exp(-y(m)/t))/(1*(1-exp(-1/t))))+((1-exp(-y(m)/mu))/(1*(1-exp(-1/mu))))-2*sin((pi/2)*y(m)));
end
for m=1:(nel+1)
  w(m)=((1-exp(-y(m)/mu))/(1-exp(-1/mu)))-(y(m)*exp(y(m)-1));
end

%%%%Derivative of the exact solution....
%  for m=1:(nel+1)
%   derv(m)=(((exp(-x(m)/t))/(t*(1-exp(-1/t))))+((exp(-x(m)/mu))/(mu*(1-exp(-1/mu))))-(pi)*cos((pi/2)*x(m)));  
%  end
%  for m=1:(nel+1)
%     derw(m)=(exp(-x(m)/mu)/((mu)*(1-exp(-1/mu)))-(exp(x(m)-1))-(x(m)*exp(x(m)-1)));
%  end
%  
for m=1:(nel+1)
    derv(m)=(exp(-y(m)/t)/(t*(1-exp(-1/t))))+(exp(-y(m)/mu)/(mu*(1-exp(-1/mu))))-(pi*cos((pi/2)*y(m)));
end
for m=1:(nel+1)
    derw(m)=(exp(-y(m)/mu)/(mu*(1-exp(-1/mu))))-((y(m)+1)*exp(y(m)-1));
end

for m=1:glodim
    dervv(m)=(exp(-x(m)/t)/(t*(1-exp(-1/t))))+(exp(-x(m)/mu)/(mu*(1-exp(-1/mu))))-(pi*cos((pi/2)*x(m)));
end
for m=1:glodim
    derww(m)=(exp(-x(m)/mu)/(mu*(1-exp(-1/mu))))-((x(m)+1)*exp(x(m)-1));
end

%%%%derivative of the numerical solution.......
for i=1:nel
    base3=zeros(1,glodim);
    for k=(locdim*(i-1)+1):locdim*i
        base3=zeros(1,glodim);
        for j=1:locdim
             %base2(locdim*(i-1)+j)=(2*nel)*(j-1)*(((2*nel)*(x(k)-((i-1)+(1/2))*(1/nel)))^(j-2));
             %base2(locdim*(i-1)+j)=(2/h(i))*(j-1)*(((2/h(i))*(x(k)-(y(i)+(h(i)/2))))^(j-2));
             base3(locdim*(i-1)+j)=(2/h(i))*(j-1);
        end
     sol3(k)=base3*ysol1; 
    end      
end

for i=1:nel
    base4=zeros(1,glodim);
    for k=(locdim*(i-1)+1):locdim*i
        base4=zeros(1,glodim);
        for j=1:locdim
            % base2(locdim*(i-1)+j)=(2*nel)*(j-1)*(((2*nel)*(x(k)-((i-1)+(1/2))*(1/nel)))^(j-2));
            % base2(locdim*(i-1)+j)=(2/h(i))*(j-1)*(((2/h(i))*(x(k)-(y(i)+(h(i)/2))))^(j-2));
             base4(locdim*(i-1)+j)=(2/h(i))*(j-1);
        end
     sol4(k)=base4*ysol2; 
    end      
end
%%%%%%L2 and Energy norm Error.........
% for m=1:nel
%     sol3(m)=sol1(2*(m-1)+1);
% end
% sol3(nel+1)=sol1(glodim);
% for m=1:nel
%     sol4(m)=sol2(2*(m-1)+1);
% end
% sol4(nel+1)=sol2(glodim);
% for m=1:nel
%     sol5(m)=((sol3(m+1)-sol3(m)))/h(m);
% end
% sol5(nel+1)=sol5(nel);
% for m=1:nel
%     sol6(m)=((sol4(m+1)-sol4(m)))/h(m);
% end
% sol6(nel+1)=sol6(nel);
%%%%%fdm error......

%  for m=1:(nel+1)
%     f(m)=(abs(v(m)-sol3(m)))^2;
%  end
%  for m=1:(nel+1)
%      g(m)=(abs(w(m)-sol4(m)))^2;
%  end 
%  for m=1:(nel+1)
%      df(m)=(abs(derv(m)-sol5(m)))^2;
%  end
%  
%  for m=1:(nel+1)
%      dg(m)=(abs(derw(m)-sol6(m)))^2;
%  end
 
 
 for m=1:glodim
    f(m)=(abs(vv(m)-sol1(m)))^2;
 end
 for m=1:glodim
     g(m)=(abs(ww(m)-sol2(m)))^2;
 end 
 for m=1:glodim
     df(m)=(abs(dervv(m)-sol3(m)))^2;
 end
 
 for m=1:glodim
     dg(m)=(abs(derww(m)-sol4(m)))^2;
 end
 %%penalty terms...
 
 for m=1:glodim
    pnl(m)=(vv(m)-sol1(m));
 end
for m=1:glodim
    pnll(m)=(ww(m)-sol2(m));
end
 penalty=0;
 penalty=((-1/2)+nel)*(pnl(1)^2);
 for m=2:2:(glodim/2)
 penalty=penalty+((-1/2)+nel)*(pnl(m+1)-pnl(m))^2;
 end
 for m=(glodim/2)+1:2:(glodim-1)
     penalty=penalty+(1/2)*(pnl(m+1)-pnl(m))^2;
 end
 penalty=penalty+(1/2)*(pnl(glodim))^2;
 
 penalty=penalty+(-1+nel)*(pnll(1)^2);
 for m=2:2:(glodim/2)
 penalty=penalty+(-1+nel)*(pnll(m+1)-pnll(m))^2;
 end
 for m=(glodim/2)+1:2:(glodim-1)
     penalty=penalty+(0)*(pnl(m+1)-pnl(m))^2;
 end
 penalty=penalty+(0)*(pnl(glodim))^2;
 err(p)=lnormspp3(y,f,g,df,dg,h,t,mu,penalty);
 
if(p>1)
    ord(p)=log2(err(p-1)/err(p));
end
 end

 err'
 ord'
% %%%%Penalty terms....
%  penalty=0;
%  penalty=(1/2)*pnl(1)^2;
%  for m=2:(glodim/2)
%  penalty=penalty+(1/2)*(pnl(m+1)-pnl(m))^2;
%  end
%  for m=(glodim/2)+1:(glodim-1)
%      penalty=penalty+((-1/2)+nel)*(pnl(m+1)-pnl(m))^2;
%  end
%  penalty=penalty+((-1/2)+nel)*(pnl(m))^2;

%%%%graphhh..................
 plot(vv,'r')
 hold on
 plot(sol1,'b')
 hold off