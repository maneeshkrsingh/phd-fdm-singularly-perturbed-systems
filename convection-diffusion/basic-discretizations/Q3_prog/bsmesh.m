function [x,h]=bsmesh(sig_o1,sig_o2,sig1,sig2,dis,N,e)
% sig_o=1;
% sig1=min(dis/2,sig_o*e*log(N));
% sig2=min((1-dis)/2,sig_o*e*log(N));
x=zeros(1,N+1);
h=zeros(1,N);

% for i=1:(T/del_t)+1
%     t(i)=(i-1)*del_t;
% end
for i=1:N/4
%      x(i)=sig_o*e*(-log(1-4*(1-1/N)*((i-1)/N)));
      x(i)=(4*(dis-sig1)*(i-1)/N);
end

for i=N/4+1:N/2+1

 %     x(i)=sig1+(4*(dis-sig1)/N)*((i-1)-N/4);
      x(i)=dis+sig_o1*e*(log(1-4*(1-1/N)*((1/2)-(i-1)/N)));
end

for i=N/2+2:3*N/4
%     x(i)=dis+sig_o*e*(-log(1-4*(1-1/N)*(((i-1)/N)-(1/2))));
      x(i)=dis+sig_o2*e*(-log(1-4*(1-1/N)*((i-1)/N-(1/2))));
end

for i=3*N/4+1:N+1
     x(i)=dis+sig2+(4*(1-dis-sig2)/N)*((i-1)-(3*N/4));
end

if(sig1==dis/2 && sig2==(1-dis)/2)
for i=1:N+1
    x(i)=(i-1)/N;
end

end
% x(1)=0;
% for i=1:N
%    x(i+1)=x(i)+1/N;
% end

% x(1)=0;
% 
% x(N+1)=1;