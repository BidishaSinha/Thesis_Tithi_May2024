function y=psdsimfunc6(AT,cytvis,f,sig,x)
% cytvis = 610600;
kap = 0.6; mu = 0;
s=size(x);len=max(size(x));
y=zeros(size(x));yy=zeros(size(x));
psd=zeros(len,2);
mini=10;maxi=1250;
for ww= 1:len;
   w= x(ww); 
   psd(ww,1)= x(ww);psd(ww,2)=0;
   p=10^-5*0.25*kap/cytvis;q=12.5*sig/cytvis;r=(f*0.25*10^8)/cytvis;
   pp=(mu/(cytvis*kap))*(0.09/16*pi); %pp:additional term from BiophysJ 108(12) 2794

     for i = mini:maxi;
          psd(ww,2)=psd(ww,2)+(1/((2*pi*w)^2+(p*i^3+(pp+q)*i+r/i)^2));
           
      end
      yy(ww)=(10^3)*(AT/(cytvis*pi)).*psd(ww,2);

%     fun=@(xx)(1./((2*pi*w).^2+(p.*xx.^3+q.*xx+ r./xx).^2));% added 2 pi in front of w...
%     yy(ww)=(10^3)*(AT/(cytvis*pi)).*(integral(fun,mini,maxi));
end
y=yy;