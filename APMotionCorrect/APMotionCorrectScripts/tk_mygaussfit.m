function [sigma,mu,A]=tk_mygaussfit(x,y,frame,ratio)

%
% [sigma,mu,A]=mygaussfit(x,y)
% [sigma,mu,A]=mygaussfit(x,y,h)
%
% this function is doing fit to the function
% y=A * exp( -(x-mu)^2 / (2*sigma^2) )
%
% the fitting is been done by a polyfit
% the lan of the data.
%
% frame decides how many frames to consider.
% frame should be a positive integer.
% if frame have not been taken it is set to be 5
% as default.
%
% ratio is the threshold which is the fraction
% from the maximum y height that the data
% is been taken from.
% ratio should be a number between 0-1.
% if ratio have not been taken it is set to be 0.5
% as default.



%% threshold
if nargin<=2, frame=5; end
if nargin<=3, ratio=0.5; end

%% cutting
y=reshape(y,1,numel(y));%make sure the dimension is right
ymax=max(y);
yindex = find(y>=ymax*ratio);
yleft = max(2,yindex(1)-frame);
yright = min(size(y,2)-1,yindex(end)+frame);
xnewbase = x(yleft-1:yright+1);
ynewbase = [ymax/100,y(yleft:yright),ymax/100];
ynew = ynewbase(ynewbase~=0);
xnew = xnewbase(ynewbase~=0);
% ynew=ones(1,size(y,2))*ymax/1000;
% ynew(yleft:yright)=y(yleft:yright);
% xnew=x;

%% fitting
ylog=log(ynew);
xlog=xnew;
p=polyfit(xlog,ylog,2);
A2=p(1);
A1=p(2);
A0=p(3);
sigma=sqrt(-1/(2*A2));
mu=A1*sigma^2;
A=exp(A0+mu^2/(2*sigma^2));

if ~isreal(sigma)
    error('sigma is not a real number')
end
if mu<x(yleft) || mu>x(yright)
    error('mu is out of range')
end

