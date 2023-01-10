function figset(hh,xlab,ylab,varargin)
% figset(hh,xlab,ylab,fntsz,color)
%
% preferred figure settings
% 
% - hh: current figure handle (e.g., gca)
% - xlab: (string) xlabel
% - ylab: (string) ylabel
% - fntsz: (numeric) font size (typically between 15 and 22)
% - col: (string) figure's background color (optional); default: none. Use
%   e.g. 'w' to obtain default appearance.
%
% GLC, 20 Nov 2007-2014

col='none'; fntsz=22;
if ~isempty(varargin) 
    fntsz=varargin{1}; 
    if nargin>4 col=varargin{2}; end
end
xlabel(xlab,'fontsize',fntsz,'fontweight','normal');
ylabel(ylab,'fontsize',fntsz,'fontweight','normal');
set(hh,'fontsize',fntsz,'color',col,'fontweight','normal');
box on;