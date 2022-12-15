function colorCodes = myColorCodes(n)

%Output: colorCodes (1 x 3 RGB vector)
%
%Usage: 
% 1) h = line(x,y,'Color',myColorCodes(7)); 
% 2) myColorCodes('demo') %Creates a figure displaying  all colors
%
%See also: PLOT, LINE, AXES

%Author: Korhan Gunel, Ph.D., mathematicians
%Adnan Menderes University, Dept. of Mathematics
%July 2017; 

if nargin == 0
   error('Must provide an input argument to myColorCodes')
end

colorOrder = ...
[ 0 0 1 % 1 BLUE
   0 1 0 % 2 GREEN (pale)
   1 0 0 % 3 RED
   0 1 1 % 4 CYAN
   1 0 1 % 5 MAGENTA (pale)
   1 1 0 % 6 YELLOW (pale)
   0 0 0 % 7 BLACK
   0 0.75 0.75 % 8 TURQUOISE
   0 0.5 0 % 9 GREEN (dark)
   0.75 0.75 0 % 10 YELLOW (dark)
   1 0.50 0.25 % 11 ORANGE
   0.75 0 0.75 % 12 MAGENTA (dark)
   0.7 0.7 0.7 % 13 GREY
   0.8 0.7 0.6 % 14 BROWN (pale)
   0.6 0.5 0.4 % 15 BROWN (dark)
]; 

% Algolia Color Codes 
algolia = ['#050f2c';
'#003666';
'#00aeff'; 
'#3369e7'; 
'#8e43e7'; 
'#b84592';
'#ff4f81';
'#ff6c5f';
'#ffc168'; 
'#2dde98';
'#1cc7d0';
'#cecece';
'#f0f0f0'];
newColors = colorCodeConverterHex2Double(algolia);
colorOrder = [colorOrder ; newColors];

if isnumeric(n) & n >= 1 & n <= length(colorOrder)
    colorCodes = colorOrder(n,:);
elseif strcmp(n,'demo')
    %GENERATE PLOT to display a sample of the line colors
    figure, axes;
    %PLOT N horizontal lines
    for n=1:length(colorOrder)
        h(n) = line([0 1],[n n],'Color',colorOrder(n,:));
    end
    set(h,'LineWidth',5)
    set(gca,'YLim',[0 n+1],'YTick',[1:n],'XTick',[])
    ylabel('Color Number');
else
    error('Invalid input to myColorCodes');
end