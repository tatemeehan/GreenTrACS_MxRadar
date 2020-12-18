function varargout=shadedErrorBarT8(x,y,errBar,tau,lineProps,transparent)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:    
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. 
%
% Outputs
% H - a structure of handles to the generated plot objects.     
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});    
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1); 
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1); 
% hold off
%
%
% Rob Campbell - November 2009


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Error checking    
narginchk(3,5)


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar) 
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:).';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:).';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<5, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<6, transparent=0; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Plot to get the parameters of the line 
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.85;
patchSaturation=0.15; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end

    
%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];
% Y indicies
yT8 = [y];
% X indicies
xT8 = [x];
% Color indicies
cT8 = [ones(1,length(x))];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];

% Number of Faces
N = length(x);
% Array of Verticies
q = (1:N-1)';
% Lower Error Verticies
vert1 = [x(:),lE(:);x(:),y(:)];
% Upper Error Verticies
vert2 = [x(:),uE(:);x(:),y(:)];
% Array of Faces to Paint
face = [q, q+1, q+N+1, q+N];
% Array of Face Colors
cT8 = [zeros(length(x),1);ones(length(x),1)];
% Create Gaussian Distribution for Colormap
sigma = mean(errBar(1,:)./tau);
gauss = normpdf(linspace(mean(y),mean(y)+.9945*sigma,256),mean(y),sigma);

% Create Linear ColorMap1
color1 = linspace(1,col(1),256);
if any(color1~=color1(1));
% Map Linear Array to Gaussian Distribution and Normalize Color Range
range = max(color1.*gauss) - min(color1.*gauss);
color1G = (color1.*gauss - min(color1.*gauss))./range;
range2 = color1(1) - color1(end); 
color1G = (color1G.*range2) + color1(end);
else
    color1G = color1;
end
% Create Linear ColorMap2
color2 = linspace(1,col(2),256);
if any(color2~=color2(1));
% Map Linear Array to Gaussian Distribution and Normalize Color Range
range = max(color2.*gauss) - min(color2.*gauss);
color2G = (color2.*gauss - min(color2.*gauss))./range;
range2 = color2(1) - color2(end); 
color2G = (color2G.*range2) + color2(end);
else
    color2G = color2;
end
% Create Linear ColorMap3
color3 = linspace(1,col(3),256);
if any(color3~=color3(1));
% Map Linear Array to Gaussian Distribution and Normalize Color Range
range = max(color3.*gauss) - min(color3.*gauss);
color3G = (color3.*gauss - min(color3.*gauss))./range;
range2 = color3(1) - color3(end); 
color3G = (color3G.*range2) + color3(end);
else
    color3G = color3;
end
% Concatenate ColorMap Array
cT8map = [color1G;color2G;color3G]';
% Brighten Values - Default no
bright = 0.00;
if bright>0
    cT8map = cT8map.^(1-bright);
elseif bright<=0
    cT8map = cT8map.^(1/(1+bright));
end

% Paint Color Patch
H.patch = patch('Faces', face, 'Vertices', vert1, 'FaceVertexCData', cT8,...
    'FaceColor', 'interp', 'EdgeColor', 'none'); hold on;
patch('Faces', face, 'Vertices', vert2, 'FaceVertexCData', cT8,...
    'FaceColor', 'interp', 'EdgeColor', 'none'); colormap(cT8map);alpha(.85)


%Make pretty edges around the patch. 
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
uistack(H.mainLine,'top')


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end
