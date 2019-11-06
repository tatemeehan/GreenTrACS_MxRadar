%% ShadedStairbar
% Works for monotonic function
x = [1:6]';
X = [x;x];%flipud(x)]
X = sort(X(:));
Xx = [X(2:end);flipud(X(2:end-2))];
y = rand(1,5);y = y(:);
% Stairs
X = X(2:end-1);
Y = [y,y]'; Y = Y(:);
% figure();stairs(x,y); hold on;
ci = 0.25.*abs(randn(5,1));
yl = y - ci;
yh = y+ci;
% stairs(x,yl);stairs(x,yh)
Yyh = [yh,yh]';Yyh = Yyh(:);
Yyl = [flipud(yl),flipud(yl)]';Yyl = Yyl(:);
Yy = [Yyh;Yyl];
figure();
patch(Xx,Yy,[0.5 0 0],'FaceAlpha',.5); hold on
plot(X,Y)