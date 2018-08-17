function [ output_args ] = irlsDix( Vrms )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% d = Vrms.^2.*(1:length(Vrms))';
% G = tril(ones(length(Vrms)));

% % A Constrained Interval Velocity Inversion
% G = tril(ones(length(mean(VrmsHVA,2))));
% d = mean(VrmsHVA,2).^2.*(1:length(mean(VrmsHVA,2)))'
% x = G\d
% while ~isempty(find(x<0))
% % Constrains Bad Horizons
% goodIx = find(x>0);
% G = tril(ones(length(goodIx)));
% d = mean(VrmsHVA(goodIx,:),2).^2.*(1:length(mean(VrmsHVA(goodIx,:),2)))'
% x = G\d
% end
% % Vx = sqrt(x)
% VirlsIntM = irls(G,d,.005000,.0005000,1,10);
% VirlsInt = sqrt(VirlsIntM);
rmvHorizon = 0;
tolii = 1;
% while ~isempty(rmvHorizon);
Zsurf = mean(ZlmoHVA(2),2)
Z = mean(ZnmoHVA,2)
H = [Zsurf;diff([Zsurf;Z])];
Zhl = [Zsurf;Z];
% clear rmvHorizon
rmvHorizon = [];
HorizonIx = 1:length(H);
tol = 1:.5:5;


% VirlsIntM = irls(G,mean(VrmsHVA,2),5000,5000,1,10);
% VirlsIntM = irls(G,d,5000,5000,1,10);
%%%%%
% Constrained RMS Gradient Inversion
do = [mean(VlmoHVA(2,:));mean(VrmsHVA,2)];
d = [mean(VlmoHVA(2,:));mean(VrmsHVA,2)]
% impose Herron Langway
phl = herronLangway( Zhl, -20, DryCrim(do(1)), 0.5);
rmsPhl = cumsum(phl)./(1:length(phl))';
Vhl = DryCrimVRMS(rmsPhl);
Bad = find(abs(Vhl-do)>0.005);
dfix = do; dfix(Bad) = []; dfix = [dfix;Vhl(end)];

G = tril(ones(length(d)));
% L2 RMS Velocity gradient
x = G\d;%x(1) = [];
xHL = G\Vhl;
HLgradient = mean(xHL(2:end));
% HLgradVar = tol(tolii).*std(xHL(2:end));
HLgradVar = 1.*std(xHL(2:end));

%d(1) = [];
% A stable gradient is -0.0015+-0.001
rx = (x-(HLgradient)); rxIx = find(abs(rx)>HLgradVar);rxIx(1) = []; goodIx = [1;find(abs(rx)<=HLgradVar)];
notsogood = goodIx((find((diff(sort(goodIx))~=1),1)+1):end);
if ~isempty(notsogood)
    goodIx(goodIx == notsogood) = [];
    rxIx = sort([notsogood;rxIx]);
end
% prx = (abs(rx)).^(-1);
% Nudge Data
vgrad = rx(rxIx).*(H(rxIx))./15;%prx(rxIx); 
[~,sortIx] = sort([goodIx;rxIx]);
newd = [d(goodIx);d(rxIx)-vgrad];
newd = newd(sortIx)
% Loop Over Residual Gradients
iter = 1
ii = 1;
while ~isempty(rxIx)%iter<=35%
% recompute gradient
Gx = tril(ones(length(newd)));
newx = Gx\newd;
%recompute layer thickness
dH = -sign(vgrad).*((d(rxIx)-newd(rxIx))./(x(rxIx)-newx(rxIx)));
H = [H(goodIx);H(rxIx)+dH];
H = H(sortIx);
if H(goodIx(end)+1)<0.05
    hIx = HorizonIx(goodIx(end)+1);
%     keyboard
    HorizonIx(goodIx(end)+1) = [];
    H(goodIx(end)+1) = [];
    newx(goodIx(end)+1) = [];
    newd(goodIx(end)+1) = [];
    x(goodIx(end)+1) = [];
    rmvHorizon(ii) = hIx;%(goodIx(end)+1);
%     keyboard
    ii = ii+1;
end
    
% Use Move out Model
% offsetArray

% A stable gradient is -0.0015+-0.001
rx = (newx-(HLgradient));rxIx = find(abs(rx)>HLgradVar & (~ismember(1:length(rx),goodIx))');
gooditr = find(abs(rx)<=HLgradVar);
goodIx = [goodIx;gooditr(~ismember(gooditr,goodIx))];
% Solve in order
notsogood = goodIx((find((diff(sort(goodIx))~=1),1)+1):end);
if ~isempty(notsogood)
    goodIx(find(ismember(goodIx,notsogood))) = [];
    rxIx = sort([notsogood;rxIx]);
end
% newprx = (abs(newrx)).^(-1);
% Nudge Data
vgrad = rx(rxIx).*(H(rxIx))./15;%newprx(newrxIx);
[~,sortIx] = sort([goodIx;rxIx]);
d = newd(sortIx);
newd = [newd(goodIx);newd(rxIx)-vgrad]; newd = newd(sortIx);
% x = newx;
iter = iter+1;
end
iter
tolii = tolii+1;
% end
% Compute Interval Velocity
rmsV = newd(sort(goodIx));
rmsV = newd
% rmsV = do
% rmsV = dfix;
% rmsV = Vhl;
rmsd = rmsV.^2.*(1:length(rmsV))';
Gint = tril(ones(length(rmsd)));
mint = Gint\rmsd;
Vl2Int = sqrt(mint)
% Run Forward
errin = 0;randn(12,1)./1000
frms = Gint*((Vl2Int.^2+errin));
frms = sqrt(frms./(1:length(rmsV))')
e = rmsV-frms


%irls solution
VirlsIntM = irls(Gint,rmsd,.005000,.0005000,1,10);
VirlsInt = sqrt(VirlsIntM);

% Estimate Depth With New Velocities
activeHorizons = HorizonIx(2:end)-1;
% Paul Michaels Thoughts
% Gotta Loop Here Bummer Dude
HorizonTime = fullFoldPicks.ReflectionFBpick(:,activeHorizons);
for jj = 1:size(HorizonTime,2)
% for jj = activeHorizons-1
    for kk = 1:size(HorizonTime,1)
    TwTime(kk,jj) = [mean(HorizonTime{kk,jj}(1:100))]-mean(DirectTo{1});
    To(kk,jj) = sqrt(TwTime(kk,jj).^2-(offsetArray(kk).^2./rmsV(jj).^2));
    Zed(kk,jj) = (To(kk,jj).*rmsV(jj))./2;
    end
    
end
Zo = [Zsurf,mean(Zed)];


vmsV
% H
% sum(H)
% % % A stable gradient is -0.0015+-0.001
% newrx = (newx-(-0.0015)); newrxIx = find(abs(newrx)>0.001 & (~ismember(1:length(newrx),goodIx))');
% newgoodIx = [goodIx;find(~ismember(find(abs(newrx)<=0.001),goodIx))];
% newprx = (abs(newrx)).^(-1);
% % Nudge Data
% vgrad = newrx(newrxIx).*(newH(newrxIx))./3;%newprx(newrxIx); 
% newnewd = [newd(newgoodIx);newd(newrxIx)-vgrad];


prx = abs(rx./-.001).^(-1);
sumprx = sum(prx); wsum = length(prx)./sumprx; wprx = prx.*wsum; 
newGx = tril(repmat(wprx',11,1));intsum = sum((1:length(prx))');
intwsum = intsum./sumprx; intwprx = prx.*intsum;
% dint = d.^2.*(wprx.*(1:length(prx))');%wprx;%ones(length(d),1)
dint = d.^2.*intwprx;
newx = newGx\dint
R=diag(prx)
d = d.^2.*(1:length(prx))'
Gx = tril(ones(length(prx)));
m2=(Gx'*R*Gx)\(Gx'*R*d);
% G*x
% x(1) = [];
% G(:,length(d)) = [];
% G(length(d),:) = [];
w = abs((norm(x-.00125)./x).^(-1))
R=diag(w.^(-1))
m2=(G'*R*G)\(G'*R*d);
d = mean(VrmsHVA,2).^2.*(1:length(mean(VrmsHVA,2)))';
  % put the weighting factors into R
%   R=diag(w)+1000;
%   R=diag(1./w)-1000;

  % find the solution to the weighted problem
%   m2=(G'*R*G)\(G'*R*d);
  m2=(G'*G+R)\(G'*d);

  VintDixW = sqrt(abs(m2))
  




% x = G\mean(VrmsHVA,2)

% VirlsInt = G*VirlsIntM;


end

