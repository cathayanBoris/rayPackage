[xtopo,ytopo,ztopoS] = getMeSmoothed(30000,-91.5,-87.5,25.8,27.8);
hg = -ztopoS;
[hxg,hyg,gradh,dx,dy] = getMeGradZ(xtopo,ytopo,hg);
[hxxg,hxyg,~,~,~] = getMeGradZ(xtopo,ytopo,hxg);
[~,hyyg,~,~,~] = getMeGradZ(xtopo,ytopo,hyg);
originLon = nanmean(xtopo);
originLat = nanmean(ytopo);
meandx = nanmean(dx); % NOTE: dx is not uniform
meandy = nanmean(dy);
xg = (xtopo - originLon) * meandx * 60;
yg = (ytopo - originLat) * meandy * 60;

%%
% figure(1)
% clf
% contourf(xg,yg,ztopoS',[-4000:100:0],'ShowText','on')

figure(2)
clf
contourf(xg,yg,hg',[0:100:4000],'ShowText','on','LabelSpacing',500);
hold on
aa = contour(xg,yg,hg',[0:50:1650],'m','ShowText','off','LineWidth',3);
bb = contour(xg(140:end),yg,hxg(140:end,:)',[0.005:0.001:1],'r','ShowText','off','LineWidth',3);
cc = contour(xg(140:end),yg,hyg(140:end,:)',[-0.001:0.001:0.1],'b','ShowText','off','LineWidth',3);
dd = contour(xg,yg,gradh',[0.04:0.001:0.1],'k','ShowText','off','LineWidth',3);

contour(xg,yg,hg',[0:100:4000],'ShowText','on','LabelSpacing',500);
colormap(flipud(colormap(parula)))
% legend(aa,"Depth Restriction")
% legend(bb,"hx Restriction")
% legend(cc,"hy Restriction")
% legend(dd,"∇h restriction")