%% INITIALISATION

addpath ../../../SW/macrosight;

initscript;
rootdir = handlesdir.pathtodir;

a = dir(fullfile(rootdir,[handlesdir.data '_mat_TRACKS_*']));
a={a.name};
a(contains(a,'NOTHING')) = [];

whichtrack = 8;

whichTrackIdx = find(contains(a, sprintf('lab%d-', whichtrack)));
b = dir(fullfile(rootdir, a{whichTrackIdx}, '*.mat'));
b={b.name};
vect = zeros(length(b),1);
for ix=1:length(b)
    c = b{ix}(end-6:end-4);
    d = str2num(c);
    jx=1;
    while isempty(d)
        jx=jx+1;
        d=str2num(c(jx:end));
    end
    vect(ix) = d;
end
clear ix jx c d;

[frnumbers, bx] = sort(vect);
b = b(bx);

subidx = 420:470;
b(~ismember(frnumbers, subidx)) = [];
bx(~ismember(frnumbers, subidx)) = [];
frnumbers(~ismember(frnumbers, subidx)) = [];
%% figure 1

fr = getdatafromhandles(handles, filenames{90});
figure(11)
imagesc(fr.X);
%regs = regionprops(fr.clumphandles.nonOverlappingClumps>0, 'BoundingBox');
hold on
for ix=1:size(regs)
    rectangle('Position',regs(ix).BoundingBox, 'linewidth', 2, 'edgecolor', 'r');
end

%%
figure(12)
for ix=1:size(regs)
    subplot(2,2,ix)
    imagesc(fr.X)
    axis(bbox2axis(regs(ix).BoundingBox))
    axis square
    axis off
end

%% figure 2. EVOLUTION OF SHAPES
ixx = fix(linspace(1,length(bx),8));
Nixx = length(ixx);

f2=figure(2);
clf
bb = bbox2axis(T(ixx,:).BoundingBox);
bb = [min(bb(:,1)) max(bb(:,2)) min(bb(:,3)) max(bb(:,4))];
for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    figure(f2)
    subplot(1,Nixx,whichfr)
    
    imagesc(frameinfo.X)
    axis(bb);
    title(sprintf('frame=%d', frnumbers(ixx(whichfr))), 'fontsize', 16);
    axis square
    axis off
end

%% figure 3. Synthetic data

pcirc = [];
pdrop = [];
pbi = [];
ptri = [];
circ = zeros(7,2);
drop = zeros(8,2);
bi = zeros(10,2);
tri = zeros(13,2);
for ix=1:200
    [~,b,pts] = genshape('circ', 0);
    pcirc = [pcirc;pts];
    circ = circ+pts;
    if ix==1
        bcirc = b{1};
    else
        bcirc = bcirc+b{1};
    end
    [~,b,pts] = genshape('drop', 0);
    pdrop = [pdrop;pts];
    drop = drop+pts;
    if ix==1
        bdrop = b{1};
    else
        bdrop = bdrop + b{1};
    end
    [~,b,pts] = genshape('bi', 0);
    pbi =[pbi;pts];
    bi = bi+pts;
    if ix==1
        bbi = b{1};
    else
        bbi = bbi+b{1};
    end
    [~,b,pts] = genshape('tri', 0);
    ptri = [ptri;pts];
    tri = tri+pts;
    if ix==1
        btri = b{1};
    else
        btri = btri + b{1};
    end
end

circ = circ./200;
drop = drop./200;
bi = bi./200;
tri = tri ./200;

bcirc = (bcirc./200) - repmat([256 256], size(bcirc,1),1);
bdrop = bdrop./200- repmat([256 256], size(bdrop,1),1);
bbi=bbi./200 - repmat([256 256], size(bbi,1),1);
btri=btri./200 - repmat([256 256], size(btri,1),1);

%%
figure(3)
clf
subplot(141)
plotBoundariesAndPoints([],[],bcirc,'-m')
scattpoinstandmeans(pcirc, circ);
title('Circle - no edges', 'fontsize', 16)
subplot(142)
plotBoundariesAndPoints([],[],bdrop,'-m')
scattpoinstandmeans(pdrop, drop);
title('Drop - 1 edge', 'fontsize', 16)
subplot(143)
plotBoundariesAndPoints([],[],bbi,'-m')
scattpoinstandmeans(pbi, bi);
title('BiDrop - 2 edges', 'fontsize', 16)
subplot(144)
plotBoundariesAndPoints([],[],btri,'-m')
scattpoinstandmeans(ptri, tri);
title('TriDrop - 3 edges', 'fontsize', 16)

%%
figure(31)
clf
subplot(241)
plotBoundariesAndPoints([],[],bcirc,'-m')
scattpoinstandmeans(pcirc, circ);
title('Circle - no edges', 'fontsize', 16)
subplot(242)
plotBoundariesAndPoints([],[],bdrop,'-m')
scattpoinstandmeans(pdrop, drop);
title('Drop - 1 edge', 'fontsize', 16)
subplot(243)
plotBoundariesAndPoints([],[],bbi,'-m')
scattpoinstandmeans(pbi, bi);
title('BiDrop - 2 edges', 'fontsize', 16)
subplot(244)
plotBoundariesAndPoints([],[],btri,'-m')
scattpoinstandmeans(ptri, tri);
title('TriDrop - 3 edges', 'fontsize', 16)


subplot(245)
[ag] = computeMultiAnglegram({bcirc});
imagesc(ag');
title('Points along boundary', 'fontsize', 16)
subplot(246)
[ag] = computeMultiAnglegram({bdrop});
imagesc(ag');
title('Points along boundary', 'fontsize', 16)
subplot(247)
[ag] = computeMultiAnglegram({bbi});
imagesc(ag');
title('Points along boundary', 'fontsize', 16)
subplot(248)
[ag] = computeMultiAnglegram({btri});
imagesc(ag');
title('Points along boundary', 'fontsize', 16)

%%
figure(32)
clf
subplot(141)
[ag] = computeMultiAnglegram({bcirc});
ag(end)=200;
imagesc(ag'); colorbar;
xlabel('Points along boundary',  'fontsize', 16);
title('Circle - no edges', 'fontsize', 16)
subplot(142)
[ag] = computeMultiAnglegram({bdrop});
ag(end)=200;
imagesc(ag');  colorbar;
xlabel('Points along boundary',  'fontsize', 16);
title('Drop - 1 edge', 'fontsize', 16)
subplot(143)
[ag] = computeMultiAnglegram({bbi});
ag(end)=200;
imagesc(ag'); colorbar;
xlabel('Points along boundary',  'fontsize', 16);
title('BiDrop - 2 edges', 'fontsize', 16)
subplot(144)
[ag] = computeMultiAnglegram({btri});
ag(end)=200;
imagesc(ag'); colorbar;
xlabel('Points along boundary',  'fontsize', 16);
title('TriDrop - 3 edges', 'fontsize', 16)
%% figure 5. EVOLUTION
ixx = fix(linspace(1,length(bx),8));
Nixx = length(ixx);

f33=figure(33);
set(gcf, 'Position', get(0, 'ScreenSize'))
subplot(3,Nixx,((1:Nixx)+Nixx));
plotleftright(frnumbers, m1);

subplot(3,Nixx,((1:Nixx)+2*Nixx));
plotleftright(frnumbers, m2);

%
bb = bbox2axis(T(ixx,:).BoundingBox);
bb = [min(bb(:,1)) max(bb(:,2)) min(bb(:,3)) max(bb(:,4))];
for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    figure(33)
    subplot(3,Nixx,whichfr)
    
    plotBoundariesAndPoints(frameinfo.X, frameinfo.initboundy, ...
        frameinfo.outboundy,'m-');
    plotBoundariesAndPoints([],[],cornies{ixx(whichfr)}, '*y');
    
    axis(bb);
    title(sprintf('tk = %d', frnumbers(ixx(whichfr))), 'fontsize',20);
    axis square
    axis off
end


m1p = m1; m1p.yleft = m1p.yleft(ixx); m1p.yright = m1p.yright(ixx);
m1p.linl = 'bo'; m1p.linr = 'kd';

m2p = m2; m2p.yleft = m2p.yleft(ixx); m2p.yright = m2p.yright(ixx);
m2p.linl = 'bo'; m2p.linr = 'kd';

%
subplot(3,Nixx,((1:Nixx)+Nixx));
hold on
plotleftright(frnumbers(ixx), m1p);
subplot(3,Nixx,((1:Nixx)+2*Nixx));
hold on
plotleftright(frnumbers(ixx), m2p);

%% figure 6. Anglegram evolution

ixx = fix(linspace(1,length(bx),10));
Nixx = length(ixx);

figure(6);
set(gcf, 'Position', get(0, 'ScreenSize'))
%
bb = bbox2axis(T(ixx,:).BoundingBox);
bb = [min(bb(:,1)) max(bb(:,2)) min(bb(:,3)) max(bb(:,4))];
for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    figure(6)
    subplot(3,Nixx,whichfr)
    plotBoundariesAndPoints(frameinfo.X, frameinfo.initboundy, ...
        frameinfo.outboundy,'m-');
    plotBoundariesAndPoints([],[],cornies{ixx(whichfr)});
    
    axis(bb);
    title(sprintf('tk = %d', frnumbers(ixx(whichfr))));
    axis square
    axis off
    
    subplot(3,Nixx, Nixx + whichfr)
    imagesc(angie{ixx(whichfr)}');
    axis image
    axis off
    
    subplot(3,Nixx, 2*Nixx+whichfr)
    plot(1:64, anglesumve{ixx(whichfr)},...
        minlocations{ixx(whichfr)}, minval{ixx(whichfr)},'dr')
    plotHorzLine(1:64, [mam{ixx(whichfr)} mam{ixx(whichfr)}-stam{ixx(whichfr)}]);
    grid on
    ylim([50 250]);
    xlim([1 64])
    axis square
    
end


%% figure 7. (OPTIONAL)

ixx = fix(linspace(1,length(bx),10));
Nixx = length(ixx);

f33=figure(33);
set(gcf, 'Position', get(0, 'ScreenSize'))
subplot(3,Nixx,((1:Nixx)+Nixx));
plotleftright(frnumbers, m1);

%
bb = bbox2axis(T(ixx,:).BoundingBox);
bb = [min(bb(:,1)) max(bb(:,2)) min(bb(:,3)) max(bb(:,4))];
for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    figure(33)
    subplot(3,Nixx,whichfr)
    
    plotBoundariesAndPoints(frameinfo.X, frameinfo.initboundy, ...
        frameinfo.outboundy,'m-');
    plotBoundariesAndPoints([],[],cornies{ixx(whichfr)});
    
    axis(bb);
    title(sprintf('tk = %d', frnumbers(ixx(whichfr))), 'fontsize', 16);
    axis off
    
    subplot(3,Nixx, 2*Nixx+whichfr)
    plot(1:64, anglesumve{ixx(whichfr)},...
        minlocations{ixx(whichfr)}, minval{ixx(whichfr)},'dr',...
        'linewidth',1.5);
    plotHorzLine(1:64, [mam{ixx(whichfr)} mam{ixx(whichfr)}-stam{ixx(whichfr)}]);
    grid on
    ylim([50 230]);
    xlim([1 64]);
    title(sprintf('Corners: %d', length(minval{ixx(whichfr)})),...
        'fontsize',20);
end


m1p = m1; m1p.yleft = m1p.yleft(ixx); m1p.yright = m1p.yright(ixx);
m1p.linl = 'bo'; m1p.linr = 'kd';

%
subplot(3,Nixx,((1:Nixx)+Nixx));
hold on
plotleftright(frnumbers(ixx), m1p);

%% dodgy fig 7

ixx = fix(linspace(1,length(bx),8));
Nixx = length(ixx);

for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    subplot(3,Nixx, 2*Nixx+whichfr)
    plot(1:64, anglesumve{ixx(whichfr)},...
        minlocations{ixx(whichfr)}, minval{ixx(whichfr)},'dr',...
        'linewidth',1.5);
    plotHorzLine(1:64, [mam{ixx(whichfr)} mam{ixx(whichfr)}-stam{ixx(whichfr)}]);
    grid on
    ylim([50 230]);
    xlim([1 64]);
    title(sprintf('Corners: %d', length(minval{ixx(whichfr)})),...
        'fontsize',20);
end

%% dodgier fig 7

ixx = fix(linspace(1,length(bx),8));
Nixx = length(ixx);

subplot(3,Nixx,((1:Nixx)+2*Nixx));
hold on
for jx=1:length(frnumbers)
    x = repmat(frnumbers(jx), size(minval{jx},1), size(minval{jx},2));
    gca;
    yyaxis left
    plot(x, minval{jx}, '*', 'markersize', 8, 'linewidth', 1.5);
    grid on;
    ylabel('Angles at corners (DEG)', 'fontsize',16);
    yyaxis right
    plot(x(1), length(minval{jx}), '+-', 'markersize', 8, 'linewidth', 1.5);
    grid on;
    ylabel('Number of corners (DEG)', 'fontsize',16);
end

%%
for whichfr = 1:Nixx
    load(fullfile(rootdir, a{whichTrackIdx}, b{ixx(whichfr)}));
    
    subplot(3,Nixx, 2*Nixx+whichfr)
    plot(1:64, anglesumve{ixx(whichfr)},...
        minlocations{ixx(whichfr)}, minval{ixx(whichfr)},'dr',...
        'linewidth',1.5);
    plotHorzLine(1:64, [mam{ixx(whichfr)} mam{ixx(whichfr)}-stam{ixx(whichfr)}]);
    grid on
    ylim([50 230]);
    xlim([1 64]);
    title(sprintf('Corners: %d', length(minval{ixx(whichfr)})),...
        'fontsize',20);
end

%%

function scattpoinstandmeans(ptscirc, circ)
marksize = 15;
scatter(ptscirc(:,1), ptscirc(:,2), marksize,'.', 'linewidth', 1.5);
hold on
scatter(circ(:,1), circ(:,2), 4*marksize,'dk', 'linewidth', 1.5);
grid on
end