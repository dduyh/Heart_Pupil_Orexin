
folderPath = 'C:\Users\yihudu\Desktop\Yihui\data\m1271\';
load([folderPath 'cluster_indices_8.mat'])
%%
A=xlsread([folderPath 'lick_8.xlsx']);
data = A(2:end,2:end);

%%
figure(1);
ht = suptitle('Facial Expressions during different Emotion States');
set(ht,'Position', [0.5 -0.05 0],'fontname','Times New Roman','fontsize',18)

h1 = subplot(3,1,1);
set(h1, 'Position', [0.298    0.89    0.44    0.0312])

mymap = [255 127 14
44 160 44
214 39 40
148 103 189
140 86 75
227 119 194
127 127 127
188 189 34
23 190 207
255 127 14
44 160 44];
mymap = mymap./255;

imagesc(cluster_indices)
colormap(h1,mymap);
axis off

h2 = subplot(3,1,2);
set(h2, 'Position', [0.298    0.83    0.44    0.0312])

emotions = zeros(192000,1);
emotions(7311:20503) = 1;
emotions(28812:34452) = -1;
emotions(34452:57914) = 1;
emotions(65217:73542) = -1;
emotions(73542:79982) = 1;
emotions(96153:103344) = -1;
emotions(103344:119141) = 1;
emotions(127185:135799) = -1;
emotions(135799:150524) = 1;
emotions(159169:192000) = -1;

mymap = [0 0 1
    0.8 0.8 0.8
    1 0 0];
imagesc(emotions(6254:157872)')
colormap(h2,mymap);
axis off

h3 = subplot(3,1,3);
set(h3, 'Position', [0.1300    0.04    0.7750    0.7750])
imagesc(data,[0.48 0.7])
% imagesc(data,[0.4 0.7])
colormap(h3,hot)
hbar = colorbar('east');
set(hbar, 'Position', [0.7675    0.04    0.013    0.77])
axis square

%%
savefig([folderPath 'cluster_index_8.fig']);
f = getframe(gcf);
imwrite(f.cdata,[folderPath 'cluster_index_8.png']);
