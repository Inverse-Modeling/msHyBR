function shat_grid = visual_no_Africa_no_sea(s, land, lon, lat, show)
%
% code to visualize the reconstructions on the map
%
n1 = length(lon); n2 = length(lat);
shat_grid = zeros(n1*n2,1);
shat_grid(land) = s;
shat_grid(setdiff(1:end,land)) = NaN;
shat_grid = flipud(reshape(shat_grid, n1,n2)');
shat_grid(:,160:end)=NaN;
%set(h, 'AlphaData', ~isnan(img_data))
%shat_grid(shat_grid==0)=NaN;
if show
    figure, h=imagesc(shat_grid);
    set(h, 'AlphaData', ~isnan(shat_grid));
end

[m,n] = size(shat_grid);
num_ticks = 4;

LAT_min = floor(min(lat(:)));
LAT_max = ceil(max(lat(:)));
LON_min = floor(min(lon(:)));
LON_max = ceil(max(lon(:)));

lat_ticks = linspace(LAT_max,LAT_min,num_ticks);
lon_ticks = linspace(LON_min,LON_max,num_ticks);
lat_ticks_location = linspace(0,m,num_ticks);
lat_ticks_location(1) = 1;
lon_ticks_location = linspace(0,n,num_ticks);
lon_ticks_location(1) = 1;

yticks(lat_ticks_location)
yticklabels(num2cell(round(lat_ticks)))

xticks(lon_ticks_location)
xticklabels(num2cell(round(lon_ticks)))

set(gca,'XAxisLocation','bottom','LineWidth',1)


% Boundary
[m,n] = size(shat_grid);
S_BW = zeros(n,m);
S_BW(land)= 1;
S_BW = flipud(S_BW');
BW = bwboundaries(S_BW,'holes');

hold on
for k = 1:length(BW)
    boundary = BW{k};
    if min(boundary(:,2))<160
    plot(boundary(:,2),boundary(:,1), 'k', 'Linewidth',3)
    end
end
hold off