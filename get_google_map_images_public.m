% Requires a Google Static Maps API Key

path = [];
path.final = ''; % put working directory here
path.code = fullfile(path.final,'Code');
path.data = fullfile(path.final,'Data');
path.plot = fullfile(path.final,'Plots');
addpath(genpath(path.code))
cd(path.final)

%% Load data table

Data = readtable(fullfile(path.data,'Data_Table_01.csv'));

% data selection
Data = Data(~isnan(Data.cs_mn_avg_ol),:);
Data = Data(~ismember(Data.state,{'CA','WA','NV','ID','UT','NM','CO','KS'}),:); % eastern US

Data.cs_mn_grd_ol = str2double(Data.cs_mn_grd_ol);
Data.gcs_mn_grd_ol = str2double(Data.gcs_mn_grd_ol);

Elementary_Only = true;
if Elementary_Only
    Data = Data(strcmp(Data.school_level,'Elementary'),:);
end


%% Plot school scatters

States = unique(Data.state);
cols = parula(length(States)+2);
cols = cols(1:length(States),:);
cols = cols(randperm(size(cols,1)),:);

figure()
set(gcf,'Position',[1, 913, 1920, 956])
geoaxes
hold on
for st = 1:length(States)
    state = States{st};
    ind = strcmp(Data.state,state);
    plt = geoscatter(Data.lat(ind),Data.lon(ind));
    plt.MarkerEdgeColor = cols(st,:);
end

filename = fullfile(path.plot,'east_coast_schools_scatter.png');
saveas(gcf,filename)

%% Plot test score variables

% Levels = {'Elementary','Middle','High'};
% cols = lines(3);
% Vars = {'cs_mn_avg_ol','cs_mn_coh_ol','cs_mn_grd_ol','cs_mn_mth_ol';...
%         'gcs_mn_avg_ol','gcs_mn_coh_ol','gcs_mn_grd_ol','gcs_mn_mth_ol'};
% 
% figure()
% set(gcf,'Position',[1, 913, 1920, 956])
% ax = cell(2,4); k=0;
% for i = 1:2
%     for j = 1:4
%         var = Vars{i,j}; k=k+1;
%         ax{i,j} = subplot(2,4,k); hold on;
%         for lev = 1:3
%             ind = strcmp(Data.school_level,Levels{lev});
%             Y = Data.(var)(ind); Y = Y(~isnan(Y));
%             yrange = [nanmin(Data.(var)), nanmax(Data.(var))];
%             dy = abs(diff(yrange))/100;
%             edges = (yrange(1):dy:yrange(2));
%             h = histcounts(Y,edges,'Normalization','probability');
%             x = edges(1:end-1)+(dx/2);
%             plt = plot(x,h,'col',cols(lev,:),'LineWidth',2);
%             title(var,'Interpreter','none')
%         end
%         if k==1
%             legend(Levels)
%         end
%     end
% end
% 
% filename = fullfile(path.plot,'test_score_variables_by_level.png');
% saveas(gcf,filename)

%%

zoom = 16; % 14, 15, 16 are reasonable options
maptype = 'satellite';
watermarkcrop = 25;
height = 227;
width = 227;
P = height*width*3; % number of pixels

Images = nan(size(Data,1),P);
foldername = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom)]);
if Elementary_Only
    foldername = [foldername '_elementary'];
end
if exist(foldername,'dir')==0
    mkdir(foldername)
end

for i = 1:size(Data,1)

        lat = Data.lat(i);
        lon = Data.lon(i);
        filename = fullfile(foldername,['school_' num2str(Data.school_code(i)) '.jpg']);

        if ~exist(filename,'file')
            [XX, YY, M, Mcolor] = get_google_map(lat,lon,'Zoom',zoom,'MapType',maptype,'Height',height+watermarkcrop,'Width',width);
            Im = reshape(Mcolor(M+1,:),[height+watermarkcrop,width,3]);
            Im = Im(1:height,:,:);
            
            % save image
            disp([num2str(i) '  ' num2str(100*i/size(Data,1)) '%'])
            
            imwrite(Im,filename,'jpg');
        end
end

%%

% Imvec = reshape(Im,[1,P]);
% Im = reshape(Imvec,[height,width,3]);

%% Save in folder with subfolders for classes

zoom = 16; % 14, 15, 16 are reasonable options
maptype = 'satellite';
watermarkcrop = 25;
height = 227;
width = 227;
P = height*width*3; % number of pixels
Elementary_Only = true;

ogfolder = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom)]);
if Elementary_Only
    ogfolder = [ogfolder '_elementary'];
end
imagelist=ls(ogfolder);
imagelist=imagelist(3:end,:);
imagelist=cellstr(imagelist);

%% Classify by test score

TestScores = nan(size(imagelist,1),1);
tempcode = strrep(cellstr(num2str(Data.school_code)),' ','');
for i = 1:size(imagelist,1)
    id = strsplit(imagelist{i},'.jpg');
    id = strsplit(id{1},'_');
    id = id{2};
    ind = strcmp(tempcode,id);
    if sum(ind)~=1
        disp('err')
    end
    TestScores(i) = Data.cs_mn_avg_ol(ind);
end

% thresh = [25,50,75,100];
% thresh = [50,100];
thresh = [1/3,2/3,1].*100;
threshvalues = nan(size(thresh));
% for zoom 14, 4 classes: -0.2372    0.0197    0.2810    1.5657
% for zoom 16 elementary, 4 classes: -0.2110    0.0523    0.3145    1.5657
for i = 1:length(thresh)
    threshvalues(i) = prctile(TestScores,thresh(i));
end
class = nan(size(TestScores));
for i = 1:length(TestScores)
    class(i) = find(TestScores(i)<=threshvalues,1);
end

newfolder = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom) '_' num2str(length(thresh)) '_classes']);
if Elementary_Only
    newfolder = [newfolder '_elementary'];
end

for i = 1:size(imagelist,1)
    Im = imread(fullfile(ogfolder,imagelist{i}));
    sf = ['class_' num2str(class(i))];
    subfolder = fullfile(newfolder,sf);
    if exist(subfolder,'dir')==0
        mkdir(subfolder)
    end
    fname = strsplit(imagelist{i},'.jpg');
    fname = [fname{1} '.jpg'];
    filename = fullfile(subfolder,fname);
    imwrite(Im,filename,'jpg');

    if mod(i,1000)==0
        disp(i)
    end
end

%% Classify by urbanicity

Urban = cell(size(imagelist,1),1);
tempcode = strrep(cellstr(num2str(Data.school_code)),' ','');
for i = 1:size(imagelist,1)
    id = strsplit(imagelist{i},'.jpg');
    id = strsplit(id{1},'_');
    id = id{2};
    ind = strcmp(tempcode,id);
    if sum(ind)~=1
        disp('err')
    end
    Urban{i} = Data.urbanicity{ind};
end

newfolder = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom) '_' 'urbanicity' '_classes']);
if Elementary_Only
    newfolder = [newfolder '_elementary'];
end

for i = 1:size(imagelist,1)
    Im = imread(fullfile(ogfolder,imagelist{i}));
    sf = ['class_' num2str(Urban{i})];
    subfolder = fullfile(newfolder,sf);
    if exist(subfolder,'dir')==0
        mkdir(subfolder)
    end
    fname = strsplit(imagelist{i},'.jpg');
    fname = [fname{1} '.jpg'];
    filename = fullfile(subfolder,fname);
    imwrite(Im,filename,'jpg');
    if mod(i,1000)==0
        disp(i)
    end
end

%% Make folder of images not uploaded to gogle drive


%% Extract Elementary zoom 14 from base dataset and save in subfolders

zoom = 14; % 14, 15, 16 are reasonable options
maptype = 'satellite';
watermarkcrop = 25;
height = 227;
width = 227;
P = height*width*3; % number of pixels
Elementary_Only = true;

ogfolder = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom)]);
imagelist=ls(ogfolder);
imagelist=imagelist(3:end,:);
imagelist=cellstr(imagelist);

keepind = zeros(size(imagelist));
tempcode = strrep(cellstr(num2str(Data.school_code)),' ','');
for i = 1:length(imagelist)
    id = strsplit(imagelist{i},'.jpg');
    id = strsplit(id{1},'_');
    id = id{2};
    keepind(i) = sum(strcmp(tempcode,id));
end
imagelist = imagelist(keepind==1);

TestScores = nan(size(imagelist,1),1);
tempcode = strrep(cellstr(num2str(Data.school_code)),' ','');
for i = 1:size(imagelist,1)
    id = strsplit(imagelist{i},'.jpg');
    id = strsplit(id{1},'_');
    id = id{2};
    ind = strcmp(tempcode,id);
    if sum(ind)~=1
        disp('err')
    end
    TestScores(i) = Data.cs_mn_avg_ol(ind);
end

thresh = [25,50,75,100];
threshvalues = nan(size(thresh));
% for zoom 14, 4 classes: -0.2372    0.0197    0.2810    1.5657
% for zoom 16 elementary, 4 classes: -0.2110    0.0523    0.3145    1.5657
% for zoom 14 elementary, 4 classes: -0.2292    0.0354    0.3050    1.5657
for i = 1:length(thresh)
    threshvalues(i) = prctile(TestScores,thresh(i));
end
class = nan(size(TestScores));
for i = 1:length(TestScores)
    class(i) = find(TestScores(i)<=threshvalues,1);
end

newfolder = fullfile(path.data,['google_maps_' maptype '_' num2str(height) '_' num2str(width) '_' num2str(zoom) '_' num2str(length(thresh)) '_classes']);
if Elementary_Only
    newfolder = [newfolder '_elementary'];
end

for i = 1:size(imagelist,1)
    Im = imread(fullfile(ogfolder,imagelist{i}));
    sf = ['class_' num2str(class(i))];
    subfolder = fullfile(newfolder,sf);
    if exist(subfolder,'dir')==0
        mkdir(subfolder)
    end
    fname = strsplit(imagelist{i},'.jpg');
    fname = [fname{1} '.jpg'];
    filename = fullfile(subfolder,fname);
    imwrite(Im,filename,'jpg');

    if mod(i,1000)==0
        disp(i)
    end
end


%%





