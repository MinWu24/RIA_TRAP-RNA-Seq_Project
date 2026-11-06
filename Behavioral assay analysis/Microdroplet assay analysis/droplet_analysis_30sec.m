

% DROPLET ASSAY ANALYSIS
% Chris Fang-Yen 4/09

%%%%%%%%%%%%%%%%%%%%%% LOAD IMAGE DATA %%%%%%%%%%%%%%%%%%%%%% 


button = length(questdlg('Load new data?','','Yes (JPG)','Yes (MAT) ','No', 'Yes (JPG)') ) ;

if button == 10
    clear all;
    [filename,pathname]  = uigetfile({'*.mat'});
    loa
    d([pathname filename]);
elseif button == 9
    clear all;
    [filename,pathname]  = uigetfile({'*.jpg'}, 'Select one image file');
    fname = [pathname filename(1:end-10)];
    answer = inputdlg({'Start frame', 'End frame'}, [pathname filename], 1,{['1'],['3600']});
    
    istart = str2num(answer{1});
    iend = str2num(answer{2});
    numframes = iend - istart + 1;

    i = istart ;
    fname2=strcat(fname, num2str(i, '%06d'), '.jpg');
    img = single(imread(fname2)); 
    
    [ysize xsize ] = size(img);

    % select cropping for each droplet

%     numroi = 12;
    answer = inputdlg('Number of ROIs','Number of ROIs',1,{'12'});
    numroi = str2num(answer{1});

    cropx1 = zeros(numroi,1);
    cropy1 = zeros(numroi,1);
    cropx2 = zeros(numroi,1);
    cropy2 = zeros(numroi,1);

    ignore_worm = zeros(numroi,1);

    for k=numroi:-1:1
        i=istart;
        fname2=strcat(fname, num2str(i, '%06d'), '.jpg');
    %       fname2=strcat(fname, num2str(i, '%03d'), '.jpg');
        img = single(imread(fname2)); 
        figure(1);clf;
        imagesc(img, [0 max(max(img))]); colormap autumn; hold on;
        title(['worm #' num2str(k) ' select ROI: upper left then lower right or press enter to ignore']);
        tmp = ginput(1);
        if length(tmp)<2
            ignore_worm(k)=1;
            text(20,20,['Ignoring worm #' num2str(k)], 'Color', 'w');
%             pause(1);
        else
            ignore_worm(k)=0;
            cropx1(k)= floor(tmp(1));
            cropy1(k)  = floor(tmp(2));
            plot([1 xsize], [cropy1(k) cropy1(k)], '-b');
            plot([cropx1(k) cropx1(k)], [1 ysize], '-b');
            [cropx2(k) cropy2(k) ] = ginput(1);
            cropx2(k) = floor(cropx2(k));
            cropy2(k) = floor(cropy2(k));
            plot([1 xsize], [cropy2(k) cropy2(k)], '-b');
            plot([cropx2(k) cropx2(k)], [1 ysize], '-b');
            text(20,20,['worm #' num2str(k)], 'Color', 'w');
            pause(.1);
        end
    end

    % FIND MINIMUM IMG %%%%%%
    figure(1); clf;
    axis image;
    title('minimum image');
    imgmin = 255*ones(size(img));
    numframes = iend - istart + 1;
    skip = round(numframes/10); 
    close
    for i=istart:skip:iend
    %     i = istart + (j - 1);
        fname2=strcat(fname, num2str(i, '%06d'), '.jpg');
    %       fname2=strcat(fname, num2str(i, '%03d'), '.jpg');
        img = single(imread(fname2)); 
        imgmin = min(img,imgmin);

        imagesc(imgmin); hold on;
        title('minimum image');
    %     title(num2str(i));
%         pause(0.1);
    end
    
    imgmin_thresh = 10;
    imgmin_bw = (imgmin > imgmin_thresh);
    se = strel('ball',3,3);

    imgmin_dilate = single(imdilate(imgmin,se)); 

%     pause(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numframes = iend-istart + 1;
    mask_multiplier = 1.1;
    img_threshold =  0.1;
    do_showimg = 0;
    showimg_multiple = 100;

    bw_Area = zeros(numroi, numframes);
    bw_Centroid = zeros(numroi, numframes,2);
    bw_Centroid_r = zeros(numroi, numframes);
    bw_Eccentricity = zeros(numroi, numframes);

    % figure(2);
    k0 = 2;  % worm to display.  0 for no display.

    fil = fspecial('disk', 5);
    % fspecial('gaussian', 10, 5);
    % se = strel('ball',5,5);
    if do_showimg
        figure(2); clf;
    end
    close
    h = waitbar(0,'Please wait...');
    tic
    for j=numframes:-1:1 %change
%         if mod(j,10)==0
%             disp(['j = ' num2str(j) '/' num2str(numframes)]);
%             figure(1); clf;
%             imagesc(img); colormap gray;
%             title(strcat(strcat(num2str(j), '/'), num2str(numframes)), 'Interpreter', 'None');
%             axis image;
%         end
        
        steps = numframes;
        % computations take place here
        waitbar((numframes-j) / steps)
        
        
        i = istart + (j - 1);
        fname2=strcat(fname, num2str(i, '%06d'), '.jpg'); 
        img = single(imread(fname2));   % load image
        img = img - mask_multiplier*imgmin_dilate; % subtract background
        img(img<0)=0;  % set negative pixels to zero

    %     if do_showimg
    %         figure(1); clf;
    %         imagesc(img); colormap gray;
    %         title(strcat(strcat(num2str(j), '/'), num2str(numframes)), 'Interpreter', 'None');
    %         axis image;
    % %     end

        for k=numroi:-1:1 %change
            if ~ignore_worm(k)
            imgcrop = img(cropy1(k):cropy2(k),cropx1(k):cropx2(k));
%             colormap gray;
                masksize = min([80, cropy2(k)-cropy1(k), cropx2(k) - cropx1(k)]);
                tmp = zeros(size(imgcrop));
                tmp(1:masksize,1:masksize) = hann(masksize) * hann(masksize)';
                mask2 = circshift(tmp, [round(-masksize/2), round(-masksize/2)]);
            imgcrop1 = imfilter(imgcrop, fil, 'same');
            imgmax = max(max(imgcrop1));
            [yi xi] = find(imgcrop1 == imgmax,1);
            imgcrop1f = imgcrop1 .* circshift(mask2, [yi, xi]);
            imgcrop2 = imgcrop .* circshift(mask2, [yi, xi]);
            imgcrop1f_bw = (imgcrop1f>imgmax*img_threshold);
            imgcrop1f_bws = bwselect(imgcrop1f_bw,xi,yi,8);
            if do_showimg && mod(j,showimg_multiple)==1
                [ii jj] = ind2sub([3,4], k);
                p = sub2ind([4,3], jj, ii);
%                 subplot(3,4,p);
%                 imagesc(imgcrop1f_bws); hold on; colormap gray;
%                 plot(xi, yi, 'or');
%                 axis off;
%                 if k==1
%                     title(num2str(j));
%                 end
            end
            STATS = regionprops(bwlabel(imgcrop1f_bws), 'Area', 'Centroid', 'Eccentricity');
            bw_Area(k,j) = STATS.Area;
            bw_Centroid(k,j,:) = STATS.Centroid;
            bw_Eccentricity(k,j) = STATS.Eccentricity;
%             pause(0.1);
        end
        end
    end
    toc

end
close(h) 
% %%%%%%%%%%%%%%%%%%%%%%%%% detect turns %%%%%%%%%%%%%%%%%%%%%%%%
    do_pause = 0;  % pause after each ROI analysis?
    Eccentricity_filsize=3;  % filter size for smoothing eccentricity 
    Area_filsize = 10; % filter size for smoothing Area 
    Centroid_filsize = 10; % filter size for smoothing Centroid 
    
    peak_det_abs_threshold = 0.85;  %  absolute threshold for turn detection 
    peak_det_threshold = 0.15;  % minimum size of peak in eccentricity
    
    Centroid_r_threshold = .7;  % lower = more strict
    area_threshold =0.7;  % higher = more strict

    answer = inputdlg({'peak_det_abs_threshold', 'peak_det_threshold', 'Centroid_r_threshold', 'area_threshold' }, 'Parameters', 1,...
            {num2str(peak_det_abs_threshold),num2str(peak_det_threshold),num2str(Centroid_r_threshold), num2str(area_threshold)});
    peak_det_abs_threshold = str2num(answer{1});
    peak_det_threshold = str2num(answer{2});
    Centroid_r_threshold = str2num(answer{3});
    area_threshold = str2num(answer{4});
    
    bw_Centroid_r_filtered = zeros(numroi,numframes);
    bw_Eccentricity_filtered = zeros(numroi,numframes);
    bw_Area_filtered = zeros(numroi,numframes);
    droplet_radius = zeros(numroi,numframes);
    
    for kk=1:numroi  % Smoothing
        if ~ignore_worm(kk)
            img_radius = img(cropy1(kk):cropy2(kk),cropx1(kk):cropx2(kk));
            tmp1 = size(img_radius,1)/2 * ones(1,numframes);
            tmp2 = size(img_radius,2)/2 * ones(1,numframes);
            droplet_radius(kk) = (size(img_radius,1)+ size(img_radius,2))/4;
            bw_Centroid_r(kk,:) = ((bw_Centroid(kk,:,1)-tmp2).^2 + ((bw_Centroid(kk,:,2)-tmp1)).^2).^0.5/droplet_radius(kk);
            bw_Centroid_r_filtered(kk,:) = smooth(bw_Centroid_r(kk,:), Centroid_filsize);
            bw_Eccentricity_filtered(kk,:) = smooth(bw_Eccentricity(kk,:), Eccentricity_filsize);
            bw_Area_filtered(kk,:) = smooth(bw_Area(kk,:), Area_filsize);
            mean_Area(kk) = mean(bw_Area(kk,:));
            median_Area(kk) = median(bw_Area(kk,:));
        end
    end

    turndata = zeros(numroi, numframes);
    data_invalid = zeros(numroi, numframes);
    turndata_valid = zeros(numroi, numframes);
    
    for k=1:numroi  % detect turns
        if ignore_worm(k)
            data_invalid(k,:)=1;
        else
            [maxtab mintab] = peakdet(bw_Eccentricity_filtered(k,:),peak_det_threshold);
            if size(mintab,1) > 0
            flag = zeros(size(mintab,1));
                for t=1:size(mintab,1)
                     if bw_Eccentricity_filtered(k,mintab(t,1)) < peak_det_abs_threshold
                         turndata(k,mintab(t,1))=1;
                     end
                end
            end

            data_invalid(k,:) = data_invalid(k,:) | (bw_Centroid_r_filtered(k,:)> Centroid_r_threshold);
            data_invalid(k,:) = data_invalid(k,:) | (bw_Area_filtered(k,:) < area_threshold* median_Area(k));

            turndata_valid(k,:) = turndata(k,:) .* ~data_invalid(k,:);

%             figure(100+k);clf;
%             subplot(311);
%             plot([1:numframes], bw_Eccentricity_filtered(k,:)); hold on;
%             plot(find( turndata_valid(k,:)), .50*ones(size(find( turndata_valid(k,:)))), 'or');
%             plot(find(data_invalid(k,:))', .4*ones(size(find(data_invalid(k,:))))', '.r', 'MarkerSize', 1);

    %         plot([1:numframes], data_invalid(k,:), '-r');

%             xlim([0 numframes]);
%             ylim([.2 1]);
%             title({[fname ', worm #' num2str(k)];['red line = exclude frames'];['bw_Eccentricity_filtered: ' num2str(sum(turndata(k,:))) ' turns detected']}, 'Interpreter', 'None');

%             subplot(312);
%             plot([1:numframes], bw_Area_filtered(k,:)); hold on;
%             plot([1:numframes], area_threshold* median_Area(k) * ones(1,numframes), '-c');
%             xlim([0 numframes]);
%             title(['bw_Area w/ threshold'], 'Interpreter', 'None');
%              subplot(313);

%              plot([1:numframes], bw_Centroid_r_filtered(k,:)); hold on;
%             plot([1:numframes], Centroid_r_threshold * ones(1,numframes), '-c');
%             xlim([0 numframes]);
%             title(['bw_Centroid_r_filtered w/ threshold'], 'Interpreter', 'None');

            if do_pause 
               pause; 
            end
        end
    end
    

% calculate group choice indices 
halfperiod = 300; % changed by HeLiu. it was 300;-2015/12/23

tmp = 0:numframes-1;
frames1 = (floor(tmp/halfperiod)/2 == floor(floor(tmp/halfperiod)/2));
% frames1([1:300, 601:900])=1;   % frames with stimulus #1
frames2 = ~frames1;  % frames with stimulus #2

num_datasets=1;
num_turns_all = zeros(num_datasets,numroi,2);
choice_index_all = zeros(num_datasets,numroi);

% grouping = 1: '1-3, 4-6, 7-9, 10-12'
% 2: '1-6, 7-12'
% 3: '1-6, 10-12, 7-9, 13-15'  (15 worms)

answer = 5;

if numroi ==12
    grouping = length(questdlg('Select grouping','','1-3, 4-6, 7-9, 10-12','1-6, 7-12', '1-4, 5-8, 9-12','1-3, 4-6, 7-9, 10-12') );
end

if numroi== 15
    grouping = length(questdlg('Select grouping','','1-6, 10-12, 7-9, 13-15',  '1-6, 10-12, 7-9, 13-15') );
end

switch grouping
    case 20 % '1-3, 4-6, 7-9, 10-12'
        num_groups = 4;
        worms_per_group = 3;
        group = zeros(worms_per_group,num_groups);
        group(:,1) = 1:3;  % identify worms in each group
        group(:,2) = 4:6;
        group(:,3) = 7:9;
        group(:,4) = 10:12;
    case 9      % '1-6, 7-12'
        num_groups = 2;
        worms_per_group = 6;
        group = zeros(worms_per_group,num_groups);
        group(:,1) = 1:6;  % identify worms in each group
        group(:,2) = 7:12;  
    case 22  % '1-6, 10-12, 7-9, 13-15'
        num_groups = 4;
        worms_per_group = 6;  % actually 6 and 3
        group = zeros(worms_per_group,num_groups);
        group(:,1) = 1:6;  % identify worms in each group
        group(:,2) = [7:9 7:9];
        group(:,3) = [10:12 10:12];
        group(:,4) = [13:15 13:15];       
    case 14  % '1-4, 5-8, 9-12'
        num_groups = 3;
        worms_per_group = 4;  
        group = zeros(worms_per_group,num_groups);
        group(:,1) = 1:4;  % identify worms in each group
        group(:,2) = 5:8;
        group(:,3) = 9:12;
%     case 5  % 'other'
%         answer = inputdlg('Number of groups','Number of groups',1,{'4'});
%         num_groups = str2num(answer{1});
%         answer = inputdlg('Max worms per group','Max worms per group',1,{'3'});
%         worms_per_group = str2num(answer{1});
%         for j=1:num_groups
%             for k = 1:worms_per_group
%                 prompt{worms_per_group*(j-1)+k} = ['group ' num2str(j) ' worm ' num2str(k)];
%                 def{worms_per_group*(j-1)+k} = num2str(worms_per_group*(j-1)+k);
%             end
%         end
%         dlg_title = 'Group assignment';
%         num_lines = num_groups*worms_per_group;
%         answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        
end

choice_index_grouped = zeros(num_datasets,num_groups);

valid_group_frames = zeros(num_groups,2);
turn_data_all = zeros(num_datasets, numroi, numframes);
turn_data_all(1,:,:) = turndata_valid;
turn_rate_grouped = zeros(num_datasets,num_groups,2);

for h=1:num_datasets
    num_turns_all(h,:,1) = sum(squeeze(turn_data_all(h,:,:)) .* repmat(frames1, [numroi 1]),2);
    num_turns_all(h,:,2) = sum(squeeze(turn_data_all(h,:,:)) .* repmat(frames2, [numroi 1]),2);
%     turnrate_all(h,:,1) = num_turns_al
% l(h,:,1)/sum(~data_invalid(find(frames1)),2);
%     turnrate_all(h,:,2) = num_turns_all(h,:,2)/sum(~data_invalid(find(frames2)),2);
%     choice_index_all(h,:) = (num_turns_all(h,:,1) - num_turns_all(h,:,2))./(num_turns_all(h,:,1) + num_turns_all(h,:,2)); % not normalized by valid frame number
    for k=1:num_groups
        num_turns_grouped(h,k,1) = sum(num_turns_all(h,group(:,k),1));
        num_turns_grouped(h,k,2) = sum(num_turns_all(h,group(:,k),2));
        valid_group_frames(k,1) = sum(sum(~data_invalid(group(:,k),:) .* repmat(frames1, [worms_per_group 1])));
        valid_group_frames(k,2) = sum(sum(~data_invalid(group(:,k),:) .* repmat(frames2, [worms_per_group 1])));
        turn_rate_grouped(h,k,1) =  num_turns_grouped(h,k,1) /  valid_group_frames(k,1) ;
        turn_rate_grouped(h,k,2) =  num_turns_grouped(h,k,2) /  valid_group_frames(k,2) ;
        
        choice_index_grouped(h,k) = ( turn_rate_grouped(h,k,1) -  turn_rate_grouped(h,k,2))/ ...
                                            ( turn_rate_grouped(h,k,1) + turn_rate_grouped(h,k,2));
    end
end
%% HeLiu get individual choice index
turns_01=num_turns_all(1,:,1);
turns_02=num_turns_all(1,:,2);

choice_index_individual=(turns_01-turns_02)./(turns_01+turns_02)

%%

figure(10);clf;
subplot(311)
bar(sum(turndata_valid,2));
xlabel('worm');
ylabel('total number of turns');
title({[fname]}, 'Interpreter', 'None');
subplot(312);
% bar([choice_index_grouped ; manual_turn_index_grouped]');
bar([choice_index_grouped]', 'r');
switch grouping
    case 20 % '1-3, 4-6, 7-9, 10-12'
    xlabel('group (#1 = worms 1-3, #2 = worms 4-6, #3 = worms 7-9, #4 = worms 10-12)');
    case 9      % '1-6, 7-12'
    xlabel('group (#1 = worms 1-6, #2 = worms 7-12)');
    case 22  % '1-6, 10-12, 7-9, 13-15'
    xlabel('group (#1 = worms 1-6, #2 = worms 7-9, #3 = worms 10-12, #4 = worms 13-15)');
    case 14  % '1-4, 5-8, 9-12'
    xlabel('group (#1 = worms 1-4, #2 = worms 5-8, #3 = worms 9-12)');
end
ylabel('choice index');
title([ 'choice index = ' num2str(choice_index_grouped)], 'Interpreter', 'None');
subplot(313)
bar([choice_index_individual]', 'r');
ylabel('choice index individual');
% title([ 'choice index = ' num2str(choice_index_individual)], 'Interpreter', 'None');

if length(questdlg('Save this data?'))==3
    max_num_turns=max(sum(num_turns_all,3));
    [file1 path1] =  uiputfile({'*.xls'}, 'Choose filename for turn data (.xls)', [fname '.xls']);
    pause(0.1);
    [file2 path2] =  uiputfile({'*.txt'}, 'Choose filename for turn rate / choice index (.txt)', [path1 file1(1:end-4) '.txt']);
    pause(0.1);
    [file3 path3] =  uiputfile({'*.mat'}, 'Choose filename for MAT file (.mat)', [path1 file1(1:end-4) '.mat']);
    pause(0.1);
    tmp = zeros(max_num_turns,numroi);
    for k=1:numroi
        turns=find( turndata_valid(k,:));
        tmp(1:length(turns),k)=turns;
    end
    xlswrite([path1 file1],[1:numroi; tmp]);
    disp(['Saved turn data in ' [path1 file1]]); 
    dlmwrite([path2 file2],[10*turn_rate_grouped(:,:,1); 10*turn_rate_grouped(:,:,2);choice_index_grouped ], 'delimiter', ' ', 'newline', 'pc');
    disp([10*turn_rate_grouped(:,:,1); 10*turn_rate_grouped(:,:,2);choice_index_grouped ]);
    disp(['Saved turn rate and choice indices in ' [path2 file2]]); 
    save([path3 file3]);
    disp(['Saved workspace in ' [path3 file3]]); 
    %all the turn during the recording
    file4='turndata_all.xlsx';
    for i=1:10:(length(turndata)-9)
        for j=1:12
            turn_all(j,(i+9)/10)=sum(turndata_valid(j,i:i+9));
        end
    end
    turn_all=turn_all';
    
    xlswrite([path1 file4],turn_all);
    file5='choice_index_individual';
    xlswrite([path1 file5],choice_index_individual);
    %%
    msgbox('Data saved successfully')
end
