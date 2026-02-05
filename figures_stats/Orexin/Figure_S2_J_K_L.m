%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

sema_folder = {'m2070\Apr_15_2025';
    'm2072\Apr_15_2025';
    'm2151\Apr_15_2025';
    'm2152\Apr_15_2025';
    'm2154\Apr_15_2025'};

glu_folder = {'m2070\May_01_2025';
    'm2072\May_01_2025';
    'm2151\May_01_2025';
    'm2152\May_01_2025'};

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration_sema = 4500;   % 1260 seconds.

trace_duration_glu = 3600;   % 1260 seconds.

pre_duration = 1300;

fpass_sema=[6 11.5;
    6 12;
    7 12;
    6 12;
    6 11];

fpass_glu=[5 10.2;
    5 10.5;
    5 10.5;
    5 10.5];

trials_num_sema = size(sema_folder,1);
trials_num_glu = size(glu_folder,1);

HR_sema = NaN(trials_num_sema,trace_duration_sema*Sample_Rate);    
HR_glu = NaN(trials_num_glu,trace_duration_glu*Sample_Rate);

%%

for I=1:size(sema_folder,1)

    Data_Folder = [Directory sema_folder{I} '\'];
    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals

    Injection_onset = find(datas(:,4),1);
    trace_start = Injection_onset-pre_duration*Sample_Rate;
    trace_end = trace_start+trace_duration_sema*Sample_Rate-1;

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    fpass=fpass_sema(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;
%     figure;
%     findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:size(ECG_raw,2),'nearest','extrap');
    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);
    HR_sema(I,:) = heartRate_bpm_interp_smooth;

end

%%

for I=1:size(glu_folder,1)

    Data_Folder = [Directory glu_folder{I} '\'];
    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals

    Injection_onset = find(datas(:,4),1);
    trace_start = Injection_onset-pre_duration*Sample_Rate;
    trace_end = trace_start+trace_duration_glu*Sample_Rate-1;

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    if I==2

        fpass1 = [8 12];
        fpass2 = [3 7];
        fpass3 = [7 11];

        % Injection_onset = find(datas(data_onset:data_offset,4),1);
        Injection_onset = 1600000-trace_start+1;

        bufferLen = round(300 * Sample_Rate);

        segment1 = ECG_raw(1:Injection_onset + bufferLen);
        segment2 = ECG_raw(Injection_onset - bufferLen + 1:Injection_onset + round(400 * Sample_Rate) + bufferLen);
        segment3 = ECG_raw(Injection_onset - bufferLen + round(400 * Sample_Rate) + 1:end);

        ECG_seg1 = bandpass(segment1, fpass1, Sample_Rate);
        ECG_seg2 = bandpass(segment2, fpass2, Sample_Rate);
        ECG_seg3 = bandpass(segment3, fpass3, Sample_Rate);

        valid_seg1 = ECG_seg1(1:Injection_onset);
        valid_seg2 = ECG_seg2(bufferLen+1:round(400 * Sample_Rate) + bufferLen);
        valid_seg3 = ECG_seg3(bufferLen+1:end);

        ECG_Bandpass = [valid_seg1, valid_seg2, valid_seg3];

    else

        fpass=fpass_glu(I,:);

        ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    end


    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:size(ECG_raw,2),'nearest','extrap');
    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);
    HR_glu(I,:) = heartRate_bpm_interp_smooth;

end

%%

xlims = (-pre_duration*Sample_Rate+1:(trace_duration_glu-pre_duration)*Sample_Rate)/Sample_Rate;

fig = figure(1);
set(fig, 'Position', [2655 665 2017 593]);

axes('Position', [0.1300 0.1100 0.2134 0.8099]);
hold on
patch('XData',[0, 0, 60, 60],'YData',[300 700 700 300],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.5);
% patch('XData',[0, 0, 60, 60],'YData',[450 600 600 450],'EdgeColor','none','FaceColor','#F6D7B4');

patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_glu)+std(HR_glu)/sqrt(trials_num_glu) fliplr(mean(HR_glu)-std(HR_glu)/sqrt(trials_num_glu))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',1);
plot(xlims,mean(HR_glu),'Color',[255 128 128]./255,'LineWidth',2)
% line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);

% xticks(-200:200:800);
% yticks(400:50:780);

xlim([-pre_duration trace_duration_glu-pre_duration])
% xlim([-1300 3200])
ylim([300 650])
% title({'Heart Rate Before/After Propranolol Injection',''},'FontSize',20,'FontWeight','bold','color','k')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

%%

xlims = (-pre_duration*Sample_Rate+1:(trace_duration_sema-pre_duration)*Sample_Rate)/Sample_Rate;

figure(1);

axes('Position', [0.4108 0.1100 0.2134 0.8099]);
hold on
patch('XData',[0, 0, 60, 60],'YData',[400 780 780 400],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.5);
% patch('XData',[0, 0, 60, 60],'YData',[450 600 600 450],'EdgeColor','none','FaceColor','#F6D7B4');

patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_sema)+std(HR_sema)/sqrt(trials_num_sema) fliplr(mean(HR_sema)-std(HR_sema)/sqrt(trials_num_sema))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',1);
plot(xlims,mean(HR_sema),'Color',[255 128 128]./255,'LineWidth',2)
% line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);

% xticks(-200:200:800);
% yticks(400:50:780);

xlim([-pre_duration trace_duration_sema-pre_duration])
ylim([450 650])
% title({'Heart Rate Before/After Propranolol Injection',''},'FontSize',20,'FontWeight','bold','color','k')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

%%

peak_HR_glu = min(HR_glu(:,1360*Sample_Rate+1:2360*Sample_Rate),[],2) - mean(HR_glu(:,pre_duration*Sample_Rate-1000*Sample_Rate+1:pre_duration*Sample_Rate),2);

peak_HR_sema = max(HR_sema(:,end-1000*Sample_Rate+1:end),[],2) - mean(HR_sema(:,pre_duration*Sample_Rate-1000*Sample_Rate+1:pre_duration*Sample_Rate),2);

HR_mean = [mean(peak_HR_glu) mean(peak_HR_sema)];
HR_sem = [std(peak_HR_glu)/sqrt(trials_num_glu) std(peak_HR_sema)/sqrt(trials_num_sema)];

figure(1);

axes('Position', [0.6916 0.1100 0.1 0.8099]);
hold on

b = bar(HR_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [222, 202, 224] / 255;
b.CData(2,:) = [198, 214, 234] / 255;

for k = 1:size(peak_HR_glu,1)
    plot(1,peak_HR_glu(k),'marker','o','markersize',5,...
        'markeredgecolor','#8B5C9E','markerfacecolor','#8B5C9E',...
        'linestyle','none');
end

for k = 1:size(peak_HR_sema,1)
    plot(2,peak_HR_sema(k),'marker','o','markersize',5,...
        'markeredgecolor','#719DC9','markerfacecolor','#719DC9',...
        'linestyle','none');
end

errorbar(1:2,HR_mean,HR_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

line([1 2], [200, 200], 'Color', 'k', 'LineWidth', 2);
text(1.5, 210, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2])
xticklabels({'GLUC','SEMA'})

ylim([-500 250])
xlim([0.3 2.7])

title({'Δ Heart Rate'},'FontSize',18,'FontWeight','bold')
ylabel('Δ bpm','FontSize',15,'FontWeight','bold');

hold off

[h_HR_glu_sema, p_HR_glu_sema, ~, stats_HR_glu_sema] = ttest2(peak_HR_glu,peak_HR_sema,'Tail','left')



