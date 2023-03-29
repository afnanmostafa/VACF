%% VACF plotter from LAMMPS dump (.xyz) file
% Afnan Mostafa
% 03/28/2023
% This code takes LAMMPS dump file, start and final timestep, and timestep
% size as inputs and gives the vacf considering all time lags (change lags
% variable if you want fewer lags)

%% clear space
% clearvars -except vx vy vz %% debug purposes; keep commented out
clear
clc
close
rng('shuffle');

%% system inputs (change as you go)
file = "dumpdos.nve";
start_step = 200000;
final_step = 208000;
timestep = 0.001; %% in ps (compatible with LAMMPS metal units)

%% read file
final_time = timestep*final_step;
A = regexp(fileread(file),'\n','split');
atomsLine = find(contains(A,'NUMBER OF ATOMS'));
totAtoms = A{1,atomsLine(1)+1};
totAtoms = str2double(totAtoms);
init_header = 0;
whichline = find(contains(A,num2str(final_step)));
whichline0 = find(contains(A,num2str(start_step)));
start = whichline0(1)+7; %% 7 lines between TIMESTEP and data in input file
ending = whichline+7;    %% 7 lines between TIMESTEP and data in input file
timesteps = find(contains(A,'ITEM: ATOMS'));

%% init variables
atomNo = cell(1,1);
vx = cell(1,1);
vy = cell(1,1);
vz = cell(1,1);
avg_v0 = zeros(1,1);
avg_v = zeros(1,1);
mean_vacf = zeros(1,1);

%% getting vx, vy, vz from the file
for i = 1:length(timesteps)
    fid = fopen(file);
    s = textscan(fid,'%d %d %f %f %f %f %f %f',totAtoms,'headerlines',timesteps(i));
    fclose(fid);
    atomNo{1,i} = s{1};
    vx{1,i} = s{6};
    vy{1,i} = s{7};
    vz{1,i} = s{8};
end

%% lags
lags = (1:length(vx)); %% total lags = lags (1->N) + 0th lag; N = length(vx) 

%% lag = 0
for lx = 1:length(vx)
    dotProduct = vx{1,lx}.*vx{1,lx}+vy{1,lx}.*vy{1,lx}+vz{1,lx}.*vz{1,lx};
    avg_v0(lx,1) = mean(dotProduct);
end
vacf{1,1} = avg_v0./avg_v0(1);
% plot((0:final_time),vacf{1,1}) %% evolution of vacf vs lag=0
% hold on

%% lag 1 --> final lag
for ly = 1:length(lags)
    for lx = 1:length(vx)
        if lx <= ly
            dotProduct = vx{1,lx}.*vx{1,1}+vy{1,lx}.*vy{1,1}+vz{1,lx}.*vz{1,1};
            avg_v(lx,1) = mean(dotProduct);
        elseif lx > ly
            dotProduct = vx{1,lx}.*vx{1,(lx-ly)}+vy{1,lx}.*vy{1,(lx-ly)}+vz{1,lx}.*vz{1,(lx-ly)};
            avg_v(lx,1) = mean(dotProduct);
        end
    end
    vacf{1,ly+1} = avg_v./avg_v(1);   %% normalization
    
    %% uncomment following line if you want all vacf plots
%     plot((0:final_time),vacf{1,ly+1}) %% evolution of all vacf vs lags
end
diff = final_step*timestep - start_step*timestep;
time = linspace(0, diff,length(vacf{1,1}));

%% spatial averaging across all time lags
cum_vacf = 0;
for g = 1:length(lags)
    for h = 1:length(vacf)
        cum_vacf = cum_vacf + vacf{1,h}(g);
    end
    mean_vacf(g,1) = cum_vacf/length(vacf);
    cum_vacf = 0;
end

%% plotting vacf with samples + fluctuation line about vacf = 0
p = plot(time,mean_vacf,'b-','DisplayName','Normalized VACF','LineWidth',2.5);
hold on
% psamples = plot(time(1:20:end),mean_vacf(1:20:end),'bo','DisplayName','Normalized VACF','LineWidth',2,...
%     'MarkerSize',6,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0 0 0]);
flucLine = yline(0,'Parent',gca,'Color',[0 0 0],'FontWeight','bold','LineStyle',...
    '--',...
    'LineWidth',1.5,...
    'FontName','Garamond',...
    'FontSize',18,...
    'Label',{'Fluctuations about 0'});

%% plotting features
legend([p,flucLine],'Normalized VACF','Stable Fluctuation Line','Location','northeast')
set(gca,'FontName','Garamond','FontSize',20,'FontWeight','bold',...
    'LineWidth',2,'XMinorTick','off','YMinorTick','off');
set(legend,...
    'LineWidth',2,...
    'FontSize',20,...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
    'FontWeight','bold');
ylabel('Normalized Velocity Auto-Correlation Function (VACF)','FontWeight','bold',...
    'FontSize',24,...
    'FontName','Garamond');
xlabel('Time (ps)',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Garamond');
%% close and get some rest
fclose('all');
%% END