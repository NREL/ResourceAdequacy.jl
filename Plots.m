clear,clc,close all

LoadData = csvread('LA_residential_enduses.csv',1,1);

Fans = LoadData(:,1);
Pumps = LoadData(:,2);
Heating = LoadData(:,3);
Cooling = LoadData(:,4);
InteriorLighting = LoadData(:,5);
ExteriorLighting = LoadData(:,6);
WaterSystems = LoadData(:,7);
InteriorEquipment = LoadData(:,8);

TotalLoad = sum(LoadData,2);

%% Make a yearly plot of the data

figure(1)
set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5]);

n_rows = 2;
n_columns = 3;

subplot(n_rows,n_columns,1)
hold on
plot(Cooling)
plot(Heating)
plot(TotalLoad)
axis([1000 1240 0 3000])
xlabel('Hour of Year')
ylabel('Load (MW)')
title('Winter')

subplot(n_rows,n_columns,2)
hold on
plot(Cooling)
plot(Heating)
plot(TotalLoad)
axis([6500 6740 0 6000])
xlabel('Hour of Year')
ylabel('Load (MW)')
title('Peak Load')

subplot(n_rows,n_columns,4)
hold on
plot(Cooling)
plot(Heating)
plot(TotalLoad)
axis([3600 3840 0 4500])
xlabel('Hour of Year')
ylabel('Load (MW)')
title('Minimum Load')

subplot(n_rows,n_columns,5)
hold on
plot(Cooling)
plot(Heating)
plot(TotalLoad)
axis([4100 4340 0 4000])
xlabel('Hour of Year')
ylabel('Load (MW)')
title('Summer')


temp = plot(1:10,1:10,1:10,2:11,1:10,3:12);
c = get(temp,'Color');


h = zeros(3, 1);
h(1) = plot(NaN,NaN,'Color',[0 0.4470 0.7410]);
h(2) = plot(NaN,NaN,'Color',[0.85 0.325 0.0980]);
h(3) = plot(NaN,NaN,'Color',[0.9290 0.694 0.125]);


hL = legend(h, 'Cooling','Heating','Total Load');


newPosition = [0.65 0.4 0.15 0.15];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

%suptitle('TITLE')

saveas(gcf,'LoadData.png')


figure
h = area([Cooling Heating TotalLoad-Heating-Cooling]);
legend('Cooling','Heating', 'Total Load')
h(1).FaceColor = [0 0.4470 0.7410];
h(2).FaceColor = [0.85 0.325 0.0980];
h(3).FaceColor = [0.9290 0.694 0.125];

%Get a Winter plot
axis([200 320 0 3500]) %200 = 8am, 204 = noon
xticks([204 216 228 240 252 264 276 288 300 312])
xticklabels({'12PM','12AM','12PM','12AM','12PM','12AM','12PM','12AM','12PM','12AM'})
xtickangle(45)
xlabel('Hour of Day')
ylabel('Demand (MW)')
title('Typical Winter Load')
saveas(gcf,'LoadData_stackplot_winter.png')

%Get a summer plot
axis([6530 6650 0 6000]) %6528 is midnight
xticks([6540 6552 6564 6576 6588 6600 6612 6624 6636 6648])
xticklabels({'12PM','12AM','12PM','12AM','12PM','12AM','12PM','12AM','12PM','12AM'})
title('Peak Load')
saveas(gcf,'LoadData_stackplot_peak.png')

%% Box and whisker plots
close all

nonDispatchableLoads = Fans + Pumps + InteriorLighting + ExteriorLighting + WaterSystems + InteriorEquipment;
dispatchableLoads = Cooling + Heating;

f2 = figure(2);
h1 = subplot(1,2,1);
boxplot([dispatchableLoads,nonDispatchableLoads],'Labels',{'Controllable','Uncontrollable'},'Whisker',6)
ylim([-100 6000])
ylabel('Power (MW)')


h2 = subplot(1,2,2);
boxplot(TotalLoad,'Labels',{'Total'},'Whisker',6)
ylim([-100 6000])
set(gca,'YTickLabel',{' '})
%suptitle('????')


set(h1, 'Position', [h1.Position(1) + 0.02 h1.Position(2) h1.Position(3) h1.Position(4)]);
set(h2, 'Position', [h2.Position(1) - 0.02 h2.Position(2) h1.Position(3) h2.Position(4)]); %Not a type-o. Setting the width to the same as h1
set(f2, 'Position', [680 707 633 271]);

saveas(gcf,'LoadDataBoxPlot.png')


%% Plot the results
clear,clc,close all

UnservedDRdata = csvread('C:\Users\aklem\Desktop\PrasConfPaperResults\UnservedHours_DR.csv',1,0);
UnservedBaseCasedata = csvread('C:\Users\aklem\Desktop\PrasConfPaperResults\UnservedHoursBaseCase.csv',1,0);

%Should I do a sorted bar graph, or a histogram? or what?
histogram(UnservedBaseCasedata)
xlabel('Hours of lost load per MC iteration')
ylabel('Number of iterations with this amount')

figure
bar(sort(UnservedBaseCasedata))

f = figure
xlabel('Equivalent Firm Capacity (MW)')
ylabel('Mean annual LOLE (Hours)')
title('Equivalent Firm Capacity Determination')
hold on


EFC_results = csvread('C:\Users\aklem\Desktop\PrasConfPaperResults\EFC_final_results.csv',1,0);


xi = 0:0.1:max(EFC_results(:,1));
EFC_interp = interp1(EFC_results(:,1),EFC_results(:,2),xi);


meanAnnualUnservedLoadHours = mean(UnservedDRdata);

titles = {'\lambda = 0.25','\lambda = 0.5','\lambda = 0.75','\lambda = 1'};
y_axes = [4 1 0.3 0.1];
for i = 1:length(meanAnnualUnservedLoadHours)-1
    h{i} = subplot(length(meanAnnualUnservedLoadHours)-1,1,i);
    p{i} = get(h{i}, 'Position');
    hold on
    plot([0 1*max(EFC_results(:,1))], meanAnnualUnservedLoadHours(i+1)*[1,1],'--k','Linewidth',1)
    plot(xi,EFC_interp,'-','Linewidth',2,'Color',[0 0.4470 0.7410])
    plot(EFC_results(:,1),EFC_results(:,2),'o','Linewidth',2,'Color',[0 0.4470 0.7410])
    axis([0 4000 0 y_axes(i)])
    title(titles{i},'FontSize',14)
    
    [temp1, temp2] = min(abs(EFC_interp - meanAnnualUnservedLoadHours(i+1)));
    %plot([EFC_interp(temp2) 0],[EFC_interp(temp2) meanAnnualUnservedLoadHours(i+1)],'--k')
    plot([xi(temp2) xi(temp2)],[0 EFC_interp(temp2)],'--k')
    plot(xi(temp2),EFC_interp(temp2),'rx','MarkerSize',18,'Linewidth',2)
end

set(f, 'Position',[0 0 958 952]);
h3=axes('position',[0.10 0.11 0.9 0.8],'visible','off');
h_label=ylabel('Annual Average LOLE (Hrs)','visible','on','FontSize',20);
xlabel('Equivalent Firm Capacity (MW)','visible','on','FontSize',16);
saveas(gcf,'EFC_results.png')





