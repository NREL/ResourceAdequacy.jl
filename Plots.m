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


%% Box and whisker plots
close all
figure(2)
hold on
%boxplot(TotalLoad,Cooling)
nonDispatchableLoads = Fans + Pumps + InteriorLighting + ExteriorLighting + WaterSystems + InteriorEquipment;
dispatchableLoads = Cooling + Heating;
boxplot([dispatchableLoads,nonDispatchableLoads,TotalLoad],'Labels',{'dis','nondisp','{Total Load}'},'Whisker',6)
title('Compare Random Data from Different Distributions')






