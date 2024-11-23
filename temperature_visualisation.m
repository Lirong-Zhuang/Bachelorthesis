% 文件路径
slot_file = 'D:\BACoding\2dsolution\daten\2d24112303_Slot_Phi.dat';
yoke_file = 'D:\BACoding\2dsolution\daten\2d24112303_Yoke_Phi.dat';
tooth_file = 'D:\BACoding\2dsolution\daten\2d24112303_Tooth_Phi.dat';

% 定义读取选项
opts = detectImportOptions(slot_file, 'FileType', 'text');
opts.VariableNamesLine = 1; % 表头在第一行
opts.DataLine = 2;          % 数据从第二行开始
opts.SelectedVariableNames = [1, 2]; % 选择需要的两列
opts = setvartype(opts, 'double');   % 将变量设置为数值类型

% 读取三个文件数据
slot_data = readmatrix(slot_file, opts);
tooth_data = readmatrix(tooth_file, opts);
yoke_data = readmatrix(yoke_file, opts);


% 提取时间和温度
time_slot = slot_data(:, 1);
temp_slot = slot_data(:, 2);

time_yoke = yoke_data(:, 1);
temp_yoke = yoke_data(:, 2);

time_tooth = tooth_data(:, 1);
temp_tooth = tooth_data(:, 2);

% 创建折线图
figure;
plot(time_slot, temp_slot, '-', 'LineWidth', 0.8, 'DisplayName', 'Slot');
hold on;
plot(time_yoke, temp_yoke, '-', 'LineWidth', 0.8, 'DisplayName', 'Yoke');
plot(time_tooth, temp_tooth, '-', 'LineWidth', 0.8, 'DisplayName', 'Tooth');
hold off;

% 设置图像标题和轴标签
title('Temperature Changes at Different Observers');
xlabel('Time (s)');
ylabel('Temperature (°C)');

% 添加图例和网格
legend('Location', 'Best');
grid on;

% 美化图形
set(gca, 'FontSize', 12); % 设置字体大小
