clc;
clear;
% 加载和查看数据
load('data\imu_data.mat');
load('data\ground_truth_data.mat');
load('data\leg_data.mat');
%% 主函数
%将timeseries转为矩阵 
a_R = [-1,0,0;0,-1,0;0,0,1];
% a_R = eye(3);
acce_ = transpose(imu_data{:, 2:4}); %线加速度
gyro_ = transpose(imu_data{:, 5:7}); %角速度
position_true = transpose(ground_truth_data{:, 2:4}); %实际位置
R_true =  transpose(ground_truth_data{:, 5:8}); %实际姿态
p_ = transpose((ground_truth_data{1, 2:4}))+[0;0;0.34];   %初始位置
v_ = [0; 0; 0];   %初始速度
R_ = quaternionToRotationMatrix(ground_truth_data{1, 5:8}); %获得初始姿态矩阵
ba_ = [1; 1; 1] * 3.552e-04 * 0;   %线加速度偏移
bg_ = [1; 1; 1] * 1.655e-05 * 0;   %角速度偏移
gravity_ = [0; 0; -9.8]; %重力加速度  
% 确定传感器时间向量
dt = 1/765;      % 采样时间间隔
end_time = 52580*dt;    %采样时间为8s
time = 0:dt:end_time; %从第1个数据开始
% 确定真值采样时间
dt_true = 1/63;      % 采样时间间隔
end_time = 4400*dt_true;    %采样时间为8s
time_true = 0:dt_true:end_time; %从第1个数据开始
% 初始化数组以存储每个时间步的结果
num_steps = length(time);
num_steps_true = length(time_true);
p_history = zeros(3, num_steps);
v_history = zeros(3, num_steps);
R_history = zeros(3, 3, num_steps);
o_history = zeros(3, num_steps); 
o_true = zeros(3, num_steps_true); % 初始化用于存储对比的真实欧拉角数组
position_true_resampled = position_true(:,1:num_steps_true); %采样真实位置
% velocity_true_resampled = velocity_true(:,0:num_steps_true); %采样真实速度
R_true_resampled = R_true(:,1:num_steps_true); %采样真实姿态
q_history = zeros(4, num_steps); %四元数

%% 采样腿部数据
s_p0_array = zeros(3, num_steps); %初始化腿部距离身体的位置
s_p1_array = zeros(3, num_steps);
s_p2_array = zeros(3, num_steps);
s_p3_array = zeros(3, num_steps);

d0 = [0;0;0]; %腿在inertial frame下的位置
d1 = [0;0;0];
d2 = [0;0;0];
d3 = [0;0;0];

s_p0_speed = transpose(leg_data{1:num_steps, 5:7});   %腿0速度(body)
s_p1_speed = transpose(leg_data{1:num_steps, 15:17}); %腿1速度
s_p2_speed = transpose(leg_data{1:num_steps, 25:27}); %腿2速度
s_p3_speed = transpose(leg_data{1:num_steps, 35:37}); %腿3速度

d0_history = transpose(leg_data{1:num_steps, 2:4});   %腿0位置(body)
d1_history = transpose(leg_data{1:num_steps, 12:14}); %腿1位置
d2_history = transpose(leg_data{1:num_steps, 22:24}); %腿2位置
d3_history = transpose(leg_data{1:num_steps, 32:34}); %腿3位置

contact_0 = transpose(leg_data{1:num_steps, 11}); %腿0接触判断
contact_1 = transpose(leg_data{1:num_steps, 21}); %腿1接触判断
contact_2 = transpose(leg_data{1:num_steps, 31}); %腿2接触判断
contact_3 = transpose(leg_data{1:num_steps, 41}); %腿3接触判断

v_contact_history = zeros(3, num_steps); %腿的历史速度
d_history = zeros(3, num_steps); 
%% 初始化的卡尔曼滤波参数
A = [eye(3),eye(3)*dt,zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),eye(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)];     % 状态转移矩阵
R = [eye(3)*0.01,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),eye(3)*0.01,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3, 4);
     zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),eye(4)*0.001];     % 观测噪声协方差矩阵 %25
Q = [eye(3)*0.00001,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),eye(3)*0.00001,zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01];     % 状态转移噪声协方差矩阵
P = [ eye(3)*0.001,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3);
      zeros(3),eye(3)*0.001,zeros(3),zeros(3),zeros(3),zeros(3);
      zeros(3),zeros(3),eye(3)*0.001,zeros(3),zeros(3),zeros(3);
      zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3);
      zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3);
      zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.001];     % 状态转移噪声协方差矩阵
 
big_noise = 100; %抬腿大噪声
normal_noise = 0.001; %正常噪声
%% 真*主函数
for i = 1:num_steps_true
    o_true(:, i) =  quaternionToEuler(transpose( R_true_resampled(:,i) ));
end

for i = 1:num_steps
     rotationMatrix = eulerAnglesToRotationMatrix((gyro_(1:3, i) - bg_) * dt);  % 欧拉角转旋转矩阵
     R_ = R_ * rotationMatrix;% 假设欧拉角顺序为 ZYX
%     R_ = quaternionToRotationMatrix(transpose( R_true_resampled(:,i) )); %用真值测试
     p_ = p_ + v_ * dt + 0.5 * gravity_ * dt * dt + 0.5 * R_ * (acce_(1:3, i) - ba_) * dt * dt; %计算位置
     v_ = v_ + R_ * (acce_(1:3, i) - ba_) * dt + gravity_ * dt; % 计算速度

     %观测矩阵(接触版)
     H = [-1*eye(3),zeros(3),eye(3),zeros(3),zeros(3),zeros(3);
          -1*eye(3),zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
          -1*eye(3),zeros(3),zeros(3),zeros(3),eye(3),zeros(3);
          -1*eye(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3);
          zeros(3),contact_0(i)*eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
          zeros(3),contact_1(i)*eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
          zeros(3),contact_2(i)*eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
          zeros(3),contact_3(i)*eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
          0,0,0,   0,0,0, 0,0,contact_0(i), 0,0,0, 0,0,0, 0,0,0;
          0,0,0,   0,0,0, 0,0,0, 0,0,contact_1(i), 0,0,0, 0,0,0;
          0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,contact_2(i), 0,0,0;
          0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,contact_3(i);];     % 观测矩阵
    
%      %观测矩阵(噪声版)
%      H = [-1*eye(3),zeros(3),eye(3),zeros(3),zeros(3),zeros(3);
%           -1*eye(3),zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
%           -1*eye(3),zeros(3),zeros(3),zeros(3),eye(3),zeros(3);
%           -1*eye(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3);
%           zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
%           zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
%           zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
%           zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
%           0,0,0,   0,0,0, 0,0,contact_0(i), 0,0,0, 0,0,0, 0,0,0;
%           0,0,0,   0,0,0, 0,0,0, 0,0,contact_1(i), 0,0,0, 0,0,0;
%           0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,contact_2(i), 0,0,0;
%           0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,contact_3(i);];     % 观测矩阵
     
     %观测值
     delta_d0_world = R_*d0_history(:,i); %世界坐标系下腿部的位置差
     delta_d1_world = R_*d1_history(:,i);
     delta_d2_world = R_*d2_history(:,i);
     delta_d3_world = R_*d3_history(:,i);
     Z = [delta_d0_world;  %世界坐标系下腿部的位置差
          delta_d1_world;
          delta_d2_world;
          delta_d3_world;
          -R_*(cross((gyro_(1:3, i) - bg_),d0_history(:,i)) + s_p0_speed(:, i)); %世界坐标系下的速度
          -R_*(cross((gyro_(1:3, i) - bg_),d1_history(:,i)) + s_p1_speed(:, i));
          -R_*(cross((gyro_(1:3, i) - bg_),d2_history(:,i)) + s_p2_speed(:, i));
          -R_*(cross((gyro_(1:3, i) - bg_),d3_history(:,i)) + s_p3_speed(:, i));
          0;   %世界坐标系下腿部的z方向位置
          0;
          0;
          0]; 
           
      
%       %根据抬腿调整噪声
%        contacts = [contact_0(i), contact_1(i), contact_2(i), contact_3(i)];
%        [Q, R] = adjustNoise(Q, R, contacts, big_noise, normal_noise);
   
      % 卡尔曼滤波
      % 预测步骤
      x_hat_minus = [p_;v_;d0;d1;d2;d3];
      P_minus = A * P * A' + Q;
      % 更新步骤
      K = P_minus * H' / (H * P_minus * H' + R);
      x_hat = x_hat_minus + K * (Z - H * x_hat_minus);
%       x_hat = x_hat_minus ;
      P = (eye(18) - K * H) * P_minus;
      % 存储滤波后的状态
      p_history(:, i) = x_hat(1:3,1);
      v_history(:, i) = x_hat(4:6,1);
      R_history(:, :, i) = R_;
      
      %把R转为欧拉角
     o_history(:, i) = a_R * rotMatrixToWorldEulerAngles(R_);
     q_history(:, i) = eulerToQuaternion(o_history(:, i)); %记录每一次的四元数
%      o_true(:, i) =  quaternionToEuler(transpose( R_true_resampled(:,i) ));
     %更新数据
     p_ = x_hat(1:3,1);
     v_ = x_hat(4:6,1);
     d0 = x_hat(7:9,1);
     d1 = x_hat(10:12,1);
     d2 = x_hat(13:15,1);
     d3 = x_hat(16:18,1);
     
     % 存储结果
     p_history(:, i) = a_R*(p_+[0;0;0.026]);
     v_history(:, i) = v_;
     R_history(:, :, i) = R_;
     d_history(:, i) = d0;
end

%% 需要打印的数据
ground_truth_data{1:num_steps_true,1} = transpose(time_true);
outcome_ground_truth = ground_truth_data{1:num_steps_true,1:8};
outcome_kf = transpose([time;p_history(:,:);q_history(:,:)]);
% 将 outcome_luen_luen 保存为 TXT 文件
writematrix(outcome_kf, 'KF', 'Delimiter', ' ');
% 将 outcome_ground_truth 保存为 TXT 文件
writematrix(outcome_ground_truth, 'Ground_Truth.txt', 'Delimiter', ' ');

% %% 打印均方误差
% % 计算均方误差（MSE）
% mse_x = mean((p_history(1, :) - position_true_resampled(1, :)).^2);
% mse_y = mean((p_history(2, :) - position_true_resampled(2, :)).^2);
% mse_z = mean((p_history(3, :) - position_true_resampled(3, :)).^2);
% 
% % 计算均方根误差（RMSE）
% rmse_x = sqrt(mse_x);
% rmse_y = sqrt(mse_y);
% rmse_z = sqrt(mse_z);
% 
% % 计算平均绝对误差（MAE）
% mae_x = mean(abs(p_history(1, :) - position_true_resampled(1, :)));
% mae_y = mean(abs(p_history(2, :) - position_true_resampled(2, :)));
% mae_z = mean(abs(p_history(3, :) - position_true_resampled(3, :)));
% 
% % 显示结果
% fprintf('X方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_x, rmse_x, mae_x);
% fprintf('Y方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_y, rmse_y, mae_y);
% fprintf('Z方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_z, rmse_z, mae_z);

% %% 打印均欧拉角方误差
% % 计算均方误差（MSE）
% mse_x = mean((o_history(3, :) - o_true(1, :)).^2);
% mse_y = mean((o_history(2, :) - o_true(2, :)).^2);
% mse_z = mean((o_history(1, :) - o_true(3, :)).^2);
% 
% % 计算均方根误差（RMSE）
% rmse_x = sqrt(mse_x);
% rmse_y = sqrt(mse_y);
% rmse_z = sqrt(mse_z);
% 
% % 计算平均绝对误差（MAE）
% mae_x = mean(abs(o_history(3, :) - o_true(1, :)));
% mae_y = mean(abs(o_history(2, :) - o_true(2, :)));
% mae_z = mean(abs(o_history(1, :) - o_true(3, :)));
% 
% % 显示结果
% fprintf('Roll方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_x, rmse_x, mae_x);
% fprintf('Pitch方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_y, rmse_y, mae_y);
% fprintf('Yaw方向上的MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_z, rmse_z, mae_z);

% 绘制位置图
figure;

subplot(3, 1, 1);
plot(time, p_history(1, :), 'b', time_true, position_true_resampled(1, :), 'r');
title('Position - X');
legend('Computed', 'True');

subplot(3, 1, 2);
plot(time, p_history(2, :), 'b', time_true, position_true_resampled(2, :), 'r');
title('Position - Y');
legend('Computed', 'True');

subplot(3, 1, 3);
plot(time, p_history(3, :), 'b', time_true, position_true_resampled(3, :), 'r');
title('Position - Z');
legend('Computed', 'True');

% 绘制姿态图
figure;
subplot(3, 1, 1);
plot(time, o_history(1, :), 'b',time_true, o_true(1, :), 'r');
title('Orientation - Roll');
legend('Computed', 'True');

subplot(3, 1, 2);
plot(time, o_history(2, :), 'b',time_true, o_true(2, :), 'r');
title('Orientation - Pitch');
legend('Computed', 'True');

subplot(3, 1, 3);
plot(time, o_history(3, :), 'b',time_true, o_true(3, :), 'r');
title('Orientation - Yaw');
legend('Computed', 'True');

% 绘制二维轨迹图
figure;

% 绘制每条线并分别设置属性
plot(position_true_resampled(1, :), position_true_resampled(2, :), 'r-', 'LineWidth', 1.5); % 红色，实线
hold on;
plot(p_history(1, :), p_history(2, :), 'b', 'LineWidth', 1.5); % 蓝色，虚线

% 添加标题和图例
title('2D Trajectory - XY Plane');
xlabel('X Position');
ylabel('Y Position');
legend('True', 'Computed');
% 释放 hold
hold off;

%% 调整噪声参数
    % 输入:
    %   Q - 状态转移噪声协方差矩阵
    %   R - 观测噪声协方差矩阵
    %   contacts - 接触状态数组
    %   big_noise - 大噪声值
    %   normal_noise - 正常噪声值
function [Q, R] = adjustNoise(Q, R, contacts, big_noise, normal_noise)
    % 调整噪声参数
    % 输入:
    %   Q - 状态转移噪声协方差矩阵
    %   R - 观测噪声协方差矩阵
    %   contacts - 接触状态数组
    %   big_noise - 大噪声值
    %   normal_noise - 正常噪声值
    % 输出:
    %   Q - 调整后的状态转移噪声协方差矩阵
    %   R - 调整后的观测噪声协方差矩阵
    for j = 0:3
        if contacts(j+1) == 0 
            Q(7+3*j:9+3*j, 7+3*j:9+3*j) = eye(3) * big_noise;
            R(1+3*j:3+3*j, 1+3*j:3+3*j) = eye(3) * big_noise;
            R(13+3*j:15+3*j, 13+3*j:15+3*j) = eye(3) * big_noise;
            R(25+j,25+j) = big_noise;
        else 
            Q(7+3*j:9+3*j, 7+3*j:9+3*j) = eye(3) * normal_noise;
            R(1+3*j:3+3*j, 1+3*j:3+3*j) = eye(3) * normal_noise;
            R(13+3*j:15+3*j, 13+3*j:15+3*j) = eye(3) * normal_noise;
            R(25+j,25+j) = normal_noise;
        end
    end
end
%% 欧拉角函数
function R = eulerAnglesToRotationMatrix(angles)
    % angles 是一个包含绕X轴、Y轴和Z轴的欧拉角的 1x3 向量

    % 提取欧拉角
    roll = angles(1); %绕x轴旋转
    pitch = angles(2);   %绕y轴旋转
    yaw = angles(3);  %绕z轴旋转

    % 构建绕X轴的旋转矩阵
    Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];

    % 构建绕Y轴的旋转矩阵
    Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
    
    % 构建绕Z轴的旋转矩阵
    Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];

    % 构建最终的旋转矩阵
    R = Rz * Ry * Rx;
end

%% 旋转矩阵转欧拉角
function euler_angles = rotMatrixToWorldEulerAngles(rotation_matrix)
    % rotation_matrix 是一个3x3的旋转矩阵

    % 将旋转矩阵转为欧拉角，使用 'ZYX' 旋转顺序
    euler_angles = transpose(rotm2eul(rotation_matrix, 'ZYX'));
    
    %把弧度以roll pitch yaw 的顺序发布
    yaw = euler_angles(1);
    pitch = euler_angles(2);
    roll = euler_angles(3);
    euler_angles = [roll;pitch;yaw];
    
    % 将弧度转为角度
    euler_angles = rad2deg(euler_angles);
    
end

%% 四元数转欧拉角
function gyro = quaternionToEuler(quaternion)
    %QUATERNIONTOEULER Converts quaternion to Euler angles (roll, pitch, yaw).
    % Input:
    %   x, y, z, w - Quaternion components
    % Output:
    %   roll - Rotation around x-axis
    %   pitch - Rotation around y-axis
    %   yaw - Rotation around z-axis
    x = quaternion(1);
    y = quaternion(2);
    z = quaternion(3);
    w = quaternion(4);
    % Normalize the quaternion
    norm_q = sqrt(x^2 + y^2 + z^2 + w^2);
    x = x / norm_q;
    y = y / norm_q;
    z = z / norm_q;
    w = w / norm_q;
    
    % Roll (x-axis rotation)
    sinr_cosp = 2 * (w * x + y * z);
    cosr_cosp = 1 - 2 * (x * x + y * y);
    roll = atan2(sinr_cosp, cosr_cosp);
    
    % Pitch (y-axis rotation)
    sinp = 2 * (w * y - z * x);
    if abs(sinp) >= 1
        pitch = sign(sinp) * pi / 2; % use 90 degrees if out of range
    else
        pitch = asin(sinp);
    end
    
    % Yaw (z-axis rotation)
    siny_cosp = 2 * (w * z + x * y);
    cosy_cosp = 1 - 2 * (y * y + z * z);
    yaw = atan2(siny_cosp, cosy_cosp);
    
    gyro = transpose([roll, pitch, yaw]);
  
    % 将弧度转为角度
    gyro = rad2deg(gyro);
end

%% 四元数转R矩阵
function R = quaternionToRotationMatrix(quaternion)
    %QUATERNIONTOROTATIONMATRIX Converts quaternion to rotation matrix.
    % Input:
    %   x, y, z, w - Quaternion components
    % Output:
    %   R - 3x3 rotation matrix
    x = quaternion(1);
    y = quaternion(2);
    z = quaternion(3);
    w = quaternion(4);
    % Normalize the quaternion
    norm_q = sqrt(x^2 + y^2 + z^2 + w^2);
    x = x / norm_q;
    y = y / norm_q;
    z = z / norm_q;
    w = w / norm_q;
    
    % Compute the rotation matrix elements
    R = [1 - 2*y^2 - 2*z^2, 2*x*y - 2*z*w, 2*x*z + 2*y*w;
         2*x*y + 2*z*w, 1 - 2*x^2 - 2*z^2, 2*y*z - 2*x*w;
         2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x^2 - 2*y^2];
end

%% R矩阵转四元数
function quaternion = rotationMatrixToQuaternion(R)
    %ROTATIONMATRIXTOQUATERNION Converts rotation matrix to quaternion.
    % Input:
    %   R - 3x3 rotation matrix
    % Output:
    %   x, y, z, w - Quaternion components
    
    % Ensure the matrix is orthogonal and its determinant is 1
    if abs(det(R) - 1) > 1e-6
        error('Input matrix is not a valid rotation matrix');
    end
    
    % Compute the trace of the matrix
    trace_R = trace(R);
    
    if trace_R > 0
        S = sqrt(trace_R + 1.0) * 2; % S = 4 * w
        w = 0.25 * S;
        x = (R(3,2) - R(2,3)) / S;
        y = (R(1,3) - R(3,1)) / S;
        z = (R(2,1) - R(1,2)) / S;
    elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2; % S = 4 * x
        w = (R(3,2) - R(2,3)) / S;
        x = 0.25 * S;
        y = (R(1,2) + R(2,1)) / S;
        z = (R(1,3) + R(3,1)) / S;
    elseif R(2,2) > R(3,3)
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2; % S = 4 * y
        w = (R(1,3) - R(3,1)) / S;
        x = (R(1,2) + R(2,1)) / S;
        y = 0.25 * S;
        z = (R(2,3) + R(3,2)) / S;
    else
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2; % S = 4 * z
        w = (R(2,1) - R(1,2)) / S;
        x = (R(1,3) + R(3,1)) / S;
        y = (R(2,3) + R(3,2)) / S;
        z = 0.25 * S;
    end
    quaternion = [x, y, z, w];
end

%% 欧拉角转四元数
function q = eulerToQuaternion(euler)
    % 输入: euler - [roll; pitch; yaw] 欧拉角（度）
    % 输出: q - 四元数 [qx, qy, qz, qw]
    
    % 将度转换为弧度
    roll = deg2rad(euler(1));
    pitch = deg2rad(euler(2));
    yaw = deg2rad(euler(3));
    
    % 计算半角
    cy = cos(yaw * 0.5);
    sy = sin(yaw * 0.5);
    cp = cos(pitch * 0.5);
    sp = sin(pitch * 0.5);
    cr = cos(roll * 0.5);
    sr = sin(roll * 0.5);

    % 计算四元数
    qw = cr * cp * cy + sr * sp * sy;
    qx = sr * cp * cy - cr * sp * sy;
    qy = cr * sp * cy + sr * cp * sy;
    qz = cr * cp * sy - sr * sp * cy;

    % 返回四元数
    q = [qx; qy; qz; qw];
end