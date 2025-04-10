clc;
clear;
% ���غͲ鿴����
load('data\imu_data.mat');
load('data\ground_truth_data.mat');
load('data\leg_data.mat');
%% ������
%��timeseriesתΪ���� 
a_R = [-1,0,0;0,-1,0;0,0,1];
% a_R = eye(3);
acce_ = transpose(imu_data{:, 2:4}); %�߼��ٶ�
gyro_ = transpose(imu_data{:, 5:7}); %���ٶ�
position_true = transpose(ground_truth_data{:, 2:4}); %ʵ��λ��
R_true =  transpose(ground_truth_data{:, 5:8}); %ʵ����̬
p_ = transpose((ground_truth_data{1, 2:4}))+[0;0;0.34];   %��ʼλ��
v_ = [0; 0; 0];   %��ʼ�ٶ�
R_ = quaternionToRotationMatrix(ground_truth_data{1, 5:8}); %��ó�ʼ��̬����
ba_ = [1; 1; 1] * 3.552e-04 * 0;   %�߼��ٶ�ƫ��
bg_ = [1; 1; 1] * 1.655e-05 * 0;   %���ٶ�ƫ��
gravity_ = [0; 0; -9.8]; %�������ٶ�  
% ȷ��������ʱ������
dt = 1/765;      % ����ʱ����
end_time = 52580*dt;    %����ʱ��Ϊ8s
time = 0:dt:end_time; %�ӵ�1�����ݿ�ʼ
% ȷ����ֵ����ʱ��
dt_true = 1/63;      % ����ʱ����
end_time = 4400*dt_true;    %����ʱ��Ϊ8s
time_true = 0:dt_true:end_time; %�ӵ�1�����ݿ�ʼ
% ��ʼ�������Դ洢ÿ��ʱ�䲽�Ľ��
num_steps = length(time);
num_steps_true = length(time_true);
p_history = zeros(3, num_steps);
v_history = zeros(3, num_steps);
R_history = zeros(3, 3, num_steps);
o_history = zeros(3, num_steps); 
o_true = zeros(3, num_steps_true); % ��ʼ�����ڴ洢�Աȵ���ʵŷ��������
position_true_resampled = position_true(:,1:num_steps_true); %������ʵλ��
% velocity_true_resampled = velocity_true(:,0:num_steps_true); %������ʵ�ٶ�
R_true_resampled = R_true(:,1:num_steps_true); %������ʵ��̬
q_history = zeros(4, num_steps); %��Ԫ��

%% �����Ȳ�����
s_p0_array = zeros(3, num_steps); %��ʼ���Ȳ����������λ��
s_p1_array = zeros(3, num_steps);
s_p2_array = zeros(3, num_steps);
s_p3_array = zeros(3, num_steps);

d0 = [0;0;0]; %����inertial frame�µ�λ��
d1 = [0;0;0];
d2 = [0;0;0];
d3 = [0;0;0];

s_p0_speed = transpose(leg_data{1:num_steps, 5:7});   %��0�ٶ�(body)
s_p1_speed = transpose(leg_data{1:num_steps, 15:17}); %��1�ٶ�
s_p2_speed = transpose(leg_data{1:num_steps, 25:27}); %��2�ٶ�
s_p3_speed = transpose(leg_data{1:num_steps, 35:37}); %��3�ٶ�

d0_history = transpose(leg_data{1:num_steps, 2:4});   %��0λ��(body)
d1_history = transpose(leg_data{1:num_steps, 12:14}); %��1λ��
d2_history = transpose(leg_data{1:num_steps, 22:24}); %��2λ��
d3_history = transpose(leg_data{1:num_steps, 32:34}); %��3λ��

contact_0 = transpose(leg_data{1:num_steps, 11}); %��0�Ӵ��ж�
contact_1 = transpose(leg_data{1:num_steps, 21}); %��1�Ӵ��ж�
contact_2 = transpose(leg_data{1:num_steps, 31}); %��2�Ӵ��ж�
contact_3 = transpose(leg_data{1:num_steps, 41}); %��3�Ӵ��ж�

v_contact_history = zeros(3, num_steps); %�ȵ���ʷ�ٶ�
d_history = zeros(3, num_steps); 
%% ��ʼ���Ŀ������˲�����
A = [eye(3),eye(3)*dt,zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),eye(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)];     % ״̬ת�ƾ���
R = [eye(3)*0.01,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),eye(3)*0.01,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3),zeros(3, 4);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.1,zeros(3, 4);
     zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),zeros(4,3),eye(4)*0.001];     % �۲�����Э������� %25
Q = [eye(3)*0.00001,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),eye(3)*0.00001,zeros(3),zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3);
     zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01];     % ״̬ת������Э�������
P = [ eye(3)*0.001,zeros(3),zeros(3),zeros(3),zeros(3),zeros(3);
      zeros(3),eye(3)*0.001,zeros(3),zeros(3),zeros(3),zeros(3);
      zeros(3),zeros(3),eye(3)*0.001,zeros(3),zeros(3),zeros(3);
      zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3),zeros(3);
      zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.01,zeros(3);
      zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),eye(3)*0.001];     % ״̬ת������Э�������
 
big_noise = 100; %̧�ȴ�����
normal_noise = 0.001; %��������
%% ��*������
for i = 1:num_steps_true
    o_true(:, i) =  quaternionToEuler(transpose( R_true_resampled(:,i) ));
end

for i = 1:num_steps
     rotationMatrix = eulerAnglesToRotationMatrix((gyro_(1:3, i) - bg_) * dt);  % ŷ����ת��ת����
     R_ = R_ * rotationMatrix;% ����ŷ����˳��Ϊ ZYX
%     R_ = quaternionToRotationMatrix(transpose( R_true_resampled(:,i) )); %����ֵ����
     p_ = p_ + v_ * dt + 0.5 * gravity_ * dt * dt + 0.5 * R_ * (acce_(1:3, i) - ba_) * dt * dt; %����λ��
     v_ = v_ + R_ * (acce_(1:3, i) - ba_) * dt + gravity_ * dt; % �����ٶ�

     %�۲����(�Ӵ���)
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
          0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,contact_3(i);];     % �۲����
    
%      %�۲����(������)
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
%           0,0,0,   0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,contact_3(i);];     % �۲����
     
     %�۲�ֵ
     delta_d0_world = R_*d0_history(:,i); %��������ϵ���Ȳ���λ�ò�
     delta_d1_world = R_*d1_history(:,i);
     delta_d2_world = R_*d2_history(:,i);
     delta_d3_world = R_*d3_history(:,i);
     Z = [delta_d0_world;  %��������ϵ���Ȳ���λ�ò�
          delta_d1_world;
          delta_d2_world;
          delta_d3_world;
          -R_*(cross((gyro_(1:3, i) - bg_),d0_history(:,i)) + s_p0_speed(:, i)); %��������ϵ�µ��ٶ�
          -R_*(cross((gyro_(1:3, i) - bg_),d1_history(:,i)) + s_p1_speed(:, i));
          -R_*(cross((gyro_(1:3, i) - bg_),d2_history(:,i)) + s_p2_speed(:, i));
          -R_*(cross((gyro_(1:3, i) - bg_),d3_history(:,i)) + s_p3_speed(:, i));
          0;   %��������ϵ���Ȳ���z����λ��
          0;
          0;
          0]; 
           
      
%       %����̧�ȵ�������
%        contacts = [contact_0(i), contact_1(i), contact_2(i), contact_3(i)];
%        [Q, R] = adjustNoise(Q, R, contacts, big_noise, normal_noise);
   
      % �������˲�
      % Ԥ�ⲽ��
      x_hat_minus = [p_;v_;d0;d1;d2;d3];
      P_minus = A * P * A' + Q;
      % ���²���
      K = P_minus * H' / (H * P_minus * H' + R);
      x_hat = x_hat_minus + K * (Z - H * x_hat_minus);
%       x_hat = x_hat_minus ;
      P = (eye(18) - K * H) * P_minus;
      % �洢�˲����״̬
      p_history(:, i) = x_hat(1:3,1);
      v_history(:, i) = x_hat(4:6,1);
      R_history(:, :, i) = R_;
      
      %��RתΪŷ����
     o_history(:, i) = a_R * rotMatrixToWorldEulerAngles(R_);
     q_history(:, i) = eulerToQuaternion(o_history(:, i)); %��¼ÿһ�ε���Ԫ��
%      o_true(:, i) =  quaternionToEuler(transpose( R_true_resampled(:,i) ));
     %��������
     p_ = x_hat(1:3,1);
     v_ = x_hat(4:6,1);
     d0 = x_hat(7:9,1);
     d1 = x_hat(10:12,1);
     d2 = x_hat(13:15,1);
     d3 = x_hat(16:18,1);
     
     % �洢���
     p_history(:, i) = a_R*(p_+[0;0;0.026]);
     v_history(:, i) = v_;
     R_history(:, :, i) = R_;
     d_history(:, i) = d0;
end

%% ��Ҫ��ӡ������
ground_truth_data{1:num_steps_true,1} = transpose(time_true);
outcome_ground_truth = ground_truth_data{1:num_steps_true,1:8};
outcome_kf = transpose([time;p_history(:,:);q_history(:,:)]);
% �� outcome_luen_luen ����Ϊ TXT �ļ�
writematrix(outcome_kf, 'KF', 'Delimiter', ' ');
% �� outcome_ground_truth ����Ϊ TXT �ļ�
writematrix(outcome_ground_truth, 'Ground_Truth.txt', 'Delimiter', ' ');

% %% ��ӡ�������
% % ���������MSE��
% mse_x = mean((p_history(1, :) - position_true_resampled(1, :)).^2);
% mse_y = mean((p_history(2, :) - position_true_resampled(2, :)).^2);
% mse_z = mean((p_history(3, :) - position_true_resampled(3, :)).^2);
% 
% % �����������RMSE��
% rmse_x = sqrt(mse_x);
% rmse_y = sqrt(mse_y);
% rmse_z = sqrt(mse_z);
% 
% % ����ƽ��������MAE��
% mae_x = mean(abs(p_history(1, :) - position_true_resampled(1, :)));
% mae_y = mean(abs(p_history(2, :) - position_true_resampled(2, :)));
% mae_z = mean(abs(p_history(3, :) - position_true_resampled(3, :)));
% 
% % ��ʾ���
% fprintf('X�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_x, rmse_x, mae_x);
% fprintf('Y�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_y, rmse_y, mae_y);
% fprintf('Z�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_z, rmse_z, mae_z);

% %% ��ӡ��ŷ���Ƿ����
% % ���������MSE��
% mse_x = mean((o_history(3, :) - o_true(1, :)).^2);
% mse_y = mean((o_history(2, :) - o_true(2, :)).^2);
% mse_z = mean((o_history(1, :) - o_true(3, :)).^2);
% 
% % �����������RMSE��
% rmse_x = sqrt(mse_x);
% rmse_y = sqrt(mse_y);
% rmse_z = sqrt(mse_z);
% 
% % ����ƽ��������MAE��
% mae_x = mean(abs(o_history(3, :) - o_true(1, :)));
% mae_y = mean(abs(o_history(2, :) - o_true(2, :)));
% mae_z = mean(abs(o_history(1, :) - o_true(3, :)));
% 
% % ��ʾ���
% fprintf('Roll�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_x, rmse_x, mae_x);
% fprintf('Pitch�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_y, rmse_y, mae_y);
% fprintf('Yaw�����ϵ�MSE: %.4f, RMSE: %.4f, MAE: %.4f\n', mse_z, rmse_z, mae_z);

% ����λ��ͼ
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

% ������̬ͼ
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

% ���ƶ�ά�켣ͼ
figure;

% ����ÿ���߲��ֱ���������
plot(position_true_resampled(1, :), position_true_resampled(2, :), 'r-', 'LineWidth', 1.5); % ��ɫ��ʵ��
hold on;
plot(p_history(1, :), p_history(2, :), 'b', 'LineWidth', 1.5); % ��ɫ������

% ��ӱ����ͼ��
title('2D Trajectory - XY Plane');
xlabel('X Position');
ylabel('Y Position');
legend('True', 'Computed');
% �ͷ� hold
hold off;

%% ������������
    % ����:
    %   Q - ״̬ת������Э�������
    %   R - �۲�����Э�������
    %   contacts - �Ӵ�״̬����
    %   big_noise - ������ֵ
    %   normal_noise - ��������ֵ
function [Q, R] = adjustNoise(Q, R, contacts, big_noise, normal_noise)
    % ������������
    % ����:
    %   Q - ״̬ת������Э�������
    %   R - �۲�����Э�������
    %   contacts - �Ӵ�״̬����
    %   big_noise - ������ֵ
    %   normal_noise - ��������ֵ
    % ���:
    %   Q - �������״̬ת������Э�������
    %   R - ������Ĺ۲�����Э�������
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
%% ŷ���Ǻ���
function R = eulerAnglesToRotationMatrix(angles)
    % angles ��һ��������X�ᡢY���Z���ŷ���ǵ� 1x3 ����

    % ��ȡŷ����
    roll = angles(1); %��x����ת
    pitch = angles(2);   %��y����ת
    yaw = angles(3);  %��z����ת

    % ������X�����ת����
    Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];

    % ������Y�����ת����
    Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
    
    % ������Z�����ת����
    Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];

    % �������յ���ת����
    R = Rz * Ry * Rx;
end

%% ��ת����תŷ����
function euler_angles = rotMatrixToWorldEulerAngles(rotation_matrix)
    % rotation_matrix ��һ��3x3����ת����

    % ����ת����תΪŷ���ǣ�ʹ�� 'ZYX' ��ת˳��
    euler_angles = transpose(rotm2eul(rotation_matrix, 'ZYX'));
    
    %�ѻ�����roll pitch yaw ��˳�򷢲�
    yaw = euler_angles(1);
    pitch = euler_angles(2);
    roll = euler_angles(3);
    euler_angles = [roll;pitch;yaw];
    
    % ������תΪ�Ƕ�
    euler_angles = rad2deg(euler_angles);
    
end

%% ��Ԫ��תŷ����
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
  
    % ������תΪ�Ƕ�
    gyro = rad2deg(gyro);
end

%% ��Ԫ��תR����
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

%% R����ת��Ԫ��
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

%% ŷ����ת��Ԫ��
function q = eulerToQuaternion(euler)
    % ����: euler - [roll; pitch; yaw] ŷ���ǣ��ȣ�
    % ���: q - ��Ԫ�� [qx, qy, qz, qw]
    
    % ����ת��Ϊ����
    roll = deg2rad(euler(1));
    pitch = deg2rad(euler(2));
    yaw = deg2rad(euler(3));
    
    % ������
    cy = cos(yaw * 0.5);
    sy = sin(yaw * 0.5);
    cp = cos(pitch * 0.5);
    sp = sin(pitch * 0.5);
    cr = cos(roll * 0.5);
    sr = sin(roll * 0.5);

    % ������Ԫ��
    qw = cr * cp * cy + sr * sp * sy;
    qx = sr * cp * cy - cr * sp * sy;
    qy = cr * sp * cy + sr * cp * sy;
    qz = cr * cp * sy - sr * sp * cy;

    % ������Ԫ��
    q = [qx; qy; qz; qw];
end