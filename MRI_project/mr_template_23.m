% BME/EECS516
% MRI Project Template
clc 
clear all
% Other m-files required: ift2, ift, ft2, ft, blochsim_516
% Subfunctions: none
% MAT-files required: object18.mat

% Oct 2022; Last revision: Oct-30-2022

%% Select whether to load complex 2D object or create simple point object
complexobj = 0;
if complexobj
   % 2D Object for reconstruction
  load object23;
else
  % Single point object at (x,y,z) = (2,2,0) cm;
  % Point object has T1 of 1000 ms, T2 of 100 ms
  obj_x = 4;
  obj_y = 4;
  obj_z = 0;
  obj_T1 = 1000;
  obj_T2 = 100;
end

FOVx = 16;
FOVy = 12;
Nx = 80;
Ny = 60;
T_read = 8;
T_y = 2;
obj_n = length(obj_x); % Determine number of objects

%% Define simulation constants
% Physical constants
gambar = 42570;               % Gamma/2pi in kHz/T
gam = gambar*2*pi;            % Gamma in kiloradians/T

% Simulation values
dt = 0.05;                    % Time step for simulation, ms (50 us step size)
te = 10.0;                    % Echo time, ms
endtime = 16;                 % Total runtime of simulation, ms
time = [0:dt:endtime]';       % Vector containing each time step, ms (size #timepoints x 1)
totalTimepoints = length(time);          % Number of time points for simulation

% Initialize B vectors, the effective (x,y,z) applied magnetic field
% Vectors define applied magnetic field at time tp_n for object obj_n

bx = zeros([totalTimepoints obj_n]);
by = zeros([totalTimepoints obj_n]);
bz = zeros([totalTimepoints obj_n]);

% Define a 90 RF pulse
rf90pw = 3;                  % Pulse width in ms
sincper = rf90pw/4;          % in ms (this is the sinc stretch parameter)
                             % e.g. sinc(time/sincper) as shown below
rf_timepoints = rf90pw/dt;   % Number of simulation steps for RF
rf_time = [-(rf_timepoints-1)/2:(rf_timepoints-1)/2]'.*dt; % Time vector for creating sinc, centered at 0
rf_shape = hanning(rf_timepoints).*sinc(rf_time./sincper); % Sinc waveform shape with hanning window, with amplitude 1
rf_simulation = @(t)interp1(rf_time,rf_shape,t,'spline');
rf_amplitude90 = integral(rf_simulation,0,3);
rf_amplitude90 = pi/2/rf_amplitude90/gam;

%rf_amplitude90 = pi/2/gam/rf90pw         % REPLACE 0 with amplitude of the RF pulse here, in T

% Scale rf_shape by a_rf90 (amplitude), then fill the remainder of the time with zeros
b1_90 = rf_amplitude90.*[rf_shape; zeros([totalTimepoints-rf_timepoints 1])];
m0 = [0; 0; 1];
gz = 0;
omega_shift = gam * gz * 0;
%% part 1
bx = b1_90;
[mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1,obj_T2,dt);
figure (1)
subplot(3,1,1)
plot(time,mx);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
title('(x,y,z) = (4,4,0); (T1, T2) = (1000, 100)');

subplot(3,1,2)
plot(time,my);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,1,3)
plot(time,mz);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);

%% part 2
obj_T1_1 = 20;
obj_T2_2 = 10;
[mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1_1,obj_T2_2,dt);

figure (2)
subplot(3,1,1)
plot(time,mx);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
title('(x,y,z) = (4,4,0); (T1, T2) = (20, 10)');

subplot(3,1,2)
plot(time,my);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,1,3)
plot(time,mz);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);


%% Create gradients
% Create gz
rf90bw = 1 / sincper;    %bandwith of RF
slThick = 1;             % Slick thickness in cm
gz1_a = 2*pi*rf90bw/gam/slThick;                % REPLACE 0 with amplitude of gz1 in T/cm
gz1_pw = rf90pw;           % Match the width of gz1 to the RF pulse
gz2_a = -gz1_a;             % REPLACE 0 with amplitude of gz2 in T/cm
gz2_pw = rf90pw/2;
% Create gz with positve area gz1_a*gz1_pw, followed by negative area gz2_a*gz2pw
% gz step size is dt, with amplitude values in T/cm
gz =  (time < gz1_pw) .* gz1_a ...
       + (time >= gz1_pw).*(time < (gz1_pw+gz2_pw)) .* gz2_a;
bz_3 = gz * obj_z;

%% part 3
bz_3 = gz * obj_z;
[mx,my,mz] = blochsim_516(m0,bx,by,bz_3,obj_T1,obj_T2,dt);
figure (3)
subplot(3,1,1)
plot(time,mx);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
title('(x,y,z) = (4,4,0); (T1, T2) = (1000, 100); with Gz');

subplot(3,1,2)
plot(time,my);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,1,3)
plot(time,mz);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);

%% part 4
obj_z_41 = 0.2;
obj_z_42 = 1;

bz_41 = gz * obj_z_41;
bz_42 = gz * obj_z_42;

[mx_41,my_41,mz_41] = blochsim_516(m0,bx,by,bz_41,obj_T1,obj_T2,dt);
[mx_42,my_42,mz_42] = blochsim_516(m0,bx,by,bz_42,obj_T1,obj_T2,dt);
figure (4)
subplot(3,2,1)
plot(time,mx_41);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
subtitle('z = 0.2');

subplot(3,2,3)
plot(time,my_41);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,2,5)
plot(time,mz_41);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);

subplot(3,2,2)
plot(time,mx_42);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
subtitle('z = 1');

subplot(3,2,4)
plot(time,my_42);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,2,6)
plot(time,mz_42);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);

%% part 5

% Create gx
k_x = 1/FOVx;
gx_5b = 2 * pi * Nx * k_x / gam / T_read;
gx_5a = -2 * gx_5b;

t_gx_1 = 5.5;
t_gx_2 = 7.5;
t_gx_3 = 15.5;

gx_5 =  (time >= t_gx_1).*(time < t_gx_2) .* gx_5a ...
       + (time >= t_gx_2).*(time < t_gx_3) .* gx_5b;

bx_5 = b1_90;
by_5 = zeros([totalTimepoints obj_n]);
bz_5 = gx_5 * obj_x;

[mx,my,mz] = blochsim_516(m0,bx_5,by_5,bz_5,obj_T1,obj_T2,dt);
figure (5)
subplot(3,1,1)
plot(time,mx);
xlabel('time (ms)');
ylabel('Mx');
axis([0 endtime -1 1]);
title('(x,y,z) = (4,4,0); (T1, T2) = (1000, 100); with Gz, Gx');

subplot(3,1,2)
plot(time,my);
xlabel('time (ms)');
ylabel('My');
axis([0 endtime -1 1]);

subplot(3,1,3)
plot(time,mz);
xlabel('time (ms)');
ylabel('Mz');
axis([0 endtime -1 1]);

%% part 6
gx_func = @(t)interp1(time,gx_5,t,'linear');
gy = zeros([totalTimepoints obj_n]);
%gy_func = @(t)interp1(time,gx_5,t,'linear');
fi_r = zeros([totalTimepoints obj_n]);
s_r = zeros([totalTimepoints obj_n]);
m0_r = 1;
nread = 80;
npe = 60;
index_1 = 1;
sig = zeros([nread 1]);
for index = 1 : totalTimepoints
    fi_r(index,:) = integral(gx_func,0,(index-1) * dt);
    s_r(index,:) = m0_r * exp(-1i * 4 * gam * fi_r(index,:));
    s_r(index,:) = ift(s_r(index,:));
    if mod(index,2) == 1 && (index >= 151 && index <= 310)
        sig(index_1,:) = s_r(index,:);
        index_1 = index_1 + 1;
    end
end

fi_r = fi_r * obj_x * gam;
xpos = [-nread/2:nread/2-1]/nread.*FOVx;
ypos = [-npe/2:npe/2-1]/npe*FOVy;
figure(6)
plot(xpos,abs(ift(sig)));
%plot(fi_r ,s_r);
xlabel('x_position');
%ylabel('y_position');
ylabel('received signal');
%axis(xpos(1) xpos(79) -1 1);
title('magnitude vs spatial position');

%% part 7
% Create gy
k_y = 1/FOVy;
gy_max = pi * Ny * k_y / gam / T_y;
%gy_func = @(t)interp1(time,gy_7,t,'linear');
%fi_total = integral(gy_func,0,endtime);
%fi_total = fi_total * gam;
delta_gy = 2 * pi / (gam * 2 *FOVy);
% fi_x = zeros([totalTimepoints 1]);
% fi_y = zeros([npe 1]);
% s_r_7 = zeros([totalTimepoints npe]);
% sig_7 = zeros([nread npe]);
% index_1 = 1;
% for pe = 1:npe
%     gy_7 = (time >= t_gx_1).*(time < t_gx_2) .* (delta_gy * (npe - 1) - gy_max);
%     gy_func_7 = @(t)interp1(time, gy_7, t,'linear');
%     fi_y(pe,:) = obj_y * integral(gy_func_7,0,endtime);
%     for index = 1 : totalTimepoints
%         gx_func_7 = @(t)interp1(time,gx_5,t,'linear');
%         fi_x(index,:) = obj_x * integral(gx_func_7,0,(index-1) * dt);
%         s_r_7(index,pe) = m0_r * exp(-1i * gam * (fi_x(index,:) + fi_y(pe,:)));
%         s_r_7(index,pe) = ift2(s_r_7(index,pe));
%         if (index >= 151 && index <= 310) && mod(index,2) == 1 && index_1 <= nread
%             sig_7(index_1,pe) = s_r_7(index,pe);
%             index_1 = index_1 + 1;
%         end
%     end
% end
by_7 = zeros([totalTimepoints obj_n]);
sig_7 = zeros([nread npe]);
index_1 = 1;
M_7 = zeros([totalTimepoints npe]);
for pe = 1:npe
    gy_7 = (time >= t_gx_1).*(time < t_gx_2) .* (delta_gy * (pe - 1) - gy_max);
    bz_7 = gx_5 * obj_x + gy_7 * obj_y + gz * obj_z;
    [mx_7,my_7,mz_7] = blochsim_516(m0,bx_5,by_7,bz_7,obj_T1,obj_T2,dt);
    M_7(:,pe)= mx_7 + 1i * my_7;
    for index = 1 : totalTimepoints
        %M_7(index,pe) = ift2(M_7(index,pe));
        if (index >= 151 && index <= 310) && mod(index,2) == 1 && index_1 <= nread
        sig_7(index_1,pe) = M_7(index,pe);
        index_1 = index_1 + 1;
        end
    end
end


% show images for parts 7-10
delta_kx = gambar * gx_5b * dt;
delta_ky = gambar * delta_gy * T_y;
W_kx = delta_kx * Nx;
W_ky = delta_ky * Ny;
kxpos = linspace(-W_kx/2, W_kx/2, Nx);  % vector of kx locations
kypos = linspace(-W_ky/2, W_ky/2, Ny);  % vector of ky locations

figure(7)
imagesc(kxpos,kypos,real(sig_7)); colormap gray; axis('image'); axis('xy')
xlabel('x (units?)');
ylabel('y (units?)');
title('REAL M(kx,ky)')

disp 'Press any key to continue...'; pause

figure(7)
imagesc(kxpos,kypos,imag(sig_7)); colormap gray; axis('image'); axis('xy')
xlabel('x (units?)');
ylabel('y (units?)');
title('IMAG M(kx,ky)')
disp 'Press any key to continue...'; pause

%% part 8
for pe = 1:npe
    for index = 1 : nread
        sig_7(index,pe) = ift2(sig_7(index,pe));
    end
end

figure(8)
imagesc(xpos,ypos,abs(sig_7)); colormap gray; axis('image'); axis('xy')
xlabel('x (units?)');
ylabel('y (units?)');
title('abs(image(x,y))')

disp 'Press any key to continue...'; pause

%% part 9
by_9 = zeros([totalTimepoints obj_n]);
sig_9 = zeros([nread npe]);
sig_9 = zeros([nread npe]);
index_1 = 1;
M_9 = zeros([totalTimepoints npe]);
for pe = 1:npe
    gy_9 = (time >= t_gx_1).*(time < t_gx_2) .* (delta_gy * (npe - 1) - gy_max);
    bz_9 = gx_5 * obj_x + gy_9 * (obj_y + 10) + gz * obj_z;
    [mx_9,my_9,mz_9] = blochsim_516(m0,bx_5,by_7,bz_7,obj_T1,obj_T2,dt);
    M_9(:,pe)= mx_9 + 1i * my_9;
    for index = 1 : totalTimepoints
        M_9(index,pe) = ift2(M_9(index,pe));
        if (index >= 151 && index <= 310) && mod(index,2) == 1 && index_1 <= nread
        sig_9(index_1,pe) = M_9(index,pe);
        index_1 = index_1 + 1;
        end
    end
end
figure(9)
imagesc(xpos,ypos,abs(sig_9)); colormap gray; axis('image'); axis('xy')
xlabel('x (units?)');
ylabel('y (units?)');
title('abs(image(x,y))')

disp 'Press any key to continue...'; pause

%% part 10
load object23.mat
mx_10sum = zeros([totalTimepoints 1]);
my_10sum = zeros([totalTimepoints 1]);
mz_10sum = zeros([totalTimepoints 1]);
obj_n = 2850;
m0_10 = [zeros([2 obj_n]);
         ones([1 obj_n])];
M_10 = zeros([totalTimepoints npe]);
sig_10 = zeros([totalTimepoints npe]);
% for slice = 1:2  % slice loop
    %omega_shift = gam * gz1_a * slice;
    for pe = 1:npe
        bx_10 = b1_90 * ones([1 obj_n]);
        by_10 = b1_90 * 0 * ones([1 obj_n]);
        gx_10 = gx_5;
        gy_10 = (time >= t_gx_1).*(time < t_gx_2) .* (delta_gy * (npe - 1) - gy_max);
        gz_10 = gz;
        bz_10 = gx_10 * obj_x + gy_10 * obj_y + gz * obj_z;
        [mx_10,my_10,mz_10] = blochsim_516(m0_10,bx_10,by_10,bz_10,obj_T1,obj_T2,dt);
        for num = 1 : obj_n
            mx_10sum = mx_10sum + mx_10(:,obj_n);
            my_10sum = my_10sum + my_10(:,obj_n);
            mz_10sum = mz_10sum + mz_10(:,obj_n);
        end
        M_10(:,pe)= mx_10sum + 1i * my_10sum;
        for index = 1 : totalTimepoints
            M_10(index,pe) = ift2(M_10(index,pe));
            if (index >= 151 && index <= 310) && mod(index,2) == 1 && index_1 <= nread
            sig_10(index_1,pe) = M_10(index,pe);
            index_1 = index_1 + 1;
            end
        end
    end
%end

figure(10)
imagesc(xpos,ypos,abs(sig_10)); colormap gray; axis('image'); axis('xy')
xlabel('x (units?)');
ylabel('y (units?)');
title('abs(image(x,y))')

disp 'Press any key to continue...'; pause

imagesc(kxpos,kypos,abs(M_10)); colormap gray; axis('image'); axis('xy')
xlabel('kx (units?)');
ylabel('ky (units?)');
title('abs(M(kx,ky))')


% %% Loop through s lices and phase encodes
% for slice = 1:2  % slice loop
%     
%     modify B1 here for other slices 
%     
%     bx = b1_90 * sin(omega_shift * ) ones([1 obj_n]);




% end  % end of slice loop

