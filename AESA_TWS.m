clear;
clc;
close all;

%--------------------------------------------------------------------------
% CONSTANTS
%--------------------------------------------------------------------------
c = 3e8;
%--------------------------------------------------------------------------
% RADAR PARAMETERS
%--------------------------------------------------------------------------
Pt = 50e3;
fc = 3e9;
lambda = c/fc;
Tp = 2e-6;
B = 5e6;
Gt = 1;
Gr = 1;
L_dB = 6;
L = 10^(L_dB/10);
sigma = 1;
PRI = 1e-3;
PRF = 1/PRI;
Npulses = 32;
fs = 20e6;
%--------------------------------------------------------------------------
% URA PARAMETERS
%--------------------------------------------------------------------------
Mx = 8;
My = 8;
M = Mx*My;
dx = lambda/2;
dy = lambda/2;
%--------------------------------------------------------------------------
% URA GEOMETRY
%--------------------------------------------------------------------------
[mx_grid,my_grid] = meshgrid(0:Mx-1,0:My-1);
x_pos = mx_grid*dx;
y_pos = my_grid*dy;
x_vec = x_pos(:);
y_vec = y_pos(:);
%--------------------------------------------------------------------------
% FIELD OF VIEW
%--------------------------------------------------------------------------
az_min = 30;
az_max = 150;
el_min = 30;
el_max = 150;
az_min_deg = rad2deg(az_min);
az_max_deg = rad2deg(az_max);
el_min_deg = rad2deg(el_min);
el_max_deg = rad2deg(el_max);
%--------------------------------------------------------------------------
% TARGET
%--------------------------------------------------------------------------
R0 = 68e3;
v0 = 17;
tau = 2*R0/c;
fd = 2*v0/lambda;
az_target = 85;
el_target = 85;
%--------------------------------------------------------------------------
% TARGET ANGLES
%--------------------------------------------------------------------------
az0_target = deg2rad(az_target);
el0_target = deg2rad(el_target);
ux_target = cos(el0_target)*cos(az0_target);
uy_target = cos(el0_target)*sin(az0_target);
k0 = 2*pi/lambda;
steer_target = exp(-1j*k0*(x_vec*ux_target + y_vec*uy_target));
%--------------------------------------------------------------------------
% BEAM STEERING ANGLES
%--------------------------------------------------------------------------
az_beam = 90;
el_beam = 90;
az0_beam = deg2rad(az_beam);
el0_beam = deg2rad(el_beam);
ux_beam = cos(el0_beam)*cos(az0_beam);
uy_beam = cos(el0_beam)*sin(az0_beam);
steer_beam = exp(-1j*k0*(x_vec*ux_beam + y_vec*uy_beam));
%--------------------------------------------------------------------------
% TX BEAM RESPONSE TOWARD TARGET
%--------------------------------------------------------------------------
AF_target = sum(conj(steer_beam).*steer_target)/M;
Gt_eff = Gt*abs(AF_target)^2;
Gr_eff = Gr*abs(AF_target)^2;
%--------------------------------------------------------------------------
% RADAR EQUATION
%--------------------------------------------------------------------------
Pr = Pt*Gt_eff*Gr_eff*lambda^2*sigma/((4*pi)^3*R0^4*L);
Ar = sqrt(Pr);
beam_loss_dB_one_way = 20*log10(abs(AF_target)+eps);
beam_loss_dB_two_way_power = 40*log10(abs(AF_target)+eps);
disp(['One-way beam field loss: ' num2str(beam_loss_dB_one_way) ' dB'])
disp(['Two-way beam power loss: ' num2str(beam_loss_dB_two_way_power) ' dB'])
%--------------------------------------------------------------------------
% SNR
%--------------------------------------------------------------------------
SNR_dB = 0;
noise_var = Ar^2/(10^(SNR_dB/10));
%--------------------------------------------------------------------------
% CLUTTER CONFIGURATION
%--------------------------------------------------------------------------
clutter_mode = 'deterministic';
SCR_dB = SNR_dB;
clutter_power = Ar^2/(10^(SCR_dB/10));
%--------------------------------------------------------------------------
% TIME AXES
%--------------------------------------------------------------------------
Nfast = round(PRI*fs);
Npulse = round(Tp*fs);
t = (0:Npulse-1)/fs;
delay = round(tau*fs);
%--------------------------------------------------------------------------
% TRANSMITTED CHIRP
%--------------------------------------------------------------------------
k = B/Tp;
tx = exp(1j*pi*k*t.^2);
Es = sum(abs(tx).^2);
%--------------------------------------------------------------------------
% MATCHED FILTER
%--------------------------------------------------------------------------
mf = conj(fliplr(tx));
%--------------------------------------------------------------------------
% STORAGE
%--------------------------------------------------------------------------
rx = zeros(M,Npulses,Nfast);
mf_out = zeros(M,Npulses,Nfast+Npulse-1);
RD_cube = zeros(M,Npulses,Nfast+Npulse-1);
%--------------------------------------------------------------------------
% CLUTTER MAP
%--------------------------------------------------------------------------
Nclutter = 50;
R_clutter_min = 20e3;
R_clutter_max = 80e3;
Rc0 = zeros(Nclutter,1);
az0_clutter = zeros(Nclutter,1);
el0_clutter = zeros(Nclutter,1);
Ac0 = zeros(Nclutter,1);
vc0 = zeros(Nclutter,1);
for kc = 1:Nclutter
    Rc0(kc) = R_clutter_min + (R_clutter_max-R_clutter_min)*rand;
    az0_clutter(kc) = az_min + (az_max-az_min)*rand;
    el0_clutter(kc) = el_min + (el_max-el_min)*rand;
    Ac0(kc) = sqrt(clutter_power/Nclutter);
    vc0(kc) = -2 + 4*rand;
end
%--------------------------------------------------------------------------
% MAIN LOOP
%--------------------------------------------------------------------------
for m = 1:M
    for p = 1:Npulses
        %------------------------------------------------------------------
        % TARGET ECHO
        %------------------------------------------------------------------
        echo = zeros(1,Nfast);
        t_echo = (0:Npulse-1)/fs + (p-1)*PRI;
        echo(delay+1:delay+Npulse) = Ar*steer_target(m)*tx.*exp(1j*2*pi*fd*t_echo);
        %------------------------------------------------------------------
        % AWGN
        %------------------------------------------------------------------
        noise = sqrt(noise_var/2)*(randn(1,Nfast)+1j*randn(1,Nfast));
        %------------------------------------------------------------------
        % CLUTTER
        %------------------------------------------------------------------
        clutter = zeros(1,Nfast);
        for kc = 1:Nclutter
            Rc = Rc0(kc) + 5*randn;
            az_c = az0_clutter(kc) + deg2rad(0.1)*randn;
            el_c = el0_clutter(kc) + deg2rad(0.1)*randn;
            vc = vc0(kc) + 0.05*randn;
            switch clutter_mode
                case 'rayleigh'
                    Ac = sqrt(clutter_power/(2*Nclutter))*sqrt(randn^2 + randn^2);
                case 'weibull'
                    k_weib = 1.5;
                    lambda_weib = sqrt((clutter_power/Nclutter)/gamma(1+2/k_weib));
                    u = rand;
                    Ac = lambda_weib*(-log(1-u))^(1/k_weib);
                case 'deterministic'
                    Ac = Ac0(kc)*(1+0.005*randn);
            end
            tau_c = 2*Rc/c;
            delay_c = round(tau_c*fs);
            fd_c = 2*vc/lambda;
            t_clutter = (0:Npulse-1)/fs + (p-1)*PRI;
            ux_c = cos(el_c)*cos(az_c);
            uy_c = cos(el_c)*sin(az_c);
            steer_c = exp(-1j*(2*pi/lambda)*(x_vec(m)*ux_c + y_vec(m)*uy_c));
            if delay_c+Npulse <= Nfast
                clutter(delay_c+1:delay_c+Npulse) = clutter(delay_c+1:delay_c+Npulse) + Ac*steer_c*tx.*exp(1j*2*pi*fd_c*t_clutter);
            end
        end
        %------------------------------------------------------------------
        % RECEIVED SIGNAL
        %------------------------------------------------------------------
        rx(m,p,:) = echo + clutter + noise;
        %------------------------------------------------------------------
        % MATCHED FILTER
        %------------------------------------------------------------------
        mf_out(m,p,:) = conv(squeeze(rx(m,p,:)).',mf);
    end
end
%--------------------------------------------------------------------------
% SENSOR-DOMAIN RANGE-DOPPLER CUBE
%--------------------------------------------------------------------------
RD_cube = fftshift(fft(mf_out,Npulses,2),2);
%--------------------------------------------------------------------------
% TX-STEERED QUASI-COHERENT SENSOR INTEGRATION
%--------------------------------------------------------------------------
RD_beam = zeros(Npulses,size(RD_cube,3));
for m = 1:M
    RD_beam = RD_beam + squeeze(RD_cube(m,:,:))*conj(steer_beam(m));
end
RD_beam_power = abs(RD_beam).^2;
RD_beam_dB = 10*log10(RD_beam_power/max(RD_beam_power(:)) + eps);
%--------------------------------------------------------------------------
% RANGE-DOPPLER AXES
%--------------------------------------------------------------------------
range_axis = (0:size(RD_beam,2)-1)*c/(2*fs);
fd_axis = (-Npulses/2:Npulses/2-1)*(PRF/Npulses);
vel_axis = fd_axis*lambda/2;
%--------------------------------------------------------------------------
% RANGE-DOPPLER ESTIMATION
%--------------------------------------------------------------------------
[~,idx] = max(RD_beam_power(:));
[dop_bin,range_bin] = ind2sub(size(RD_beam_power),idx);
R_est = range_axis(range_bin);
v_est = vel_axis(dop_bin);
%--------------------------------------------------------------------------
% DISPLAY ESTIMATES
%--------------------------------------------------------------------------
disp('--------------------------------------------')
disp(['True range:        ' num2str(R0) ' m'])
disp(['Estimated range:   ' num2str(R_est) ' m'])
disp(['True velocity:     ' num2str(v0) ' m/s'])
disp(['Estimated velocity:' num2str(v_est) ' m/s'])
disp('--------------------------------------------')
%--------------------------------------------------------------------------
% BEAMFORMED RANGE-DOPPLER MAP
%--------------------------------------------------------------------------
figure;
imagesc(range_axis/1000,vel_axis,RD_beam_dB)
axis xy
xlabel('Range (km)')
ylabel('Radial Velocity (m/s)')
title(['Beamformed Range-Doppler Map: beam = ' num2str(az_beam) '°, ' num2str(el_beam) '°; target = ' num2str(az_target) '°, ' num2str(el_target) '°'])
colorbar
hold on
plot(R_est/1000,v_est,'rs','MarkerSize',18,'LineWidth',2)
plot(R0/1000,v0,'wo','MarkerSize',12,'LineWidth',2)
legend('Estimated Target','True Target')
%
%--------------------------------------------------------------------------
% ARRAY GEOMETRY
%--------------------------------------------------------------------------
figure
plot3(x_vec/lambda,y_vec/lambda,zeros(size(x_vec)),'ko','MarkerFaceColor','k','MarkerSize',6)
hold on
quiver3(0,0,0,5,0,0,'k','LineWidth',2)
quiver3(0,0,0,0,5,0,'k','LineWidth',2)
quiver3(0,0,0,0,0,5,'k','LineWidth',2)
text(5.2,0,0,'x')
text(0,5.2,0,'y')
text(0,0,5.2,'z / broadside')
grid on
axis equal
xlabel('x/\lambda')
ylabel('y/\lambda')
zlabel('z/\lambda')
title('8x8 URA Geometry')
view(40,30)
%
%--------------------------------------------------------------------------
% BEAM PATTERN IN DIRECTION-COSINE SPACE
%--------------------------------------------------------------------------
FOV_broadside = deg2rad(60);
u_vec = linspace(-sin(FOV_broadside),sin(FOV_broadside),201);
v_vec = linspace(-sin(FOV_broadside),sin(FOV_broadside),201);
[U,V] = meshgrid(u_vec,v_vec);
visible = U.^2 + V.^2 <= sin(FOV_broadside)^2;
AF = nan(size(U));

for ia = 1:size(U,1)
    for ja = 1:size(U,2)
        if visible(ia,ja)
            steer_scan = exp(-1j*k0*(x_vec*U(ia,ja) + y_vec*V(ia,ja)));
            AF(ia,ja) = abs(steer_beam'*steer_scan)/M;
        end
    end
end

AF = AF/max(AF(:),[],'omitnan');
W = sqrt(max(0,1-U.^2-V.^2));
X = AF.*U;
Y = AF.*V;
Z = AF.*W;

az0_beam = deg2rad(az_beam);
el0_beam = deg2rad(el_beam);
ux_beam_plot = cos(el0_beam)*cos(az0_beam);
uy_beam_plot = cos(el0_beam)*sin(az0_beam);
uz_beam_plot = sin(el0_beam);

az0_target = deg2rad(az_target);
el0_target = deg2rad(el_target);
ux_target_plot = cos(el0_target)*cos(az0_target);
uy_target_plot = cos(el0_target)*sin(az0_target);
uz_target_plot = sin(el0_target);

figure
surf(X,Y,Z,20*log10(AF+eps),'EdgeColor','none','FaceAlpha',0.45)
hold on
quiver3(0,0,0,1.2,0,0,'k','LineWidth',2)
quiver3(0,0,0,0,1.2,0,'k','LineWidth',2)
quiver3(0,0,0,0,0,1.2,'k','LineWidth',2)
quiver3(0,0,0,ux_beam_plot,uy_beam_plot,uz_beam_plot,'b','LineWidth',3,'MaxHeadSize',0.5)
quiver3(0,0,0,ux_target_plot,uy_target_plot,uz_target_plot,'r','LineWidth',3,'MaxHeadSize',0.5)
text(1.3,0,0,'x')
text(0,1.3,0,'y')
text(0,0,1.3,'z / broadside')
text(1.05*ux_beam_plot,1.05*uy_beam_plot,1.05*uz_beam_plot,sprintf('Beam (%d°, %d°)',az_beam,el_beam),'Color','b','FontWeight','bold')
text(1.05*ux_target_plot,1.05*uy_target_plot,1.05*uz_target_plot,sprintf('Target (%d°, %d°)',az_target,el_target),'Color','r','FontWeight','bold')
grid on
axis equal
xlabel('u_x')
ylabel('u_y')
zlabel('u_z / broadside')
title('8x8 URA Beam Pattern in Broadside-Centered FoV')
colorbar
clim([-40 0])
view(45,25)
%
% %--------------------------------------------------------------------------
% % CA-CFAR PARAMETERS
% %--------------------------------------------------------------------------
% Ntrain_r = 20;
% Ntrain_d = 2;
% Nguard_r = 4;
% Nguard_d = 1;
% Pfa = 1e-6;
% %--------------------------------------------------------------------------
% % CA-CFAR INITIALIZATION
% %--------------------------------------------------------------------------
% [Ndop,Nrange] = size(RD_beam_power);
% CFAR_map = zeros(Ndop,Nrange);
% CFAR_threshold = zeros(Ndop,Nrange);
% Ntrain_total = (2*(Ntrain_d+Nguard_d)+1)*(2*(Ntrain_r+Nguard_r)+1) - (2*Nguard_d+1)*(2*Nguard_r+1);
% alpha_cfar = Ntrain_total*(Pfa^(-1/Ntrain_total)-1);
% %--------------------------------------------------------------------------
% % 2D CA-CFAR OVER RANGE-DOPPLER MAP
% %--------------------------------------------------------------------------
% for d = Ntrain_d+Nguard_d+1:Ndop-(Ntrain_d+Nguard_d)
%     for r = Ntrain_r+Nguard_r+1:Nrange-(Ntrain_r+Nguard_r)
%         d_min = d - Ntrain_d - Nguard_d;
%         d_max = d + Ntrain_d + Nguard_d;
%         r_min = r - Ntrain_r - Nguard_r;
%         r_max = r + Ntrain_r + Nguard_r;
%         window = RD_beam_power(d_min:d_max,r_min:r_max);
%         guard_d_min = Ntrain_d + 1;
%         guard_d_max = Ntrain_d + 2*Nguard_d + 1;
%         guard_r_min = Ntrain_r + 1;
%         guard_r_max = Ntrain_r + 2*Nguard_r + 1;
%         window(guard_d_min:guard_d_max,guard_r_min:guard_r_max) = 0;
%         noise_level = sum(window(:))/Ntrain_total;
%         CFAR_threshold(d,r) = alpha_cfar*noise_level;
%         if RD_beam_power(d,r) > CFAR_threshold(d,r)
%             CFAR_map(d,r) = 1;
%         end
%     end
% end
% %--------------------------------------------------------------------------
% % CFAR DETECTION LIST
% %--------------------------------------------------------------------------
% [detections_d,detections_r] = find(CFAR_map == 1);
% Ndet = length(detections_d);
% disp('--------------------------------------------')
% disp(['Number of CA-CFAR detections: ' num2str(Ndet)])
% disp('--------------------------------------------')
% %--------------------------------------------------------------------------
% % CA-CFAR PLOT
% %--------------------------------------------------------------------------
% figure;
% imagesc(range_axis/1000,vel_axis,RD_beam_dB)
% axis xy
% xlabel('Range (km)')
% ylabel('Radial Velocity (m/s)')
% title('CA-CFAR Detections on Beamformed Range-Doppler Map')
% colorbar
% hold on
% plot(range_axis(detections_r)/1000,vel_axis(detections_d),'rx','MarkerSize',10,'LineWidth',2)
% plot(R0/1000,v0,'wo','MarkerSize',12,'LineWidth',2)
% legend('CA-CFAR Detections','True Target')
%
%--------------------------------------------------------------------------
% SENSOR-DOMAIN RANGE-DOPPLER CUBE WITHOUT MTD
%--------------------------------------------------------------------------
RD_cube = fftshift(fft(mf_out,Npulses,2),2);
%--------------------------------------------------------------------------
% TRUE MTD PROCESSOR: 3-PULSE CANCELLER BEFORE DOPPLER FFT
%--------------------------------------------------------------------------
mtd_out = zeros(size(mf_out));
for m = 1:M
    for r = 1:size(mf_out,3)
        for p = 3:Npulses
            mtd_out(m,p,r) = mf_out(m,p,r) - 2*mf_out(m,p-1,r) + mf_out(m,p-2,r);
        end
    end
end
%--------------------------------------------------------------------------
% SENSOR-DOMAIN RANGE-DOPPLER CUBE AFTER MTD
%--------------------------------------------------------------------------
RD_MTD_cube = fftshift(fft(mtd_out,Npulses,2),2);
%--------------------------------------------------------------------------
% TX-STEERED QUASI-COHERENT SENSOR INTEGRATION: WITHOUT MTD
%--------------------------------------------------------------------------
RD_beam = zeros(Npulses,size(RD_cube,3));
for m = 1:M
    RD_beam = RD_beam + squeeze(RD_cube(m,:,:))*conj(steer_beam(m));
end
RD_beam_power = abs(RD_beam).^2;
%--------------------------------------------------------------------------
% TX-STEERED QUASI-COHERENT SENSOR INTEGRATION: WITH MTD
%--------------------------------------------------------------------------
RD_MTD_beam = zeros(Npulses,size(RD_MTD_cube,3));
for m = 1:M
    RD_MTD_beam = RD_MTD_beam + squeeze(RD_MTD_cube(m,:,:))*conj(steer_beam(m));
end
RD_MTD_power = abs(RD_MTD_beam).^2;
%--------------------------------------------------------------------------
% COMMON NORMALIZATION FOR FAIR COMPARISON
%--------------------------------------------------------------------------
RD_beam_dB = 10*log10(RD_beam_power/max(RD_beam_power(:)) + eps);
RD_MTD_dB = 10*log10(RD_MTD_power/max(RD_beam_power(:)) + eps);
%--------------------------------------------------------------------------
% CA-CFAR PARAMETERS
%--------------------------------------------------------------------------
Ntrain_r = 20;
Ntrain_d = 4;
Nguard_r = 2;
Nguard_d = 1;
Pfa = 1e-6;
[Ndop,Nrange] = size(RD_beam_power);
Ntrain_total = (2*(Ntrain_d+Nguard_d)+1)*(2*(Ntrain_r+Nguard_r)+1) - (2*Nguard_d+1)*(2*Nguard_r+1);
alpha_cfar = Ntrain_total*(Pfa^(-1/Ntrain_total)-1);
%--------------------------------------------------------------------------
% BRANCH 1: CA-CFAR ONLY
%--------------------------------------------------------------------------
CFAR_map = zeros(Ndop,Nrange);
CFAR_threshold = zeros(Ndop,Nrange);
for d = Ntrain_d+Nguard_d+1:Ndop-(Ntrain_d+Nguard_d)
    for r = Ntrain_r+Nguard_r+1:Nrange-(Ntrain_r+Nguard_r)
        d_min = d - Ntrain_d - Nguard_d;
        d_max = d + Ntrain_d + Nguard_d;
        r_min = r - Ntrain_r - Nguard_r;
        r_max = r + Ntrain_r + Nguard_r;
        window = RD_beam_power(d_min:d_max,r_min:r_max);
        guard_d_min = Ntrain_d + 1;
        guard_d_max = Ntrain_d + 2*Nguard_d + 1;
        guard_r_min = Ntrain_r + 1;
        guard_r_max = Ntrain_r + 2*Nguard_r + 1;
        window(guard_d_min:guard_d_max,guard_r_min:guard_r_max) = 0;
        noise_level = sum(window(:))/Ntrain_total;
        CFAR_threshold(d,r) = alpha_cfar*noise_level;
        if RD_beam_power(d,r) > CFAR_threshold(d,r)
            CFAR_map(d,r) = 1;
        end
    end
end
[detections_d,detections_r] = find(CFAR_map);
%--------------------------------------------------------------------------
% BRANCH 2: MTD + CA-CFAR
%--------------------------------------------------------------------------
CFAR_MTD_map = zeros(Ndop,Nrange);
CFAR_MTD_threshold = zeros(Ndop,Nrange);
for d = Ntrain_d+Nguard_d+1:Ndop-(Ntrain_d+Nguard_d)
    for r = Ntrain_r+Nguard_r+1:Nrange-(Ntrain_r+Nguard_r)
        d_min = d - Ntrain_d - Nguard_d;
        d_max = d + Ntrain_d + Nguard_d;
        r_min = r - Ntrain_r - Nguard_r;
        r_max = r + Ntrain_r + Nguard_r;
        window = RD_MTD_power(d_min:d_max,r_min:r_max);
        guard_d_min = Ntrain_d + 1;
        guard_d_max = Ntrain_d + 2*Nguard_d + 1;
        guard_r_min = Ntrain_r + 1;
        guard_r_max = Ntrain_r + 2*Nguard_r + 1;
        window(guard_d_min:guard_d_max,guard_r_min:guard_r_max) = 0;
        noise_level = sum(window(:))/Ntrain_total;
        CFAR_MTD_threshold(d,r) = alpha_cfar*noise_level;
        if RD_MTD_power(d,r) > CFAR_MTD_threshold(d,r)
            CFAR_MTD_map(d,r) = 1;
        end
    end
end
[detections_MTD_d,detections_MTD_r] = find(CFAR_MTD_map);
%--------------------------------------------------------------------------
% BRANCH 3: VELOCITY-GATED MTD + CA-CFAR
%--------------------------------------------------------------------------
v_min_detect = 15;
v_max_detect = 20;
CFAR_GATED_map = zeros(Ndop,Nrange);
CFAR_GATED_threshold = zeros(Ndop,Nrange);
for d = Ntrain_d+Nguard_d+1:Ndop-(Ntrain_d+Nguard_d)
    if vel_axis(d) < v_min_detect || vel_axis(d) > v_max_detect
        continue
    end
    for r = Ntrain_r+Nguard_r+1:Nrange-(Ntrain_r+Nguard_r)
        d_min = d - Ntrain_d - Nguard_d;
        d_max = d + Ntrain_d + Nguard_d;
        r_min = r - Ntrain_r - Nguard_r;
        r_max = r + Ntrain_r + Nguard_r;
        window = RD_MTD_power(d_min:d_max,r_min:r_max);
        guard_d_min = Ntrain_d + 1;
        guard_d_max = Ntrain_d + 2*Nguard_d + 1;
        guard_r_min = Ntrain_r + 1;
        guard_r_max = Ntrain_r + 2*Nguard_r + 1;
        window(guard_d_min:guard_d_max,guard_r_min:guard_r_max) = 0;
        noise_level = sum(window(:))/Ntrain_total;
        CFAR_GATED_threshold(d,r) = alpha_cfar*noise_level;
        if RD_MTD_power(d,r) > CFAR_GATED_threshold(d,r)
            CFAR_GATED_map(d,r) = 1;
        end
    end
end
[detections_GATED_d,detections_GATED_r] = find(CFAR_GATED_map);

disp(['CFAR only detections:        ' num2str(length(detections_d))])
disp(['MTD + CFAR detections:       ' num2str(length(detections_MTD_d))])
disp(['Gated MTD + CFAR detections: ' num2str(length(detections_GATED_d))])

%--------------------------------------------------------------------------
% PLOT 1: CA-CFAR ONLY
%--------------------------------------------------------------------------
figure;
imagesc(range_axis/1000,vel_axis,RD_beam_dB)
axis xy
xlabel('Range (km)')
ylabel('Radial Velocity (m/s)')
title('CA-CFAR Only')
colorbar
hold on
plot(range_axis(detections_r)/1000,vel_axis(detections_d),'rx','MarkerSize',10,'LineWidth',2)
plot(R0/1000,v0,'wo','MarkerSize',12,'LineWidth',2)
legend('CFAR detections','True Target')
%--------------------------------------------------------------------------
% PLOT 2: MTD + CA-CFAR
%--------------------------------------------------------------------------
figure;
imagesc(range_axis/1000,vel_axis,RD_MTD_dB)
axis xy
xlabel('Range (km)')
ylabel('Radial Velocity (m/s)')
title('MTD + CA-CFAR')
colorbar
hold on
plot(range_axis(detections_MTD_r)/1000,vel_axis(detections_MTD_d),'rx','MarkerSize',10,'LineWidth',2)
plot(R0/1000,v0,'wo','MarkerSize',12,'LineWidth',2)
legend('MTD + CFAR detections','True Target')
%--------------------------------------------------------------------------
% PLOT 3: VELOCITY-GATED MTD + CA-CFAR
%--------------------------------------------------------------------------
figure;
imagesc(range_axis/1000,vel_axis,RD_MTD_dB)
axis xy
xlabel('Range (km)')
ylabel('Radial Velocity (m/s)')
title('Velocity-Gated MTD + CA-CFAR')
colorbar
hold on
plot(range_axis(detections_GATED_r)/1000,vel_axis(detections_GATED_d),'rx','MarkerSize',10,'LineWidth',2)
plot(R0/1000,v0,'wo','MarkerSize',12,'LineWidth',2)
yline(v_min_detect,'w--','LineWidth',1.5)
yline(v_max_detect,'w--','LineWidth',1.5)
legend('Gated detections','True Target','Velocity limits')

%--------------------------------------------------------------------------
% SELECT DETECTION LIST FOR DOA
%--------------------------------------------------------------------------
DOA_detections_d = detections_GATED_d;
DOA_detections_r = detections_GATED_r;
%--------------------------------------------------------------------------
% LOCAL DOA SEARCH REGION
%--------------------------------------------------------------------------
HPBW_az_deg = rad2deg(0.886*lambda/(Mx*dx));
HPBW_el_deg = rad2deg(0.886*lambda/(My*dy));
DOA_margin_az_deg = 2*HPBW_az_deg;
DOA_margin_el_deg = 2*HPBW_el_deg;
az_min_deg = 30;
az_max_deg = 150;
el_min_deg = 30;
el_max_deg = 150;
az_scan_deg = (az_beam-DOA_margin_az_deg):0.25:(az_beam+DOA_margin_az_deg);
el_scan_deg = (el_beam-DOA_margin_el_deg):0.25:(el_beam+DOA_margin_el_deg);
az_scan_deg = az_scan_deg(az_scan_deg >= az_min_deg & az_scan_deg <= az_max_deg);
el_scan_deg = el_scan_deg(el_scan_deg >= el_min_deg & el_scan_deg <= el_max_deg);
%--------------------------------------------------------------------------
% DETECTION CLUSTERING FOR DOA
%--------------------------------------------------------------------------
DOA_source = 'gated';
switch DOA_source
    case 'cfar'
        raw_detections_d = detections_d;
        raw_detections_r = detections_r;
        detection_power_map = RD_beam_power;
    case 'mtd'
        raw_detections_d = detections_MTD_d;
        raw_detections_r = detections_MTD_r;
        detection_power_map = RD_MTD_power;
    case 'gated'
        raw_detections_d = detections_GATED_d;
        raw_detections_r = detections_GATED_r;
        detection_power_map = RD_MTD_power;
end
cluster_radius_d = 1;
cluster_radius_r = 3;
remaining = true(length(raw_detections_d),1);
DOA_detections_d = [];
DOA_detections_r = [];
while any(remaining)
    idx_remaining = find(remaining);
    idx0 = idx_remaining(1);
    d0 = raw_detections_d(idx0);
    r0 = raw_detections_r(idx0);
    neighbor_idx = find(abs(raw_detections_d-d0) <= cluster_radius_d & abs(raw_detections_r-r0) <= cluster_radius_r & remaining);
    powers = zeros(length(neighbor_idx),1);
    for kk = 1:length(neighbor_idx)
        powers(kk) = detection_power_map(raw_detections_d(neighbor_idx(kk)),raw_detections_r(neighbor_idx(kk)));
    end
    [~,idx_peak_local] = max(powers);
    idx_peak_global = neighbor_idx(idx_peak_local);
    DOA_detections_d = [DOA_detections_d; raw_detections_d(idx_peak_global)];
    DOA_detections_r = [DOA_detections_r; raw_detections_r(idx_peak_global)];
    remaining(neighbor_idx) = false;
end
disp(['Raw detections before clustering: ' num2str(length(raw_detections_d))])
disp(['Detections after clustering:      ' num2str(length(DOA_detections_d))])
%--------------------------------------------------------------------------
% DOA ESTIMATION USING GATED MTD + CA-CFAR DETECTIONS
%--------------------------------------------------------------------------
for kdet = 1:length(DOA_detections_d)
    dop_bin = DOA_detections_d(kdet);
    range_bin = DOA_detections_r(kdet);
    snapshot = squeeze(RD_cube(:,dop_bin,range_bin));
    P_doa = zeros(length(el_scan_deg),length(az_scan_deg));
    for iaz = 1:length(az_scan_deg)
        for iel = 1:length(el_scan_deg)
            az_test = deg2rad(az_scan_deg(iaz));
            el_test = deg2rad(el_scan_deg(iel));
            ux_test = cos(el_test)*cos(az_test);
            uy_test = cos(el_test)*sin(az_test);
            a_test = exp(-1j*k0*(x_vec*ux_test + y_vec*uy_test));
            P_doa(iel,iaz) = abs(a_test'*snapshot)^2/(a_test'*a_test);
        end
    end
    P_doa_dB = 10*log10(P_doa/max(P_doa(:)) + eps);
    [~,idx_doa] = max(P_doa(:));
    [el_idx,az_idx] = ind2sub(size(P_doa),idx_doa);
    az_est = az_scan_deg(az_idx);
    el_est = el_scan_deg(el_idx);
    disp('--------------------------------------------')
    disp(['DOA detections processed:  ' num2str(length(DOA_detections_d))])
    disp(['Detection # ' num2str(kdet)])
    disp(['Range:               ' num2str(range_axis(range_bin)) ' m'])
    disp(['Velocity:            ' num2str(vel_axis(dop_bin)) ' m/s'])
    disp(['Estimated azimuth:   ' num2str(az_est) ' deg'])
    disp(['Estimated elevation: ' num2str(el_est) ' deg'])
    disp('--------------------------------------------')
    figure;
    imagesc(az_scan_deg,el_scan_deg,P_doa_dB)
    axis xy
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(['2D Bartlett DOA - Detection #' num2str(kdet)])
    colorbar
    hold on
    plot(az_est,el_est,'rx','MarkerSize',12,'LineWidth',2)
    plot(az_target,el_target,'wo','MarkerSize',10,'LineWidth',2)
    legend('Estimated DOA','True DOA')
    
    %--------------------------------------------------------------------------
    % DIGITAL MONOPULSE REFINEMENT USING NUMERICAL DERIVATIVE BEAMS
    %--------------------------------------------------------------------------
    mono_step_deg = 2;
    a0 = local_steer_vec(az_est,el_est,k0,x_vec,y_vec);
    a_az_p = local_steer_vec(az_est+mono_step_deg,el_est,k0,x_vec,y_vec);
    a_az_m = local_steer_vec(az_est-mono_step_deg,el_est,k0,x_vec,y_vec);
    a_el_p = local_steer_vec(az_est,el_est+mono_step_deg,k0,x_vec,y_vec);
    a_el_m = local_steer_vec(az_est,el_est-mono_step_deg,k0,x_vec,y_vec);
    Sigma = a0'*snapshot;
    Delta_az = (a_az_p-a_az_m)'*snapshot;
    Delta_el = (a_el_p-a_el_m)'*snapshot;
    mono_ratio_az = real(Delta_az/Sigma);
    mono_ratio_el = real(Delta_el/Sigma);
    %--------------------------------------------------------------------------
    % MONOPULSE SLOPE CALIBRATION
    %--------------------------------------------------------------------------
    a_true_az_p = a_az_p;
    a_true_az_m = a_az_m;
    a_true_el_p = a_el_p;
    a_true_el_m = a_el_m;
    Sigma_az_p = a0'*a_true_az_p;
    Sigma_az_m = a0'*a_true_az_m;
    Sigma_el_p = a0'*a_true_el_p;
    Sigma_el_m = a0'*a_true_el_m;
    Delta_az_p = (a_az_p-a_az_m)'*a_true_az_p;
    Delta_az_m = (a_az_p-a_az_m)'*a_true_az_m;
    Delta_el_p = (a_el_p-a_el_m)'*a_true_el_p;
    Delta_el_m = (a_el_p-a_el_m)'*a_true_el_m;
    ratio_az_p = real(Delta_az_p/Sigma_az_p);
    ratio_az_m = real(Delta_az_m/Sigma_az_m);
    ratio_el_p = real(Delta_el_p/Sigma_el_p);
    ratio_el_m = real(Delta_el_m/Sigma_el_m);
    K_az = (ratio_az_p-ratio_az_m)/(2*mono_step_deg);
    K_el = (ratio_el_p-ratio_el_m)/(2*mono_step_deg);
    az_mono_est = az_est + mono_ratio_az/K_az;
    el_mono_est = el_est + mono_ratio_el/K_el;
    disp(['Coarse azimuth:              ' num2str(az_est) ' deg'])
    disp(['Coarse elevation:            ' num2str(el_est) ' deg'])
    disp(['Monopulse azimuth ratio:     ' num2str(mono_ratio_az)])
    disp(['Monopulse elevation ratio:   ' num2str(mono_ratio_el)])
    disp(['Refined monopulse azimuth:   ' num2str(az_mono_est) ' deg'])
    disp(['Refined monopulse elevation: ' num2str(el_mono_est) ' deg'])

end
%
function ratio = local_mono_ratio_az(az_true,el_true,az_center,el_center,mono_sep_deg,k0,x_vec,y_vec)
az_L = az_center - mono_sep_deg;
az_R = az_center + mono_sep_deg;
el_D = el_center - mono_sep_deg;
el_U = el_center + mono_sep_deg;
a_true = local_steer_vec(az_true,el_true,k0,x_vec,y_vec);
a_UL = local_steer_vec(az_L,el_U,k0,x_vec,y_vec);
a_UR = local_steer_vec(az_R,el_U,k0,x_vec,y_vec);
a_LL = local_steer_vec(az_L,el_D,k0,x_vec,y_vec);
a_LR = local_steer_vec(az_R,el_D,k0,x_vec,y_vec);
B_UL = a_UL'*a_true;
B_UR = a_UR'*a_true;
B_LL = a_LL'*a_true;
B_LR = a_LR'*a_true;
Sigma = B_UL + B_UR + B_LL + B_LR;
Delta_az = (B_UR + B_LR) - (B_UL + B_LL);
ratio = real(Delta_az/Sigma);
end

function ratio = local_mono_ratio_el(az_true,el_true,az_center,el_center,mono_sep_deg,k0,x_vec,y_vec)
az_L = az_center - mono_sep_deg;
az_R = az_center + mono_sep_deg;
el_D = el_center - mono_sep_deg;
el_U = el_center + mono_sep_deg;
a_true = local_steer_vec(az_true,el_true,k0,x_vec,y_vec);
a_UL = local_steer_vec(az_L,el_U,k0,x_vec,y_vec);
a_UR = local_steer_vec(az_R,el_U,k0,x_vec,y_vec);
a_LL = local_steer_vec(az_L,el_D,k0,x_vec,y_vec);
a_LR = local_steer_vec(az_R,el_D,k0,x_vec,y_vec);
B_UL = a_UL'*a_true;
B_UR = a_UR'*a_true;
B_LL = a_LL'*a_true;
B_LR = a_LR'*a_true;
Sigma = B_UL + B_UR + B_LL + B_LR;
Delta_el = (B_UL + B_UR) - (B_LL + B_LR);
ratio = real(Delta_el/Sigma);
end

function a = local_steer_vec(az_deg,el_deg,k0,x_vec,y_vec)
az = deg2rad(az_deg);
el = deg2rad(el_deg);
ux = cos(el)*cos(az);
uy = cos(el)*sin(az);
a = exp(-1j*k0*(x_vec*ux + y_vec*uy));
end
%%
%--------------------------------------------------------------------------
% 2D PLOT - ONLY DOA POINTS IN TRUE ANGULAR COORDINATES
%--------------------------------------------------------------------------
az_plot_min = az_target - 10;
az_plot_max = az_target + 10;
el_plot_min = el_target - 10;
el_plot_max = el_target + 10;


az_grid_step = 0.25;
el_grid_step = 0.25;

az_axis = az_plot_min:az_grid_step:az_plot_max;
el_axis = el_plot_min:el_grid_step:el_plot_max;

figure
plot(az_target,el_target,'bo','MarkerSize',10,'LineWidth',2)
hold on
plot(az_est,el_est,'ko','MarkerSize',9,'LineWidth',2)
plot(az_mono_est,el_mono_est,'rx','MarkerSize',12,'LineWidth',2)
grid on
axis equal
xlim([az_plot_min az_plot_max])
ylim([el_plot_min el_plot_max])
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('True DOA, Coarse DOA and Refined Monopulse DOA')
legend('True DOA','Coarse DOA','Monopulse DOA','Location','best')
%%
%--------------------------------------------------------------------------
% 3D PLOT - FOUR NUMERICAL-DERIVATIVE MONOPULSE BEAMS OVER FULL FOV
%--------------------------------------------------------------------------
FOV_broadside = deg2rad(60);
u_vec = linspace(-sin(FOV_broadside),sin(FOV_broadside),201);
v_vec = linspace(-sin(FOV_broadside),sin(FOV_broadside),201);
[U,V] = meshgrid(u_vec,v_vec);
visible = U.^2 + V.^2 <= sin(FOV_broadside)^2;

AF_AZP_full = nan(size(U));
AF_AZM_full = nan(size(U));
AF_ELP_full = nan(size(U));
AF_ELM_full = nan(size(U));

for ia = 1:size(U,1)
    for ja = 1:size(U,2)
        if visible(ia,ja)
            steer_scan = exp(-1j*k0*(x_vec*U(ia,ja) + y_vec*V(ia,ja)));
            AF_AZP_full(ia,ja) = abs(a_az_p'*steer_scan)/M;
            AF_AZM_full(ia,ja) = abs(a_az_m'*steer_scan)/M;
            AF_ELP_full(ia,ja) = abs(a_el_p'*steer_scan)/M;
            AF_ELM_full(ia,ja) = abs(a_el_m'*steer_scan)/M;
        end
    end
end

AF_norm_full = max([ ...
    max(AF_AZP_full(:),[],'omitnan'), ...
    max(AF_AZM_full(:),[],'omitnan'), ...
    max(AF_ELP_full(:),[],'omitnan'), ...
    max(AF_ELM_full(:),[],'omitnan')]);

AF_AZP_full = AF_AZP_full/AF_norm_full;
AF_AZM_full = AF_AZM_full/AF_norm_full;
AF_ELP_full = AF_ELP_full/AF_norm_full;
AF_ELM_full = AF_ELM_full/AF_norm_full;

W = sqrt(max(0,1-U.^2-V.^2));

X_AZP = AF_AZP_full.*U;
Y_AZP = AF_AZP_full.*V;
Z_AZP = AF_AZP_full.*W;

X_AZM = AF_AZM_full.*U;
Y_AZM = AF_AZM_full.*V;
Z_AZM = AF_AZM_full.*W;

X_ELP = AF_ELP_full.*U;
Y_ELP = AF_ELP_full.*V;
Z_ELP = AF_ELP_full.*W;

X_ELM = AF_ELM_full.*U;
Y_ELM = AF_ELM_full.*V;
Z_ELM = AF_ELM_full.*W;

figure
surf(X_AZP,Y_AZP,Z_AZP,20*log10(AF_AZP_full+eps),'EdgeColor','none','FaceAlpha',0.35)
hold on
surf(X_AZM,Y_AZM,Z_AZM,20*log10(AF_AZM_full+eps),'EdgeColor','none','FaceAlpha',0.35)
surf(X_ELP,Y_ELP,Z_ELP,20*log10(AF_ELP_full+eps),'EdgeColor','none','FaceAlpha',0.35)
surf(X_ELM,Y_ELM,Z_ELM,20*log10(AF_ELM_full+eps),'EdgeColor','none','FaceAlpha',0.35)

quiver3(0,0,0,1.2,0,0,'k','LineWidth',2)
quiver3(0,0,0,0,1.2,0,'k','LineWidth',2)
quiver3(0,0,0,0,0,1.2,'k','LineWidth',2)

text(1.3,0,0,'x')
text(0,1.3,0,'y')
text(0,0,1.3,'z / broadside')

grid on
axis equal
xlabel('u_x')
ylabel('u_y')
zlabel('u_z / broadside')
title('Four Digital Monopulse Derivative Beampatterns over Full FoV')
colorbar
colormap jet
clim([-40 0])
view(45,25)
%legend('+Az beam','-Az beam','+El beam','-El beam')