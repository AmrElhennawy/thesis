%%%% Plot z coordinate
figure(1)
subplot(3,2,5)
plot(DesiredTrack.Time, DesiredTrack.Data(:,7), 'r' , SystemStates.Time, SystemStates.Data(:, 5), 'b')
grid;
% title('z Controller');
xlabel('Time (sec)'); ylabel('z (m)')
% subplot(2,2,2)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,8), 'r' , SystemStates.Time, SystemStates.Data(:, 6), 'b')
% grid; legend('Desired Track','Actual')
% title('Altitude Controller'); xlabel('Time (sec)'); ylabel('Vertical Velocity in meters/sec')

%% Plot X Controller
% figure(3)
subplot(3,2,1)
plot(DesiredTrack.Time, DesiredTrack.Data(:,1), 'r' , SystemStates.Time, SystemStates.Data(:, 1), 'b')
grid; 
% legend('Desired Track','Actual')
% title('x Controller');
xlabel('Time (sec)'); ylabel('x (m)')
% subplot(2,1,2)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,2), 'r' , SystemStates.Time, SystemStates.Data(:, 2), 'b')
% grid; legend('Desired Track','Actual')
% title('X Controller'); xlabel('Time (sec)'); ylabel('x speed in m/sec')
% figure(4)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,3)), grid

%% Plot Y Controller
subplot(3,2,3)
plot(DesiredTrack.Time, DesiredTrack.Data(:,4), 'r' , SystemStates.Time, SystemStates.Data(:, 3), 'b')
grid; 
% title('y Controller');
xlabel('Time (sec)'); ylabel('y (m)')
% subplot(2,2,4)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,5), 'r' , SystemStates.Time, SystemStates.Data(:, 4), 'b')
% grid; legend('Desired Track','Actual')
% title('Y Controller'); xlabel('Time (sec)'); ylabel('Y in m/sec')

%% Plot Psi Controller
subplot(3,2,6)
psi = round(SystemStates.Data(:, 11),1)
plot(DesiredTrack.Time, DesiredTrack.Data(:,16), 'r' , SystemStates.Time,psi, 'b')
grid; 
% legend('Desired Track','Actual')
% title('Psi Controller');
xlabel('Time (sec)'); ylabel('Psi (rad)')
% subplot(2,2,4)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,17), 'r' , SystemStates.Time, SystemStates.Data(:, 12), 'b')
% grid; legend('Desired Track','Actual')
% title('Psi Controller'); xlabel('Time (sec)'); ylabel('Psi in radians/sec')

%% Plot theta Controller
% figure(2)
subplot(3,2,4)
plot(DesiredTrack.Time, DesiredTrack.Data(:,13), 'r' , SystemStates.Time, SystemStates.Data(:, 9), 'b')
grid; 
% legend('Desired Track','Actual')
% title('Theta Controller');
xlabel('Time (sec)'); ylabel('theta (rad)')
% subplot(2,2,4)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,14), 'r' , SystemStates.Time, SystemStates.Data(:, 10), 'b')
% grid; legend('Desired Track','Actual')
% title('Theta Controller'); xlabel('Time (sec)'); ylabel('theta dot in rad/sec')
%% Plot Phi Controller
subplot(3,2,2)
plot(DesiredTrack.Time, DesiredTrack.Data(:,10), 'r' , SystemStates.Time, SystemStates.Data(:, 7), 'b')
grid; legend('Desired Track','Actual')
% title('Phi Controller');
xlabel('Time (sec)'); ylabel('Phi (rad)')
% subplot(2,2,2)
% plot(DesiredTrack.Time, DesiredTrack.Data(:,11), 'r' , SystemStates.Time, SystemStates.Data(:, 8), 'b')
% grid; legend('Desired Track','Actual')
% title('Phi Controller'); xlabel('Time (sec)'); ylabel('Phi in radians/sec')

%% 3D plot
figure(3)
title('3D Tracking')
% plot3(SystemStates.Data(:, 1), SystemStates.Data(:, 3), SystemStates.Data(:, 5), 'LineStyle', '--', DesiredTrack.Data(:,1), DesiredTrack.Data(:,4), DesiredTrack.Data(:,7), 'LineStyle', '-.')
plot3(SystemStates.Data(:, 1), SystemStates.Data(:, 3), SystemStates.Data(:, 5), DesiredTrack.Data(:,1), DesiredTrack.Data(:,4), DesiredTrack.Data(:,7))
xlabel('x')
ylabel('y')
zlabel('z')
legend('Actual Track','Desired Track')
grid

%% Control Input
figure(4)
title('Control Action')
subplot(4,1,1)
plot(Us.Time, Us.Data(:,1)), grid
ylabel('U1')
% ylim([8 15])
subplot(4,1,2)
plot(Us.Time, Us.Data(:,2)), grid, ylim([-0.02 0.02])
ylabel('U2')
ylim([-2 2])
subplot(4,1,3)
plot(Us.Time, Us.Data(:,3)), grid, ylim([-0.02 0.02])
ylabel('U3')
ylim([-2 2])
subplot(4,1,4)
plot(Us.Time, Us.Data(:,4)), grid
ylabel('U4')
xlabel('Time (Sec)')
ylim([-0.5 0.5])