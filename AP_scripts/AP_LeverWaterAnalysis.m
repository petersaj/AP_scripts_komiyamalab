%% Perform analysis for LeverWater

AP_getLeftLeverTimes;

%% plot percent press v trial

figure
subplot(2,2,1)
plot(left_down,[1:length(left_down)],'k.');
title('Lever press timing')
xlabel('Time (s)')
ylabel('Lever press')
subplot(2,2,2)
left_down_diff = diff(left_down);
plot([1:length(left_down_diff)],left_down_diff,'k.')
title('Time between lever presses')
xlabel('Lever press')
ylabel('Time between lever press (s)')

%% plot rewarded lever presses
subplot(2,2,3)
rewarded = find(states_up == 42);
plot(left_up(rewarded),[1:length(rewarded)],'k.');
title('Time of rewarded lever presses')
xlabel('Time (s)')
ylabel('Lever release #')

%% plot time between rewarded lever presses

subplot(2,2,4)
plot([1:length(rewarded)-1],diff(left_up(rewarded)),'k.');
title('Time between rewarded lever presses')
xlabel('Lever release #')
ylabel('Time from last reward (s)')