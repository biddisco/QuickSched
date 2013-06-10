
%% Some global stuff
tpms = 1/2.6e6;
ms_max = 30;

%% Plot the task timelines for tasks allocation
% Load the data
tasks = importdata( 'test.dump' );
tasks(:,6) = ( tasks(:,6) - tasks(:,5) ) * tpms;
start = min( tasks(:,5) );
tasks(:,5) = ( tasks(:,5) - start ) * tpms;
nr_cores = max( tasks(:,2) ) + 1;

% Init the plot
clf;
subplot('position',[ 0.05 , 0.1 , 0.9 , 0.8 ]);
hold on;

% Plot the tasks
for k=1:size(tasks,1)
    rectangle( 'Position' , [ tasks(k,5) , tasks(k,2)+0.5 , tasks(k,6) , 1 ] , ...
        'EdgeColor' , [ 0 0.8 0 ] , 'LineWidth' , 1 , 'FaceColor' , [ tasks(k,3)/3 , tasks(k,4)/3 , 0 ] );
    text( tasks(k,5) + tasks(k,6)*0.5 , tasks(k,2)+1 , ...
        sprintf( '%i,%i' , tasks(k,3) , tasks(k,4) ) , ...
        'HorizontalAlignment' , 'Center' );
end

% Set the axes and stuff.
hold off;
xlabel('time (ms)');
ylabel('core ID');
set(gca,'YTick',1:(max(tasks(:,1))+1))
title('QuickSched tasks');
axis([ 0 , max( tasks(:,5) + tasks(:,6) ) , 0.5 , nr_cores+0.5 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 16 4 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 16 4 ] );
print -depsc2 tasks_dynamic.eps
!epstopdf tasks_dynamic.eps 



