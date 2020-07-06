% Here we will make a video based on U
close all;

fprintf('\n Assembling avi file \n');

title = sprintf('LTHeatedplateM%d-T%dto%.2fP%d.avi',mesh_choice,0,tFin,p);

framerate = 60; % frames per second
frame_step = floor( size(U,2)/(tFin*framerate));
zlims = [0,2];  % desired z limits
cameraPos = [-0.9507,  -16.7937,  5.1315]; % desired camera position
cameraAngle = 8; % desired camera angle
clims = [0,2];
pos = [50, 50, 560, 420];
View_Soln(U(:,1), Elements, Nodes, map, sgn, p, Psi,how_rough);

 
%% Set up the movie.
writerObj = VideoWriter(title); % Name it.
writerObj.FrameRate = 24; % How many frames per second.
open(writerObj); 
 
for i=1:frame_step:size(U,2) 
    % We just use pause but pretend you have some really complicated thing here...
View_Soln(U(:,i), Elements, Nodes, map, sgn, p, Psi, how_rough);
a1 = gca; % copy current axes
fid = figure('Visible','Off'); % create new figure
a2 = copyobj(a1,fid);
close(1);
fid.CurrentAxes.ZLim = zlims;
fid.CurrentAxes.CameraPosition = cameraPos;
fid.CurrentAxes.CameraViewAngle = cameraAngle;
fid.CurrentAxes.CLim = clims;
fid.Position = pos;
% figure(fid);
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 close all;
 fprintf('On frame %d of %d \n', i, size(U,2));
    pause(0.0001);
end
close(writerObj); % Saves the movie.