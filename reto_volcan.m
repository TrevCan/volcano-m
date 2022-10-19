clear all;

%Constantes
m = 15 ;
cd = 0.47;
A = 0.8;
dt = 0.001;
step = dt:dt:25;
g = 9.81;


%Functions%
p = @(y) y*-7.8*10^(-5)+1.225;
k = @(y) cd*A*p(y);
angle = @(vx, vy) atan2(vy,vx);
getay = @(y, vx, vy, v) ( -k(y)*v^2*sin(angle(vx,vy)) / m - g );
getax = @(y, vx, vy, v) (-k(y)*v^2*cos(angle(vx,vy))) / m;
getvy = @(pvy, dt, ay) pvy+ay*dt;
getvx = @(pvx, dt, ax) pvx+ax*dt;
getv = @(vx, vy) sqrt(vx^2 + vy^2);
getx = @(px, dt, vx) px+vx*dt;
gety = @(py, dt, vy) py+vy*dt;

%Constantes no tan constantes (input)
% v = zeros(1, length(step))
% ---   OBJ 1 
v(1,1) = 35;
y(1,1) = 10;

% ---   OBJ 2
v(2,1) = rand(1)*20 + 10;
y(2,1) = 10;

% ---   OBJ 3 
v(3,1) = rand(1)*20 + 10;
y(3,1) = 10;


%Variables

% ANGLES
%
% ang = zeros(1, length(step))
% ---   OBJ 1 
ang(1,1) = deg2rad(55);

% ---   OBJ 2 
ang(2,1) = deg2rad( 20 + rand(1)* 50 );

% ---   OBJ 3 
ang(3,1) = deg2rad( 20 + rand(1)* 50 );

%       VELOCIDADES
 
% vx = zeros(1, length(step))
% ---   OBJ 1
vx(1,1) = v(1,1)*cos(ang(1,1));
% vy = zeros(1, length(step))
vy(1,1) = v(1,1)*sin(ang(1,1));
x(1,1) = 0;

% ---   OBJ 2
vx(2,1) = v(2,1)*cos(ang(2,1));
% vy = zeros(1, length(step))
vy(2,1) = v(2,1)*sin(ang(2,1));
x(2,1) = 0;

% ---   OBJ 3
vx(3,1) = v(3,1)*cos(ang(3,1));
% vy = zeros(1, length(step))
vy(3,1) = v(3,1)*sin(ang(3,1));
x(3,1) = 0;

%       ACCELERATION

% ax = zeros(1, length(step))
% ay = zeros(1, length(step))
% ---   OBJ 1
ax(1,1) = getax(y(1,1), vx(1,1), vy(1,1), v(1,1));
ay(1,1) = getay(y(1,1), vx(1,1), vy(1,1), v(1,1));

% ---   OBJ 2
ax(2,1) = getax(y(2,1), vx(2,1), vy(2,1), v(2,1));
ay(2,1) = getay(y(2,1), vx(2,1), vy(2,1), v(2,1));

% ---   OBJ 3
ax(3,1) = getax(y(3,1), vx(3,1), vy(3,1), v(3,1));
ay(3,1) = getay(y(3,1), vx(3,1), vy(3,1), v(3,1));



%Volcan plot



% ax(2) = getax(y(2-1), vx(2-1), vy(2-1), v(2-1))




for i=2:1:length(step)
    % ---   OBJ 1
    ax(1, i) = getax(y(1, i-1), vx(1, i-1), vy(1, i-1), v(1, i-1));
    ay(1, i) = getay(y(1, i-1), vx(1, i-1), vy(1, i-1), v(1, i-1));
    vx(1, i) = getvx(vx(1, i-1), dt, ax(1, i));
    vy(1, i) = getvy(vy(1, i-1), dt, ay(1, i));
    ang(1, i) = angle(vx(1, i), vy(1, i));
    x(1, i) = getx(x(1, i-1), dt, vx(1, i));
    y(1, i) = gety(y(1, i-1), dt, vy(1, i));
    v(1, i) = getv(vx(1, i), vy(1, i));

    % ---   OBJ 2
    ax(2,i) = getax(y(2, i-1), vx(2, i-1), vy(2, i-1), v(2, i-1));
    ay(2, i) = getay(y(2, i-1), vx(2, i-1), vy(2, i-1), v(2, i-1));
    vx(2, i) = getvx(vx(2, i-1), dt, ax(2, i));
    vy(2, i) = getvy(vy(2, i-1), dt, ay(2, i));
    ang(2, i) = angle(vx(2, i), vy(2, i));
    x(2, i) = getx(x(2, i-1), dt, vx(2, i));
    y(2, i) = gety(y(2, i-1), dt, vy(2, i));
    v(2, i) = getv(vx(2, i), vy(2, i));

    % ---   OBJ 3
    ax(3,i) = getax(y(3, i-1), vx(3, i-1), vy(3, i-1), v(3, i-1));
    ay(3, i) = getay(y(3, i-1), vx(3, i-1), vy(3, i-1), v(3, i-1));
    vx(3, i) = getvx(vx(3, i-1), dt, ax(3, i));
    vy(3, i) = getvy(vy(3, i-1), dt, ay(3, i));
    ang(3, i) = angle(vx(3, i), vy(3, i));
    x(3, i) = getx(x(3, i-1), dt, vx(3, i));
    y(3, i) = gety(y(3, i-1), dt, vy(3, i));
    v(3, i) = getv(vx(3, i), vy(3, i));

    if y(1,i) < 0 || y(2,i) < 0 || y(3,i) < 0 
        break
    end
    
%     t=-3:0.01:0;
%     plot(t-3, exp(t)*y(1))
%     hold
%     plot(-t+3,exp(t)*y(1))
%     plot(-3:0.5:3,(0.1*sin(5*(-3:0.5:3))+y(1)))
%     plot(x(i), y(i), 'O') 
%     axis([-5 100 0 40])
%     drawnow limitrate
%     hold off

end

% plot(x, y);

% 100 puntos dentro de length(step)
% 

frameCount = round(length(x(1,:))/500):length(x(1,:));

for i=1:round(length(x(1,:))/500):length(x(1,:))
    
    t=-3:0.01:0;
    plot(t-3, exp(t)*y(1,1));
    hold;
    plot(-t+3,exp(t)*y(1,1));
    plot(-3:0.5:3,(0.1*sin(5*(-3:0.5:3))+y(1,1)));

    plot(x(1,i), y(1,i), 'Ored') ;
    plot(x(2,i), y(2,i), 'Ogreen') ;
    plot(x(3,i), y(3,i), 'Oblue') ;

    axis([-5 100 0 40]);
    drawnow;
    hold off;

    frames(i) = getframe(gcf);

end


% VELOCIDADES rocas pequenas de objetos

length_v = length(v(1,:));
% ---   OBJ 1
o1v(1,1) = rand(1)*15 + 5;
o1y(1,1) = y(1,length_v);
o1x(1,1) = x(1,length_v);
o1angle(1,1) = deg2rad( rand(1)*140 + 20 ) ;
o1vx(1,1) = o1v(1,1)*cos(o1angle(1,1));
o1vy(1,1) = o1v(1,1)*sin(o1angle(1,1));
o1ax(1,1) = getax(o1y(1,1), o1vx(1,1), o1vy(1,1), o1v(1,1));
o1ay(1,1) = getay(o1y(1,1), o1vx(1,1), o1vy(1,1), o1v(1,1));

o1v(2,1) = rand(1)*15 + 5;
o1y(2,1) = y(1,length_v);
o1x(2,1) = x(1,length_v);
o1angle(2,1) = deg2rad( rand(1)*140 + 20 ) ;
o1vx(2,1) = o1v(2,1)*cos(o1angle(2,1));
o1vy(2,1) = o1v(2,1)*sin(o1angle(2,1));
o1ax(2,1) = getax(o1y(2,1), o1vx(2,1), o1vy(2,1), o1v(2,1));
o1ay(2,1) = getay(o1y(2,1), o1vx(2,1), o1vy(2,1), o1v(2,1));


o1v(3,1) = rand(1)*15 + 5;
o1y(3,1) = y(1,length_v);
o1x(3,1) = x(1,length_v);
o1angle(3,1) = deg2rad( rand(1)*140 + 20 ) ;
o1vx(3,1) = o1v(3,1)*cos(o1angle(3,1));
o1vy(3,1) = o1v(3,1)*sin(o1angle(3,1));
o1ax(3,1) = getax(o1y(3,1), o1vx(3,1), o1vy(3,1), o1v(3,1));
o1ay(3,1) = getay(o1y(3,1), o1vx(3,1), o1vy(3,1), o1v(3,1));


% ---   OBJ 2
o2v(1,1) = rand(1)*15 + 5;
o2y(1,1) = y(2,length_v);
o2x(1,1) = x(2,length_v);
o2angle(1,1) = deg2rad( rand(1)*140 + 20 ) ;
o2vx(1,1) = o2v(1,1)*cos(o2angle(1,1));
o2vy(1,1) = o2v(1,1)*sin(o2angle(1,1));
o2ax(1,1) = getax(o2y(1,1), o2vx(1,1), o2vy(1,1), o2v(1,1));
o2ay(1,1) = getay(o2y(1,1), o2vx(1,1), o2vy(1,1), o2v(1,1));

o2v(2,1) = rand(1)*15 + 5;
o2y(2,1) = y(2,length_v);
o2x(2,1) = x(2,length_v);
o2angle(2,1) = deg2rad( rand(1)*140 + 20 ) ;
o2vx(2,1) = o2v(2,1)*cos(o2angle(2,1));
o2vy(2,1) = o2v(2,1)*sin(o2angle(2,1));
o2ax(2,1) = getax(o2y(2,1), o2vx(2,1), o2vy(2,1), o2v(2,1));
o2ay(2,1) = getay(o2y(2,1), o2vx(2,1), o2vy(2,1), o2v(2,1));

o2v(3,1) = rand(1)*15 + 5;
o2y(3,1) = y(2,length_v);
o2x(3,1) = x(2,length_v);
o2angle(3,1) = deg2rad( rand(1)*140 + 20 ) ;
o2vx(3,1) = o2v(3,1)*cos(o2angle(3,1));
o2vy(3,1) = o2v(3,1)*sin(o2angle(3,1));
o2ax(3,1) = getax(o2y(3,1), o2vx(3,1), o2vy(3,1), o2v(3,1));
o2ay(3,1) = getay(o2y(3,1), o2vx(3,1), o2vy(3,1), o2v(3,1));

% ---   OBJ 3 
o3v(1,1) = rand(1)*15 + 5;
o3y(1,1) = y(3,length_v);
o3x(1,1) = x(3,length_v);
o3angle(1,1) = deg2rad( rand(1)*140 + 20 ) ;
o3vx(1,1) = o3v(1,1)*cos(o3angle(1,1));
o3vy(1,1) = o3v(1,1)*sin(o3angle(1,1));
o3ax(1,1) = getax(o3y(1,1), o3vx(1,1), o3vy(1,1), o3v(1,1));
o3ay(1,1) = getay(o3y(1,1), o3vx(1,1), o3vy(1,1), o3v(1,1));

o3v(2,1) = rand(1)*15 + 5;
o3y(2,1) = y(3,length_v);
o3x(2,1) = x(3,length_v);
o3angle(2,1) = deg2rad( rand(1)*140 + 20 ) ;
o3vx(2,1) = o3v(2,1)*cos(o3angle(2,1));
o3vy(2,1) = o3v(2,1)*sin(o3angle(2,1));
o3ax(2,1) = getax(o3y(2,1), o3vx(2,1), o3vy(2,1), o3v(2,1));
o3ay(2,1) = getay(o3y(2,1), o3vx(2,1), o3vy(2,1), o3v(2,1));

o3v(3,1) = rand(1)*15 + 5;
o3y(3,1) = y(3,length_v);
o3x(3,1) = x(3,length_v);
o3angle(3,1) = deg2rad( rand(1)*140 + 20 ) ;
o3vx(3,1) = o3v(3,1)*cos(o3angle(3,1));
o3vy(3,1) = o3v(3,1)*sin(o3angle(3,1));
o3ax(3,1) = getax(o3y(3,1), o3vx(3,1), o3vy(3,1), o3v(3,1));
o3ay(3,1) = getay(o3y(3,1), o3vx(3,1), o3vy(3,1), o3v(3,1));


% calculate values o1..., o2..., o3...

for i=2:1:length(step)
    % ---   OBJ 1
    for index=1:1:3
        o1ax(index, i) = getax(o1y(index, i-1), o1vx(index, i-1), o1vy(index, i-1), o1v(index, i-1));
        o1ay(index, i) = getay(o1y(index, i-1), o1vx(index, i-1), o1vy(index, i-1), o1v(index, i-1));
        o1vx(index, i) = getvx(o1vx(index, i-1), dt, o1ax(index, i));
        o1vy(index, i) = getvy(o1vy(index, i-1), dt, o1ay(index, i));
        o1angle(index, i) = angle(o1vx(index, i), o1vy(index, i));
        o1x(index, i) = getx(o1x(index, i-1), dt, o1vx(index, i));
        o1y(index, i) = gety(o1y(index, i-1), dt, o1vy(index, i));
        o1v(index, i) = getv(o1vx(index, i), o1vy(index, i));
    end

    % ---   OBJ 2
    for index=1:1:3
        o2ax(index, i) = getax(o2y(index, i-1), o2vx(index, i-1), o2vy(index, i-1), o2v(index, i-1));
        o2ay(index, i) = getay(o2y(index, i-1), o2vx(index, i-1), o2vy(index, i-1), o2v(index, i-1));
        o2vx(index, i) = getvx(o2vx(index, i-1), dt, o2ax(index, i));
        o2vy(index, i) = getvy(o2vy(index, i-1), dt, o2ay(index, i));
        o2angle(index, i) = angle(o2vx(index, i), o2vy(index, i));
        o2x(index, i) = getx(o2x(index, i-1), dt, o2vx(index, i));
        o2y(index, i) = gety(o2y(index, i-1), dt, o2vy(index, i));
        o2v(index, i) = getv(o2vx(index, i), o2vy(index, i));
    end

    
    % ---   OBJ 3
    for index=1:1:3
        o3ax(index, i) = getax(o3y(index, i-1), o3vx(index, i-1), o3vy(index, i-1), o3v(index, i-1));
        o3ay(index, i) = getay(o3y(index, i-1), o3vx(index, i-1), o3vy(index, i-1), o3v(index, i-1));
        o3vx(index, i) = getvx(o3vx(index, i-1), dt, o3ax(index, i));
        o3vy(index, i) = getvy(o3vy(index, i-1), dt, o3ay(index, i));
        o3angle(index, i) = angle(o3vx(index, i), o3vy(index, i));
        o3x(index, i) = getx(o3x(index, i-1), dt, o3vx(index, i));
        o3y(index, i) = gety(o3y(index, i-1), dt, o3vy(index, i));
        o3v(index, i) = getv(o3vx(index, i), o3vy(index, i));
    end


    if o1y(1,i) < -1 || o1y(2,i) < -1 || o1y(3,i) < -1 || o2y(1,i) < -1 || o2y(2,i) < -1 || o2y(3,i) < -1 || o3y(1,i) < -1 || o3y(2,i) < -1 || o3y(3,i) < -1 
%         disp("BREAKINGGG... ! ")
%         disp(o1y(1,i))
%         disp(o1y(2,i))
%         disp(o1y(3,i))
% 
%         disp(o2y(1,i))
%         disp(o2y(2,i))
%         disp(o2y(3,i))
% 
%         disp(o3y(1,i))
%         disp(o3y(2,i))
%         disp(o3y(3,i))

        break
    end

end

for i=1:round(length(o1x(1,:))/500):length(o1x(1,:))
    
    t=-3:0.01:0;
    plot(t-3, exp(t)*y(1,1));
    hold;
    plot(-t+3,exp(t)*y(1,1));
    plot(-3:0.5:3,(0.1*sin(5*(-3:0.5:3))+y(1,1)));

    plot(o1x(1,i), o1y(1,i), '.red') ;
    plot(o1x(2,i), o1y(2,i), '.red') ;
    plot(o1x(3,i), o1y(3,i), '.red') ;

    plot(o2x(1,i), o2y(1,i), '.green') ;
    plot(o2x(2,i), o2y(2,i), '.green') ;
    plot(o2x(3,i), o2y(3,i), '.green') ;

    plot(o3x(1,i), o3y(1,i), '.blue') ;
    plot(o3x(2,i), o3y(2,i), '.blue') ;
    plot(o3x(3,i), o3y(3,i), '.blue') ;

    axis([-5 100 0 40]);
    drawnow;
    hold off;

    frames(i+frameCount) = getframe(gcf);

end

video = VideoWriter("reto_volcan_cumbias.avi");
video.FrameRate = 24;
open(video);
writeVideo(video, frames);
close(video);
