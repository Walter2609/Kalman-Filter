%% Kalman Filter für Wurfparabel
% Authors: Walter Schierl, Daniel Fetz, Philipp Wörz
clear all; close all; clc;
%% Parameter
g = 9.80665;
H = 0.005;
k = 0;
n = 134;
 
random = randi(n, [-10, 10]); 
%--------------------------------------------------------------------------
%% Für ideales System
t = [0];

x_ideal = [0];
y_ideal = [0];

v_0 = 3.8;
vx_0 = 0.5*v_0;
vy_0 = 0.86*v_0;

vx_ideal = [vx_0];
vy_ideal = [vy_0];
%--------------------------------------------------------------------------
%% Für Kalman Filter
mat_K = zeros(4,4);
% Startschätzung
mat_Pm = [ 1 0 0 0;
           0 1 0 0;
           0 0 1 0;
           0 0 0 1 ];

mat_R = [ 0.01  0    ;
          0     0.01 ];
mat_H = [ 1 0 0 0;
          0 1 0 0 ];      
      
% Startschätzung      
vec_Xm = [0 0 1 1].';
array_Xm = zeros(4, n+1);
array_Xm(1) = vec_Xm(1);
array_Xm(2) = vec_Xm(2);
array_Xm(3) = vec_Xm(3);
array_Xm(4) = vec_Xm(4);
x_kal = [0];
y_kal = [0];
x_mes = [0];
y_mes = [0];

vec_Xp = [0 0 0 0].';
vec_Y  = [0 0].';

mat_Pp = zeros(4,4);
mat_I  = [ 1 0 0 0;
           0 1 0 0;
           0 0 1 0;
           0 0 0 1 ];
       
mat_A = [ 1 0 H 0;
          0 1 0 H;
          0 0 1 0;
          0 0 0 1 ];
mat_B = [ 0     0;
          0   (0.5*H*H);
          0     0;
          0     H];

vec_u = [ 0 -g].';

mat_Q = 0.00001*mat_I;

%% For-Schleife
for k = 1 : n
    % ideales System
    t(k+1) = t(k) + H;
    
    x_ideal(k+1) = x_ideal(k) + H*vx_ideal(k);
    y_ideal(k+1) = y_ideal(k) + H*vy_ideal(k) - 0.5*g*H*H;
    
    vx_ideal(k+1) = vx_ideal(k);
    vy_ideal(k+1) = vy_ideal(k) - g*H;
    %____________________________________________________
    %Messung
    vec_Y(1) = x_ideal(k);
    noise     = (rand-0.5) /10;
    vec_Y(2) = y_ideal(k) + noise;
    
    %Wurfparabel mit Störung
    if( 0.7 < x_ideal(k) && x_ideal(k) < 0.8 )
        vec_Y(2) = 0.1 + noise;
        mat_R(1) = 1;
        mat_R(4) = 1;
    else
        mat_R(1) = 0.01;
        mat_R(4) = 0.01;
    end
    
    if( vec_Y(2) < 0)
        vex_Y(2) = 0;
    end
    x_mes(k) = vec_Y(1);
    y_mes(k) = vec_Y(2);
       
    %_______________________________________________________
    % Korrektur mit der Messung
    Inverse = inv(mat_R + mat_H*mat_Pm*mat_H');
    mat_K   = mat_Pm*mat_H' * Inverse;
    vec_Xp  = vec_Xm + mat_K * (vec_Y -mat_H*vec_Xm);
    
    mat_Pp  = (mat_I -mat_K*mat_H)*mat_Pm;
    
    %_______________________________________________________
    % Prädiktion
    vec_Xm = mat_A*vec_Xp + mat_B*vec_u;
    
    array_Xm((k*4)+1) = vec_Xm(1);
    array_Xm((k*4)+2) = vec_Xm(2);
    array_Xm((k*4)+3) = vec_Xm(3);
    array_Xm((k*4)+4) = vec_Xm(4);
    x_kal(k) = vec_Xm(1);
    y_kal(k) = vec_Xm(2);
    
    mat_Pm = mat_A*mat_Pp*mat_A' + mat_Q;  
end

% ideal 
subplot(2,1,1);
plot(x_ideal, y_ideal, 'LineWidth', 2);
axis([0, 1.4, 0, 0.8]);
title('Ideale Wurfparabel' , 'FontSize', 14);
xlabel('x/m', 'FontSize', 14);
ylabel('y/m', 'FontSize', 14);

str1 = ['v_{x0} = ' num2str(vx_0) 'm/s'];
str2 = ['v_{y0} = ' num2str(vy_0) 'm/s'];
str3 = ['t_{end} = ' num2str(t(k)) 's'];
text(0.03, 0.7, str1, 'FontSize', 12);
text(0.03, 0.6, str2, 'FontSize', 12);
text(0.03, 0.5, str3, 'FontSize', 12);

% Kalman Filter
subplot(2,1,2);
plot(x_ideal, y_ideal, '--', x_mes, y_mes, 'g+', x_kal, y_kal, 'r', 'LineWidth', 2);
axis([0, 1.4, 0, 0.8]);

xlabel('x/m', 'FontSize', 14);
ylabel('y/m', 'FontSize', 14);
title('Rekonstruierte Wurfparabel durch Kalman-Filter mit Suchfenster' , 'FontSize', 14);    

str1 = ['v_{x0} = ' num2str(vx_0) 'm/s'];
str2 = ['v_{y0} = ' num2str(vy_0) 'm/s'];
str3 = ['t_{end} = ' num2str(t(k)) 's'];
text(0.03, 0.7, str1, 'FontSize', 12);
text(0.03, 0.6, str2, 'FontSize', 12);
text(0.03, 0.5, str3, 'FontSize', 12);
    
    
    
    
    
    
    
    
    
    
    
    
      
      
      
       
       
       
       
       
       
      
      
      
      
      
       