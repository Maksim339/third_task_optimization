q1 = 0; % Угол звена 1
q2 = 0; % Угол звена 2
q3 = 0; % Угол звена 3
[d, R] = manipulator_model(q1, q2, q3);

% Вывод результатов
% disp('Положение схвата (d):');
% disp(d);
% disp('Ориентация схвата (R):');
% disp(R);
% disp([R d; 0 0 0 1])


load('calib.mat');

Q = calib.Q;
T_actual = calib.T;

for i = 1:size(Q, 1)
    q1 = Q(i, 1); % Угол звена 1
    q2 = Q(i, 2); % Угол звена 2
    q3 = Q(i, 3); % Угол звена 3

    [d, R] = manipulator_model(q1, q2, q3);
    T_calculated = [R d; 0 0 0 1];

    disp(['Матрица преобразования для набора углов ', num2str(i), ':']);
    disp(T_calculated);

    disp(['Фактическая матрица преобразования для набора углов ', num2str(i), ':']);
    disp(T_actual(:, :, i));

    % Сравнение с фактической матрицей из T
    diff = norm(T_calculated - T_actual(:, :, i), 'fro');
    disp(['Разница с фактической матрицей: ', num2str(diff)]);
end



function Rz = rotz(theta)
    Rz = [cos(theta) -sin(theta) 0 0;
          sin(theta) cos(theta) 0 0;
          0 0 1 0;
          0 0 0 1];
end

function Tz = transz(d)
    Tz = [1 0 0 0;
          0 1 0 0;
          0 0 1 d;
          0 0 0 1];
end

function Rx = rotx(alpha)
    Rx = [1 0 0 0;
          0 cos(alpha) -sin(alpha) 0;
          0 sin(alpha) cos(alpha) 0;
          0 0 0 1];
end

function Tx = transx(a)
    Tx = [1 0 0 a;
          0 1 0 0;
          0 0 1 0;
          0 0 0 1];
end

function [d, R] = manipulator_model(q1, q2, q3)
    % Звено 1
    d1 = 300;
    a1 = 0;
    alpha1 = pi/2;

    % Звено 2
    d2 = 0;
    a2 = 400;
    alpha2 = 0;

    % Звено 3
    d3 = 0;
    a3 = 500;
    alpha3 = 0;

    T1 = dh_transform(q1, d1, a1, alpha1);
    T2 = dh_transform(q2, d2, a2, alpha2);
    T3 = dh_transform(q3, d3, a3, alpha3);

    % Общая матрица преобразования
    T = T1 * T2 * T3;

    d = T(1:3, 4);
    R = T(1:3, 1:3);
end


function T = dh_transform(theta, d, a, alpha)
    % T = [cos(theta), sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
    %      sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    %      0, sin(alpha), cos(alpha), d;
    %      0, 0, 0, 1];
    T = rotz(theta)*transz(d)*rotx(alpha)*transx(a);
          
end
