
% Пример использования функции
q1 = 0; % Угол звена 1
q2 = 0; % Угол звена 2
q3 = 0; % Угол звена 3
[d, R] = manipulator_model(q1, q2, q3);

% Вывод результатов
disp('Положение схвата (d):');
disp(d);
disp('Ориентация схвата (R):');
disp(R);



% Функции для матриц преобразований
function Rz = rotz(theta)
    Rz = [cos(theta) -sin(theta) 0 0;
          -sin(theta) cos(theta) 0 0;
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

% Функция для построения модели манипулятора и вычисления прямой кинематики
function [d, R] = manipulator_model(q1, q2, q3)
    % DH параметры манипулятора
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

    % Вычисление матриц преобразований для каждого звена
    T1 = dh_transform(q1, d1, a1, alpha1);
    T2 = dh_transform(q2, d2, a2, alpha2);
    T3 = dh_transform(q3, d3, a3, alpha3);

    % Общая матрица преобразования
    T = T1 * T2 * T3;

    % Выходные данные: положение и ориентация схвата манипулятора
    d = T(1:3, 4);
    R = T(1:3, 1:3);
end

% Функция для создания матрицы преобразования Денавита-Хартенберга
function T = dh_transform(theta, d, a, alpha)
    T = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
         sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
         0, sin(alpha), cos(alpha), d;
         0, 0, 0, 1];
end

