% Загрузка данных
load('calib.mat');
Q = calib.Q;
T_actual = calib.T;

% Начальное приближение параметров
initial_params = [300, 0, pi/2, 0, 400, 0, 0, 500, 0, 0, 0, 0];

% Ограничения для параметров
lb = [290, -0.01, pi/2-0.01, -10, 390, -0.01, -10, 490, -0.01, -10, -0.01, -0.01];
ub = [310, 0.01, pi/2+0.01, 10, 410, 0.01, 10, 510, 0.01, 10, 0.01, 0.01];

% Опции для fmincon
options = optimset('TolX',1e-16,'MaxIter',2e3,'MaxFunEvals',2000);

% Выполнение оптимизации
[optimized_params, fval] = fmincon(@(params) objective_function(params, Q, T_actual), initial_params, [], [], [], [], lb, ub, [], options);

% Вывод оптимизированных параметров
disp('Оптимизированные параметры:');
disp(optimized_params);

function [d, R] = manipulator_model(q, params)
    % Извлечение параметров для каждого звена
    d1 = params(1);
    a1 = params(2);
    alpha1 = params(3);

    d2 = params(4);
    a2 = params(5);
    alpha2 = params(6);

    d3 = params(7);
    a3 = params(8);
    alpha3 = params(9);

    % Передача углов и параметров в DH преобразование
    T1 = dh_transform(q(1), d1, a1, alpha1);
    T2 = dh_transform(q(2), d2, a2, alpha2);
    T3 = dh_transform(q(3), d3, a3, alpha3);

    % Общая матрица преобразования
    T = T1 * T2 * T3;

    d = T(1:3, 4);
    R = T(1:3, 1:3);
end

function total_diff = objective_function(params, Q, T_actual)
    total_diff = 0;
    for i = 1:size(Q, 1)
        [d, R] = manipulator_model(Q(i, :), params);
        T_calculated = [R d; 0 0 0 1];

        % Вычисление нормы Фробениуса разности между матрицами
        diff = norm(T_calculated - T_actual(:, :, i), 'fro');
        total_diff = total_diff + diff;
    end
end

function T = dh_transform(theta, d, a, alpha)
    % T = [cos(theta), sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
    %      sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    %      0, sin(alpha), cos(alpha), d;
    %      0, 0, 0, 1];
    T = rotz(theta)*transz(d)*rotx(alpha)*transx(a);

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
