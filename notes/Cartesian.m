function Cartesian_3D

x_min = -5
x_max = 5
stepX = 0.25

y_min = -5
y_max = 5
stepY = 0.25

% plot set {X, Y, Z} in 3D

[X,Y] = meshgrid(x_min:stepX:x_max, y_min:stepY:y_max);

% enter your function here
Z = X.^2 - Y.^2

surf(X,Y,Z)

end

