function cylindrical

R_max = 20
stepR = 0.2
steps = 50

[R,THETA] = meshgrid(0:stepR:R_max, 0:6.28/steps:6.28);

% enter Z as a function of R and THETA here
Z = sin(R)./R

X = R.*cos(THETA)
Y = R.*sin(THETA)

surf(X,Y,Z)
 % cylindrical coordinates

end