function spherical

[THETA,PHI] = meshgrid(0:6.28/50:6.28, 0:3.14/50:3.14)

% enter R as a function of THETA and PHI
R = 2*cos(THETA)

X = R.*cos(THETA).*sin(PHI) % do not edit these
Y = R.*sin(THETA).*sin(PHI)
Z = R.*cos(PHI)

surf(X,Y,Z) % spherical coordinates
 
end
