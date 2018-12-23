function second_order_ode

t1 = 0;
t2 = 20;
step = 0.01;

%t = [t1 t2]; % either or of these two will work
t=t1:step:t2;

y_0      = 0;
yPrime_0 = 2;

[t,y]=ode45( @rhs, t, [y_0 yPrime_0] );

plot(t,y(:,1));

xlabel('t');
ylabel('y');
title('Second Order ODE');

% solve a*DDy + b*Dy + c*y = f(t)

    function ode = rhs(t,y)
        Dy = y(2); % do not change this line
        %DDy = -y(1) + square(t); % -c/a*y(1) - b/a*y(2) + f(t)
        DDy = -sin(y(1)); % pendulum motion
        %DDy = y(1)
        ode = [Dy; DDy];
    end
end