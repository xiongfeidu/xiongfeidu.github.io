function first_order_ode

    t0 = 0;
    y0 = 4;
    t_max = 4;
    trange = [t0 t_max];
    
    [t,y] = ode45( @rhs , trange, y0);

    plot(t,y)
    
    xlabel('t'); ylabel('y');
    title('First Order ODE');
    
    function f = rhs(t,y)
        f = -t/y;
    end

end