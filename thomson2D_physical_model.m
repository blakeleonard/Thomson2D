% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis

% Physical Model (Euler Method) of 2D Thomson Problem System with radius one
         

clear;


% Generate Intial Random Distribution of Charges on disk

% Initial conditions

N = 17;      % number of charges

time_0 = cputime;

U_min = 10000;

anneal_flag = 1;


% Generate N * 2 matrix of charge positions

polarposition = rand (N, 2 );

for istep = 1:N

    polarposition (istep, 1) = sqrt (polarposition(istep, 1));

    polarposition (istep, 2) = polarposition (istep, 2) * 6.2832;

    rplot(istep) = polarposition (istep, 1);

    thetaplot(istep) = polarposition (istep, 2);

end


% Plot Initial Configuration

axis([-1, 1, -1, 1]);

polar(thetaplot, rplot, '+');


% Convert from Polar to Cartesian

position = zeros(N, 2);

for istep = 1:N


    position(istep, 1 ) = polarposition(istep, 1) * cos(polarposition(istep, 2));

    position(istep, 2 ) = polarposition(istep, 1) * sin(polarposition(istep, 2));

end


% Perform Physical Coulomb Simulation

% Initial Conditions and Variable Declaration

mass = 1;

k = 0.05;

q = 1;

tau = 0.001;

plotstep = 100;

maxstep = 100000;

prevx = [0 0 0]; prevy = [0 0 0]; prevt = [0 0 0];

velocity = zeros(N,2);

e = 0.5;   % velocity change factor

A = 0.85;   % Annealing factor

anneal_count = 0;

anneal_energy = 136.3;

anneal_step = 15;

v_mag_coll = 0.1;


% Loop through max steps (Euler Method)

for hstep = 1:maxstep

    for istep = 1:plotstep

        for kstep = 1:N

            r = [position(kstep, 1), position(kstep, 2)];

            v = [velocity(kstep, 1), velocity(kstep, 2)];

            t = (istep-1)*tau;                                     % Current Time


            %* Calculate the acceleration of the particle due to the force from the other charges

            accel = [0, 0];

            for jstep = 1:N

                if ( jstep ~= kstep )

                    src_r = [position(jstep, 1), position(jstep, 2)];

                    d = sqrt( (r(1) - src_r(1))^2 + (r(2) - src_r(2))^2 );

                    accel = accel + ( ( ( k * q^2 ) / ( d ^ 3 ) ) * ( r - src_r) );

                end

            end

            accel = accel / mass;    % Convert force to acceleration


            %* Calculate the new position and velocity using Euler method

            r = r + tau*v;

            v = v + tau*accel;


            % Impose Boundary Conditions

            d = sqrt( r(1)^2 + r(2)^2 );

            if ( d > 1 )

                r = r/d;

                v_mag = sqrt(v(1)^2 + v(2)^2);

                if (v_mag > v_mag_coll );

                    v(1) =  - e* v(1);     % Play around with altering velocity ( elastic collision? )

                    v(2) =  - e * v(2);

                end

            end


            % Record Position and Velocity

            position(kstep, 1) = r(1);

            position(kstep, 2) = r(2);

            velocity(kstep, 1) = v(1);

            velocity(kstep, 2) = v(2);


            if (istep == plotstep)

                xplot(kstep) = position(kstep, 1);

                yplot(kstep) = position(kstep, 2);

            end

        end

    end


    % Plot Charges

    drawnow;

    plot(xplot, yplot, '+');

    axis([-1, 1, -1, 1]);


    % Compute Energy in units of q^2 / (4*Pi*Eps0)

    U = 0;

    for istep=1:N

        for jstep = (istep+1):N

            Q = sqrt( ( position(istep, 1) - position(jstep, 1) ) ^ 2 + ( position(istep, 2) - position(jstep, 2) ) ^ 2 );

            U = U + Q ^ (-1);

        end

    end


    if (U < U_min)

        U_min = U;

        fprintf('Best Energy: %g \n', U_min);

        disp('');

        fprintf('Time Elapsed: %g \n', cputime - time_0);

        disp('');

    end


    if ( U_min < anneal_energy & anneal_flag == 1 )    %  Annealing

        if (anneal_count == anneal_step)

            fprintf('Annealing... \n');

            for kstep = 1:N

                for lstep = 1:2

                    velocity(kstep,lstep) = A * velocity(kstep, lstep);

                end

            end

            anneal_count = 0;

        end

        anneal_count = anneal_count + 1;

    end

end

disp('');

z = input('Hit Enter to Advance');

disp('');