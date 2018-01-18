% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis

% Replica Exchange Technique applied to 2D Thomson Problem with radius one
         

clear;


% Generate Intial Random Distribution of Charges on disk for all replicas

% Initial conditions

N = 17;          % number of charges

R = 20;          % number of replicas

U_min = 10000;

time_0 = cputime;


% Generate Random Matrix for initial poisitions

rposition = rand (R, N);

thetaposition = rand (R, N);

xposition = zeros(R, N);

yposition = zeros(R, N);


for rstep = 1:R


    % Generate R * N * 2 array of charge positions

    for istep = 1:N

        rposition (rstep, istep) = sqrt (rposition(rstep, istep));

        thetaposition (rstep, istep) = thetaposition (rstep, istep) * 6.2832;

    end


    % Convert from Polar to Cartesian

    for istep = 1:N

        xposition(rstep, istep) = rposition(rstep, istep) * cos(thetaposition(rstep, istep));

        yposition(rstep, istep) = rposition(rstep, istep) * sin(thetaposition(rstep, istep));

    end


    % Compute Energy in units of q^2 / (4*Pi*Eps0)

    U_old(rstep) = 0;

    for istep=1:N

        for jstep = (istep+1):N

            Q = sqrt( ( xposition(rstep, istep) - xposition(rstep, jstep) ) ^ 2 + ( yposition(rstep, istep) - yposition(rstep, jstep) ) ^ 2 );

            U_old(rstep) = U_old(rstep) + Q ^ (-1);

        end

    end


    fprintf('Initial Energy of replica #%g: %g \n', rstep, U_old(rstep));

    disp('');

end


% Perform Metropolis

% Initial Conditions and Variable Declaration

randstep_max = 0.9;

maxstep = 100000;

exchstep = 10

beta_max = 200;

sigma = 25;


% Generate distribution of betas

for rstep = 1:R

    beta(rstep) = beta_max * exp ( (- ( rstep-(R/2) ) ^2 )/(2*sigma));

end


% Metropolis Algorithm

for istep = 1:maxstep

    change_count = 0;

    for jstep = 1:exchstep

        for rstep = 1:R


            % Move a random charge by a random amount in a random direction

            randchvect = randperm(N);

            randcharge = randchvect(1);


            rand_xstep = randstep_max * rand(1);

            if ( rand(1) < 0.5 )

                rand_xstep = rand_xstep * -1;

            end


            rand_ystep = randstep_max * rand(1);

            if ( rand(1) < 0.5 )

                rand_ystep = rand_ystep * -1;

            end


            oldx = xposition(rstep, randcharge);

            oldy = yposition(rstep, randcharge);


            xposition(rstep, randcharge) = xposition(rstep, randcharge) + rand_xstep;

            yposition(rstep, randcharge) = yposition(rstep, randcharge) + rand_ystep;


            % Impose Boundary Conditions

            d = sqrt( xposition(rstep, randcharge)^2 + yposition(rstep, randcharge)^2 );

            if ( d > 1 )

                xposition(rstep, randcharge) = xposition(rstep, randcharge) / d;

                yposition(rstep, randcharge) = yposition(rstep, randcharge) / d;

            end


            % Compute Energy in units of q^2 / (4*Pi*Eps0)

            U(rstep) = 0;

            for kstep=1:N

                for lstep = (kstep+1):N

                    Q = sqrt( ( xposition(rstep, kstep) - xposition(rstep, lstep) ) ^ 2 + ( yposition(rstep, kstep) - yposition(rstep, lstep) ) ^ 2 );

                    U(rstep) = U(rstep) + Q ^ (-1);

                end

            end


            % Compare Energies of new and old configurations & switch back if appropriate

            if ( ( exp( beta(rstep) * ( U_old(rstep) - U(rstep) ) ) ) > rand(1) )

                change_count = change_count + 1;

            else

                xposition(rstep, randcharge) = oldx;

                yposition(rstep, randcharge) = oldy;

                U(rstep) = U_old(rstep);

            end


            U_old(rstep) = U(rstep);


            % Check if Energy is minimum and archive

            if U(rstep) < U_min

                U_min = U(rstep);

                fprintf('Best Energy: %g \n', U_min);

                disp('');

                fprintf('Time Elapsed: %g \n', cputime - time_0);

                disp('');


                for qstep = 1:N

                    best_xplot(qstep) = xposition(rstep, qstep);

                    best_yplot(qstep) = yposition(rstep, qstep);

                end

            end

        end

    end


    % Swap beta values on replicas

    swapcount = 0;

    for rstep = 1:R-1

        if ( ( exp( ( beta(rstep) - beta(rstep+1) ) * ( U(rstep) - U(rstep+1) ) ) ) > rand(1) )

            betahold = beta(rstep);

            beta(rstep) = beta(rstep+1);

            beta(rstep+1) = betahold;

            swapcount = swapcount + 1;

        end

    end


    % Plot Charges

    drawnow;

    plot(best_xplot, best_yplot, '+');

    %fprintf('Swapped %g replicas on last run. \n', swapcount);

    %disp('');

end