% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis

% Metropolis Algorithm applied to 2D Thomson Problem with radius one
         

clear;

globalstep = 10;

for zstep = 1:globalstep


    % Generate Intial Random Distribution of Charges on disk

	% Initial conditions

	N = 17;      % number of charges	
	
    U_min = 10000;
    
    time_0 = cputime;
    
    end_time = 5;

	
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

		best_xplot(istep) = position(istep, 1);
        
        xplot(istep) = position(istep, 1);

		position(istep, 2 ) = polarposition(istep, 1) * sin(polarposition(istep, 2));

		best_yplot(istep) = position(istep, 2);
        
        yplot(istep) = position(istep, 2);

	end


	% Compute Energy in units of q^2 / (4*Pi*Eps0)

	U_old = 0;		

	for istep=1:N

		for jstep = (istep+1):N
			
			Q = sqrt( ( position(istep, 1) - position(jstep, 1) ) ^ 2 + ( position(istep, 2) - position(jstep, 2) ) ^ 2 );				

			U_old = U_old + Q ^ (-1);

		end

	end

	
	fprintf('Initial Energy: %g \n', U_old);
	
	disp('');


    % Perform Metropolis

	% Initial Conditions and Variable Declaration

	%beta = N / U_old;
    
    beta = 200;

	randstep_max = 0.9;
                                               
	maxstep = 10000;

	plotstep = 100;


	% Metropolis Algorithm

	for istep = 1:maxstep	

		change_count = 0;

		for jstep = 1:plotstep


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


			oldx = position(randcharge, 1 );

			oldy = position(randcharge, 2 );
			

			position(randcharge, 1 ) = position(randcharge, 1 ) + rand_xstep;

			position(randcharge, 2 ) = position(randcharge, 2 ) + rand_ystep;


			% Impose Boundary Conditions
			
			d = sqrt( position(randcharge, 1 )^2 + position(randcharge, 2 )^2 );			

			if ( d > 1 )

				position(randcharge, 1 ) = position(randcharge, 1 ) / d;

				position(randcharge, 2 ) = position(randcharge, 2 ) / d;

			end

		
			% Compute Energy in units of q^2 / (4*Pi*Eps0)

			U = 0;		

			for kstep=1:N

				for lstep = (kstep+1):N
			
					Q = sqrt( ( position(kstep, 1) - position(lstep, 1) ) ^ 2 + ( position(kstep, 2) - position(lstep, 2) ) ^ 2 );				
	
					U = U + Q ^ (-1);

				end

			end

		
			% Compare Energies of new and old configurations & switch back if appropriate

			if ( ( exp( beta * ( U_old - U ) ) ) > rand(1) )

				change_count = change_count + 1;

			else

				position(randcharge, 1 ) = oldx;

				position(randcharge, 2 ) = oldy;

				U = U_old;

            end

			U_old = U;
            
            xplot(randcharge) = position(randcharge, 1);     

			yplot(randcharge) = position(randcharge, 2);


            if ( U < U_min )

                U_min = U;
                
                fprintf('Best Energy: %g \n', U_min);
	
                disp('');
        
                fprintf('Time Elapsed: %g \n', cputime - time_0);
        
                disp('');
                
                last_time = cputime - time_0;
		
               	for mstep = 1:N

                    best_xplot(mstep) = xplot(mstep);

                    best_yplot(mstep) = yplot(mstep);

                end

            end

		end


		% Plot Charges
        
        drawnow;

		plot(best_xplot, best_yplot, '+');
        
        if ( (cputime - time_0) - last_time > end_time )
            
            run_time(zstep) = last_time;
            
            run_energy(zstep) = U_min;
            
            fprintf('Run  Energy: %g \n', U_min);
	
            disp('');
            
            fprintf('Run Time: %g \n', last_time);
	
            disp('');
            
            break;
            
        end

    end

end 


for istep = 1:globalstep
    
    fprintf('Run %g Energy: %g \n', istep, run_energy(istep));
	
    disp('');
    
    fprintf('Run %g Time: %g \n', istep, run_time(istep));
	
    disp('');
    
end