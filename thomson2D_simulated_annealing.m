% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis

% Simulated Annealing algorithm applied to 2D Thomson Problem with radius one
         

clear;


% Cycle through copies of Annealing Algorithm on different initial configurations to find lowest energy

globalstep = 10;

orig_time_0 = cputime;


for nstep=1:globalstep

	time_0 = cputime;


	% Generate Intial Random Distribution of Charges on disk with radius 1

	% Initial conditions

	N = 17;      % number of charges	
	

	% Generate N * 2 matrix of charge positions 

	polarposition = rand (N, 2 );

	for istep = 1:N	

		polarposition (istep, 1) = sqrt (polarposition(istep, 1));		

		polarposition (istep, 2) = polarposition (istep, 2) * 6.2832; 

	end


	% Convert from Polar to Cartesian

	position = zeros(N, 2);
		
	for istep = 1:N

		position(istep, 1 ) = polarposition(istep, 1) * cos(polarposition(istep, 2));

		initial_xplot(istep) = position(istep, 1);

		xplot(istep) = position(istep, 1);
		
		best_xplot(istep) = position(istep, 1);

		position(istep, 2 ) = polarposition(istep, 1) * sin(polarposition(istep, 2));

		initial_yplot(istep) = position(istep, 2);

		yplot(istep) = position(istep, 2);

		best_yplot(istep) = position(istep, 2);

	end


	% Compute Energy in units of q^2 / (4*Pi*Eps0)

	U_old = 0;		

	for kstep=1:N

		for lstep = (kstep+1):N
			
			Q = sqrt( ( position(kstep, 1) - position(lstep, 1) ) ^ 2 + ( position(kstep, 2) - position(lstep, 2) ) ^ 2 );				

			U_old = U_old + Q ^ (-1);

		end

	end


	U_min = U_old;


	%Plot Initial Configuration

	drawnow;

	plot(xplot, yplot, '+');
		
	fprintf('Initial Energy: %g \n', U_old);
	
	disp('');
	

	% Simulated Annealing Algorithm

	% Initial Conditions and Variable Declaration

	beta = N / U_old;

	F_beta = 1.15;

	ACCEP = 0.8;

	ACCEP_temp = 0.5;

	IEQ = 10;

	ISAMP = 20;

	epsilon = 0.00001;

	randstep_max = 0.9;

	F_delta = 0.9;                                        

	maxstep = 200;

	plotstep = IEQ*N;

	plotstep2 = IEQ*N + ISAMP;


	% Determine initial beta for Annealing

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

		end

		
		if ( U < U_min )

			U_min = U;

			for mstep = 1:N

				best_xplot(mstep) = xplot(mstep);

				best_yplot(mstep) = yplot(mstep);

			end

		end

			
		% Plot Charges

		drawnow;

		plot(best_xplot, best_yplot, '+');

		fprintf('Best Energy: %g \n', U_min);
	
		disp('');


		% If acceptance rate is less than ACCEP, reduce beta and repeat, otherwise break and begin tempering

		if ( ( change_count / plotstep ) < ACCEP )

			beta = beta / 2;

		else
						
			fprintf('Appropriate Acceptance Level Achieved, Begin Tempering with beta = %g \n', beta);
	
			disp('');

			initial_metro_steps = istep * plotstep;

			break;

		end

	end
		

	% Annealing Algorithm

	change_count = 0;

	for istep = 1:maxstep


		% Increase beta by F_beta factor and perform Metropolis

		beta = F_beta*beta;
		
		disp('');
		
		fprintf('BigStep %g: Beta increased by factor of %g to beta = %g \n', istep, F_beta, beta);

		disp('');
			
		if ( istep ~= 1 & ( change_count / (plotstep + ISAMP ) ) < ACCEP_temp)   

			randstep_max = randstep_max * F_delta;
			
			fprintf('Reducing max step by factor of %g to maxstep = %g\n', F_delta, randstep_max);
	
			disp('');

		end

		change_count = 0;

		U_tot = 0;

		for jstep = 1:plotstep2


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


			% Check if Energy is minimum and archive

			if ( U < U_min )

				U_min = U;
                
				fprintf('Best Energy: %g \n', U_min);
	
				disp('');
        
				fprintf('Time Elapsed: %g \n', cputime - time_0);
        
				disp('');

				for mstep = 1:N

					best_xplot(mstep) = xplot(mstep);

					best_yplot(mstep) = yplot(mstep);

				end

			end

			
			if (jstep > plotstep )

				U_tot = U_tot + U;

			end

		end


		% Plot Charges

		drawnow;

		plot(best_xplot, best_yplot, '+');


		% If avg energy being sampled is close enough to minimum recorded break out and start over with new initial config

		U_avg = U_tot / ISAMP;

		metro_steps = istep*plotstep2 + initial_metro_steps;

		%fprintf('Avg Energy after %g seconds: %g \n', cputime - time_0, U_avg);

		if (  ( abs( U_avg - U_min) ) /U_min < epsilon )

			disp('');			

			%fprintf('Avg close to minimum');

			disp('');

			break;

		end

	end

	fprintf('Run%g Energy: %g \n', nstep, U_min);
    
	fprintf('Run%g time: %g \n', nstep, cputime - time_0);

	disp('');


	run_energy(nstep) = U_min;
    
	run_time(nstep) = cputime - time_0;
	
end


for istep= 1:globalstep

	fprintf('Run %g Energy: %g \n', istep, run_energy(istep));

	fprintf('at %g seconds\n', run_time(istep) );

	disp('');

end

fprintf('Total time: %g seconds\n', cputime - orig_time_0 );