% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis

% Interacting Replica technique applied to 2D Thomson Problem with radius one         


clear;

globalstep = 10;

for zstep = 1:globalstep


    % Initial conditions & Variable Declaration

    time_0 = cputime;

	N = 17;          % number of charges	
	
	R = 20;          % number of replicas

	U_min = 10000;
	
	randstep_max = 0.9;
                                               
	maxstep = 100000;

	plotstep = 20;	

	infostep = 5;
    
    metro_steps = 20;
	
	change_count = 0;

	beta_max = 200;

	beta_stddev = 25;

	info_on_step = 1;
    
    end_time = 35;

	
	% Generate distribution of betas

	for rstep = 1:R		

		beta(rstep) = beta_max * exp ( - ( ( rstep-(R/2) ) ^ 2 ) / ( 2 * beta_stddev ) );

		fprintf('Beta(R=%g) = %g\n', rstep, beta(rstep));
	
		disp('');

	end

	
    % Generate Intial Random Distribution of Charges on Disk for all replicas

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


    % Perform Global Free Energy Alogorithm

	for istep = 1:maxstep	

		for jstep = 1:plotstep


			% Energy Minimization Step ( metro_steps on each replica per infostep on all replicas)
			
            for ystep = 1:metro_steps

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
                        
                        last = cputime - time_0;

                        for qstep = 1:N

                            best_xplot(qstep) = xposition(rstep, qstep);

                            best_yplot(qstep) = yposition(rstep, qstep);

                        end

                    end

                end

            end


			% Maximize Mutual Information Step


			% Perform Just metropolis at first, info step turning on after
			% info_on_step's

			if ( istep >= info_on_step)	
                
				if (istep == info_on_step && jstep == 1)

					fprintf('Beginning step involving info max half-step');

					disp('');
			
					disp('');

				end


				for mstep = 1:infostep


					% Select a random charge in a random replica
		
					randchvect = randperm(N);

					randcharge = randchvect(1);
			
					randrepvect = randperm(R);
		
					randreplica = randrepvect(1);

					
                    for qstep = 1:R-1

						distance_min(qstep) = 1;

					end


					% Determine Corresponding charges in all other replicas
					% by finding charge closest to random charges position

					for qstep = 1:R-1

						for ostep = 1:N

							distance = sqrt( ( xposition(randreplica, randcharge) - xposition(randrepvect(1+qstep), ostep) )^2 + ( yposition(randreplica, randcharge) - yposition(randrepvect(1+qstep), ostep) )^2 );

							if distance < distance_min(qstep)
                                
                                corres_charge(randrepvect(1+qstep)) = ostep;

								distance_min(qstep) = distance;

							end

						end

                    end


                    % Calculate the average position of corresponding
                    % charges

                    xpos_tot = 0;
                    
                    ypos_tot = 0;
                    
                    for qstep = 1:R-1
                        
                        xpos_tot = xpos_tot + xposition(randrepvect(1+qstep), corres_charge(randrepvect(1+qstep)));
                        
                        ypos_tot = ypos_tot + yposition(randrepvect(1+qstep), corres_charge(randrepvect(1+qstep)));
                        
                    end
                    
                    xpos_avg = xpos_tot / ( R - 1 );
                    
                    ypos_avg = ypos_tot / ( R - 1 );
                    
                    
                    % Move random charge to avg position

					xposition(randreplica, randcharge) = xpos_avg;

					yposition(randreplica, randcharge) = ypos_avg;


					% Impose Boundary Conditions

					d = sqrt( xposition(randreplica, randcharge)^2 + yposition(randreplica, randcharge)^2 );			

					if ( d > 1 )

						xposition(randreplica, randcharge) = xposition(randreplica, randcharge) / d;

						yposition(randreplica, randcharge) = yposition(randreplica, randcharge) / d;

					end


					% Compute Energy in units of q^2 / (4*Pi*Eps0)         

					U(randreplica) = 0;		

					for kstep=1:N

						for lstep = (kstep+1):N
			
							Q = sqrt( ( xposition(randreplica, kstep) - xposition(randreplica, lstep) ) ^ 2 + ( yposition(randreplica, kstep) - yposition(randreplica, lstep) ) ^ 2 );				
	
							U(randreplica) = U(randreplica) + Q ^ (-1);
	
						end

					end

					
					% Check if Energy is minimum and archive

					if U(randreplica) < U_min

						U_min = U(randreplica);
                        
                        fprintf('Best Energy: %g \n', U_min);
	
                        disp('');
        
                        fprintf('Time Elapsed: %g \n', cputime - time_0);
        
                        disp('');
                        
                        last = cputime - time_0;
						
                        for qstep = 1:N

							best_xplot(qstep) = xposition(randreplica, qstep);

							best_yplot(qstep) = yposition(randreplica, qstep);

						end

					end


					U_old(randreplica) = U(randreplica);

                end
			
			end

		end


		% Plot Charges

        drawnow;

		plot(best_xplot, best_yplot, '+');
        
        time = cputime - time_0;
        
        if ( time > end_time )
            
            run_time(zstep) = last;
            
            run_energy(zstep) = U_min;
            
            fprintf('Run  Energy: %g \n', U_min);
	
            disp('');
            
            fprintf('Run Time: %g \n', last);
	
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