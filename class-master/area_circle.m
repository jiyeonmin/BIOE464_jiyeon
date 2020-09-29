N = 10000;             % number of random points
Nin = 0;
rx = 2*(rand(1,N)-0.5);   % vector of N random numbers uniform in [-1,1]
ry = 2*(rand(1,N)-0.5);   % vector of N random numbers uniform in [-1,1]

rxin  = []; ryin = [];
rxout = []; ryout = [];   % initialization all vectors are empty

for i = (1:N)
    rad2 = rx(i)^2 + ry(i)^2;
    if(rad2 < 1) 
	rxin = [rxin rx(i)]; 
	ryin = [ryin ry(i)]; 
        Nin = Nin + 1;    % update the number of points inside the circle
    else
	rxout = [rxout rx(i)]; 
	ryout = [ryout ry(i)]; 
    end
end

Area = 4*Nin/N;     % Nin/N is the fraction of points inside the circle, 4 the area of the embedding square
plot(rxin,ryin,'bx',rxout,ryout,'ro');  % Plotting points inside and outside the circle
					% 'bx'= blue with symbol x
                                        % 'ro'= red with symbol o
title(sprintf(['Estimated Area =%g with %g sampling points'], Area,N));
