N =100000;          % number of random points
bins = 40;          % number of bins for the histogram (try bins=20)
rngs = rand(1,N);   % vector of N random numbers uniform in [0,1]
rg = [rngs(rngs<0.1) rngs(rngs>0.5)]; % select from rngs all elements smaller than 0.1 and larger than 0.5
                                      % [a b] concatenates the vectors a and b
subplot(1, 2, 1);
hist(rg, bins);     % plot the histogram of rg
title('1st method - Uniform in [0,0.1] U [0.5,1]')

rg = rngs*0.6+0.5;
rg = rg - floor(rg);
subplot(1, 2, 2);
hist(rg, bins);     % plot the histogram of rg
title('2nd method - Uniform in [0,0.1] U [0.5,1]')
