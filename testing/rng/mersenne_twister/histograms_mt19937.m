mtDat_continuous = csvread('mt19937_continuous.csv');
mtDat_continuous_reseed = csvread('mt19937_continuous_reseed.csv');
mtDat_reinit = csvread('mt19937_reinit.csv');
mtDat_reinit_reseed = csvread('mt19937_reinit_reseed.csv');

figure(1);
histogram(mtDat_continuous)
title('Continuous MT');

figure(2);
histogram(mtDat_reinit)
title('Reset MT');

figure(3);
histogram(mtDat_continuous_reseed)
title('Continuous MT Reseed/Reset');

figure(4);
histogram(mtDat_reinit_reseed)
title('Reset MT Reseed/Reset');