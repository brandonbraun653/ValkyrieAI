The goal of this test is to see how well the RNG object handles being accessed by a huge
number of threads (>10000) and serving requests for random numbers. In this case, I am 
trying to cause a crash. Due to the distribution test results in the mersenne_twister folder, it 
is apparent that creating one or two RNG engines per AI tuner is more important than spawning
hundreds of unique engines. This should guarantee better randomness in the AI software.