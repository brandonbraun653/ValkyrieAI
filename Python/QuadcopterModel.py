
import os
import numpy as np

import matlab.engine





if __name__ == '__main__':

    print('Starting matlab engine')
    eng = matlab.engine.start_matlab()
    eng.addpath(r'Matlab/')
    print('Done')


