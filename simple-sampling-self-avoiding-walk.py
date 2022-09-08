from random import seed
from random import random
import math
import timeit

start = timeit.default_timer()

# Seed the randomness so you get replicable results
seed(1)

# Number of steps to take
maxn = 40

# Number of iterations
samples = 1000000

# Initialise the board for our little walk
grid = [[False for i in range(2*maxn+1)] for j in range(2*maxn+1)]

def simplesamplingSAW(x,y,s):
    global iti, grid, n
    # If we have reached the maximum number of steps without crashing into our past-track, count this as a success
    if s==n:
        iti+=1
    else:
        grid[x][y] = True
        # Randomise a number between 0 and 1
        roll = random()
        # Keep heading in a random direction (North, East, South, West) if we haven't walked there yet
        if 0 <= roll < 0.25:
            if not grid[x][y+1]:
                simplesamplingSAW(x,y+1,s+1)
        elif 0.25 <= roll < 0.5:
            if not grid[x+1][y]:
                simplesamplingSAW(x+1,y,s+1)
        elif 0.5 <= roll < 0.75:
            if not grid[x][y-1]:
                simplesamplingSAW(x,y-1,s+1)
        else:
            if not grid[x-1][y]:
                simplesamplingSAW(x-1,y,s+1)
        # If we can't proceed in the randomised direction, we will backtrack and try again
        grid[x][y] = False

# Estimates of unique walks for walks from 1 to n steps:
for n in range(1,maxn+1):
    iti = 0
    # Try the walk "samples" of times
    for _ in range(1,samples):
        # Take 1st step from center of board
        simplesamplingSAW(n,n,0)
    
    # Extrapolate out to 4^n possible cases to estimate the number of self-avoiding walks
    esti = 4 ** n * iti/samples
    sd = math.sqrt(4 ** n * iti/samples * (1-iti/samples))
    print('Steps taken: {}; Estimated number of unique walks: {:.0f}; Standard deviation of estimate: {:.4f}'.format(n, esti, sd))

stop = timeit.default_timer()

print('Time: ', stop - start)