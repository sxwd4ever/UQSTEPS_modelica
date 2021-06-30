
from random import Random, random, randrange

# function to optimize: takes in a list of decision variables, returns an objective value
# this is the Rosenbrock function: http://en.wikipedia.org/wiki/Rosenbrock_function
# the global minimum is located at x = (1,1) where f(x) = 0


def my_function(x):
  return (1-x[0])**2 + 100*(x[1] - x[0]**2)**2

# function to perform (a very crude, stupid) optimization
# bounds = lower and upper bounds for each decision variable (2D list)
# NFE = number of function evaluations to perform
# f = the function to be optimized


def optimize(bounds, NFE, f, eps=1e-3):
  D = len(bounds)  # number of decision variables
  best_f = 9999.0  # initialize the "best found" - both the function value and the x values
  best_x = [None]*D
  i = 0

  for i in range(NFE):
    # use an "operator" to generate a new candidate solution
    # this is "uniform mutation" in MOEA lingo
    new_x = [bounds[d][0] + random()*(bounds[d][1] - bounds[d][0])
             for d in range(D)]
    new_f = f(new_x)
    if new_f < best_f:  # see if it's an improvement -- in multiobjective, this is the Pareto sort
      err = abs(new_f - best_f)

      best_f = new_f
      best_x = new_x

      if(err < eps):
        break

  return {'best_x': best_x, 'best_f': best_f, 'i': i}


def optimize_v2(bounds, NFE, f, eps=1e-3):
    D = len(bounds[0])  # number of decision variables
    best_f = 9999.0  # initialize the "best found" - both the function value and the x values
    best_x = [None]*D
    i = 0
    
    # 0 left, 1 right
    x_lb = [bounds[0][d]  for d in range(D)]
    x_rb = [bounds[1][d]  for d in range(D)]
    x_best = x_lb

    f_lb = f(x_lb)
    f_rb = f(x_rb)
    if (f_lb > f_rb):
        f_best = f_rb
    else:
        f_best = f_lb

    for i in range(NFE):
        # use an "operator" to generate a new candidate solution
        # this is "uniform mutation" in MOEA lingo
        x_new = [ x_lb[d] + random() * (x_rb[d] - x_lb[d]) for d in range(D)]
        f_new = f(x_new)


        # see if it's an improvement -- in multiobjective, this is the Pareto sort
        if f_new < f_best: # shrink right bound
            err = abs(f_new - f_best)

            f_best = f_new # update f_best
            x_best = x_new

            if f_lb > f_rb:
                x_lb = x_new
                f_lb = f_new
            else:
                x_rb = x_new
                f_rb = f_new    
            
            if(err < eps):
                break


    return {'x_best': x_best, 'f_best': f_best, 'i': i}

def dist(x,y):
    D = len(x)
    r = [(x[d] - y[d]) ** 2 for d in range(D)]

    r = sum(r) **.5

    return r

def cal_r_max(x, bounds):
    D = len(x)  # number of decision variables
    # 0 left, 1 right
    x_lb = [bounds[0][d]  for d in range(D)]
    x_rb = [bounds[1][d]  for d in range(D)]    

    r_dist_l = dist(x, bounds[0])
    r_dist_r = dist(x, bounds[1])

    return min([r_dist_l, r_dist_r])
    
def inbound(x, bounds):
    D = len(x)

    for i in range(D):
        if x[i] < bounds[0][i] or x[i] > bounds[1][i]:
            return False

    return True


def optimize_v3(bounds, NFE, f, eps=1e-3):
    D = len(bounds)  # number of decision variables
    best_f = 9999.0  # initialize the "best found" - both the function value and the x values
    best_x = [None]*D
    i = 0
    
    # 0 left, 1 right
    x_lb = [bounds[0][x]  for x in range(D)]
    x_rb = [bounds[1][x]  for x in range(D)]
    
    x_new = [(x_lb[d] + x_rb[d]) / 2 for d in range(D)]
    r_max = cal_r_max(x_new, bounds)

    f_lb = f(x_lb)
    f_rb = f(x_rb)
    if (f_lb > f_rb):
        f_best = f_rb
        x_best = x_rb
    else:
        f_best = f_lb
        x_best = x_lb

    for i in range(NFE):
        # use an "operator" to generate a new candidate solution
        # this is "uniform mutation" in MOEA lingo
        r_max_tmp = cal_r_max(x_new, bounds)
        if(r_max_tmp <= r_max):
            r_max = r_max_tmp

        for j in range(10):
            
            import random as rnd

            x_new_tmp = [ x_new[d] + rnd.uniform(-1,1) * r_max for d in range(D)]

            if(inbound(x_new_tmp, bounds)): # find a feasible x_new 
                x_new = x_new_tmp
                break
            else:
                continue        

        f_new = f(x_new)

        # see if it's an improvement -- in multiobjective, this is the Pareto sort
        if f_new < f_best: # shrink right bound
            err = abs(f_new - f_best)

            f_best = f_new # update f_best
            x_best = x_new

            if f_lb > f_rb:
                x_lb = x_new
                f_lb = f_new
            else:
                x_rb = x_new
                f_rb = f_new    
            
            if(err < eps):
                break


    return {'x_best': x_best, 'f_best': f_best, 'i': i}

# now let's try it...
# (the Rosenbrock problem technically doesn't have "bounds", but we'll make some up..)
eps = 5e-3
bounds = [([-1, 5]), [-1, 5]]
result = optimize(bounds, 10000, my_function, eps)
print(result)  # it should be near best_f = 0.0 and best_x = [1,1], hopefully

bounds = [(-1,-1), (5,5)]
result = optimize_v2(bounds, 10000, my_function, eps)
print(result)  # it should be near best_f = 0.0 and best_x = [1,1], hopefully

bounds = [(-1,-1), (5,5)]
result = optimize_v3(bounds, 10000, my_function, eps)
print(result)  # it should be near best_f = 0.0 and best_x = [1,1], hopefully
