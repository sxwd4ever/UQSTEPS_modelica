import numpy as np
import scipy.optimize

# ===============================================
# SETUP: define common compoments of the problem


def our_function(coeff, data):
    """
    The function we care to optimize.

    Args:
        coeff (np.ndarray): are the parameters that we care to optimize.
        data (np.ndarray): the input data
    """
    A, B, C = coeff
    x, y, z = data.T
    return (x - A + y - B) / 2 + np.sqrt(((x - A - y + B) / 2) ** 2 + C * z ** 2)


# ===============================================
# FORMULATION #1: a general minimization problem

# Here the bounds and error are all specified within the general objective function
def general_objective(coeff, data, target):
    """
    General function that simply returns a value to be minimized.
    The coeff will be modified to minimize whatever the output of this function
    may be.
    """
    # Constraints to keep coeff above 0
    if np.any(coeff < 0):
        # If any constraint is violated return infinity
        return np.inf
    # The function we care about
    prediction = our_function(coeff, data)
    # (optional) L2 regularization to keep coeff small
    # (optional) reg_amount = 0.0
    # (optional) reg = reg_amount * np.sqrt((coeff ** 2).sum())
    losses = (prediction - target) ** 2
    # (optional) losses += reg
    # Return the average squared error
    loss = losses.sum()
    return loss

# ===============================================
# FORMULATION #2: a special least squares problem

# Here all that is needeed is a function that computes the vector of residuals
# the optimization function takes care of the rest
def least_squares_residuals(coeff, data, target):
    """
    Function that returns the vector of residuals between the predicted values
    and the target value. Here we want each predicted value to be close to zero
    """
    # A, B, C = coeff
    # x, y, z = data.T
    prediction = np.matmul(coeff, data.T) # our_function(coeff, data)
    vector_of_residuals = (prediction - target.T) ** 2

    err = vector_of_residuals.sum()
    # if trial_count % 50 == 0:
    #     print(f'err={err}\n')

    # trial_count += 1
    return err

# Define some training data
data = np.array([
    [1, 0.445103675, -0.319950065],
    [1, 0.467270353, -0.171841875],
    # [1, 0.470722678, -0.28790764],
    [1, 0.370673736, 0.032916167]
])

# Define training target
# This is what we want the target function to be equal to
# target = 0

# target = np.array([-0.020626504,
# 0.662247263,
# 0.214637263,
# 0.751144649,
# ])


target = np.array([
    -0.008957977,
    0.287610332,
    # 0.093215779,
    0.326217976
])

# Make an initial guess as to the parameters
# either a constant or random guess is typically fine
num_coeff = 3
coeff_0 = np.ones(num_coeff)
# coeff_0 = np.random.rand(num_coeff)

# general_result = scipy.optimize.minimize(general_objective, coeff_0,
#                                          method='Nelder-Mead',
#                                          args=(data, target))
# # Test what the squared error of the returned result is
# coeff = general_result.x
# general_output = our_function(coeff, data)
# print('====================')
# print('general_result =\n%s' % (general_result,))
# print('---------------------')
# print('general_output = %r' % (general_output,))
# print('====================')

# Here the bounds are specified in the optimization call
bound_gt = np.full(shape=num_coeff, fill_value=-10, dtype=np.float)
bound_lt = np.full(shape=num_coeff, fill_value=10, dtype=np.float)
bounds = (bound_gt, bound_lt)

trial_count:int = 0

lst_sqrs_result = scipy.optimize.least_squares(least_squares_residuals, coeff_0,
                                               args=(data, target), bounds=bounds, ftol=5e-3)
# Test what the squared error of the returned result is
coeff = lst_sqrs_result.x
lst_sqrs_output = least_squares_residuals(coeff, data, target)
print('====================')
print('lst_sqrs_result =\n%s' % (lst_sqrs_result,))
print('---------------------')
print('lst_sqrs_output = %r' % (lst_sqrs_output,))
print('====================')