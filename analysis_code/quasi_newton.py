#----------------------------------------------------------------------------#
#   My Quasi-Newton Methods
#   
#
#   Last modified: Mon 14 Apr 15:53:15 2014
#----------------------------------------------------------------------------#

"""
My Quasi-Newton methods module. Only for when you have an analytic expression
for the jacobian.
WARNING: in Line search the quadratic model looks like it works, but the cubic
    model is untested. I am pretty sure it gives a lambda between 0 and 1, but
    if it gives funny behaviour probably best to output what is happening.
"""



# MODULES
from scipy import *
from scipy import linalg

def newton(solve_eq, find_jac, x, NRtol=10e-10):
    """
    My standard Newton method for use with a function which finds the jacobian.
    """

    finCond = False
    while not finCond:
        F_x0 = solve_eq(x)
        # Print norm and check if we can exit yet.
        L2 = linalg.norm(F_x0)
        print """|F_xn| = {0:10.5g}""".format(L2)
        if L2 < NRtol:
            print "Solution found!"
            break 

        # Calculate the newton step, dx, and take it
        j_x0 = find_jac(x)
        dx = linalg.solve(j_x0, -F_x0)
        x = x + dx
    
    return x 

def line_search(solve_eq, find_jac, x, alpha=1e-4, NRtol=10e-10):
    """
    My implementation of the line search global quasi Newton Raphson method for
    when you have a function to calculate the analytic jacobian. If you don't
    have this just use the scipy Newton Krylov solver (has a finite difference
    approximation to the jacobian built in).  This method only takes a Newton
    step if this decreases the residual squared.  This function is known as the
    master function f in this code. 

    Parameters:
        alpha: sets the size of average rate of decrease of f should be at least
               this fraction of the initial rate of decrease of f. Never greater
               than 1. Ideally as low as possible. Numerical recipies thinks
               1e-4 is fine.
        NRtol: The stopping criteria on the L2 norm of the residual.
    
    WARNING: Could clash with global variables if you have been silly enough to
    choose rubbish names.
    """

    alpha = 1e-4
    NRtol = 10e-10

    finCond = False
    while not finCond:
        # Calculate the newton step, dx
        F_x0 = solve_eq(x)
        j_x0 = find_jac(x)
        dx = linalg.solve(j_x0, -F_x0)

        # Define the master function
        f_x0 = 0.5*dot(conj(F_x0), F_x0)

        # Decide whether to take the Newton Step by Armijo line search method
        # First initialise variables so that first iteration happens 
        lam = 1 
        lamPrev = 1  
        slope_x0dx = -dot(conj(F_x0), F_x0)
        f_xn = f_x0 + alpha*lam*slope_x0dx + 1
        f_lam2Prev = 0 # Doesn't matter, will be set before it is used
        f_lamPrev = 0 
        counter = 0

        # Now choose a lambda and see if it is good.
        while f_xn >= f_x0 + alpha*lam*slope_x0dx:

            if counter == 1:
                # set lambda by a quadratic model for the residual master function f
                lam = - slope_x0dx / 2*(f_xn - f_x0 - slope_x0dx)
                
                # impose upper and lower bounds on lambda 
                if lam > 0.5:
                    lam = 0.5
                if lam < 0.1:
                    lam = 0.1
                print "square model lambda =", lam

            elif counter > 1:
                # set lambda by a cubic model for the residual master function f
                abmat = zeros((2,2))
                abmat[0,0] = 1/(lamPrev*lamPrev)
                abmat[0,1] = -1/(lam2Prev*lam2Prev)
                abmat[1,0] = -lam2Prev/(lamPrev*lamPrev)
                abmat[1,1] = lamPrev/(lam2Prev*lam2Prev)

                f3vec = zeros(2)
                f3vec[0] = f_lamPrev - f_x0 - slope_x0dx*lamPrev
                f3vec[1] = f_lam2Prev - f_x0 - slope_x0dx*lam2Prev

                abvec = (1./(lamPrev-lam2Prev)) * dot(abmat, f3vec)
                aaa = abvec[0]
                bbb = abvec[1]
                lam = (- bbb + sqrt(bbb**2 - 3*aaa*slope_x0dx)) / 3*aaa

                # impose upper and lower bounds on lambda 
                if lam > 0.5*lamPrev:
                    lam = 0.5*lamPrev
                if lam < 0.1*lamPrev:
                    lam = 0.1*lamPrev
                print "cubic model lambda", lam

            # calculate the residual and master function so we can see if the
            # step was a good one.
            F_xn = solve_eq(x + lam*dx)
            f_xn = 0.5*dot(conj(F_xn), F_xn)
            #print """   |F_xn| = """, linalg.norm(F_xn) 

            # update old values for cubic method
            lam2Prev = lamPrev
            lamPrev = lam
            f_lam2Prev = f_lamPrev
            f_lamPrev = f_xn

            counter += 1
        
        print " loop counter of last step = ", counter-1
        print " lambda of the step = ", lam
        # change x to the value at the step we just took
        x = x + lam*dx

        # Print norm and check if we can exit yet.
        L2 = linalg.norm(F_xn)
        print """|F_xn| = {0:10.5g}""".format(L2)
        if L2 < NRtol:
            print "Solution found!"
            finCond = True

    return x

def main():
    print 'Ideally, I guess I would write some simple unit tests to put here'
    print """
    Be warned, in Line search the quadratic model looks like it works, but the cubic
    model is untested. I am pretty sure it gives a lambda between 0 and 1, but
    if it gives funny behaviour probably best to output what is happening.
    """

if __name__=='__main__':
    main()
