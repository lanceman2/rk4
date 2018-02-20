// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
//
// A C++ template class to solve ODEs (ordinary differential equations)
// via 4th order Runge–Kutta.
//

#ifndef __RK4__HPP__
#define __RK4__HPP__

#include <inttypes.h>
#include <stdlib.h>


// Pick your template floating point types:
// example: float, double, or long double
//
// X = degrees-of-freedom, R^N, phase space, state that evolves in time
// T = time
// STEP = time step
//
template <class X, class T, class STEP>
class RK4
{
    public:

        RK4(uint32_t n/*dimensions, number of degrees-of-freedom*/, STEP tStep);
        virtual ~RK4(void);

        void go(X *x/*user state*/, T from, T to/*to time*/);

    protected:

        // The users' function that computes derivatives for the set of n
        // first order ODEs.
        virtual void derivatives(T t, const X *x, X *xDot) = 0;

        STEP tStep; // This may need to change based on ODEs parameters.

    private:


        uint32_t n; // number of degrees of freedom, 2 for sine wave.

        X *k1, *k2, *k3, *k4, *k_2, *k_3, *k_4;
};

template <class X, class T, class STEP>
RK4<X, T, STEP>::RK4(uint32_t n_in, STEP tStep_in):
    tStep(tStep_in), n(n_in)
{
    k1 = (X *) malloc(sizeof(X)*n*7);
    if(!k1) throw "malloc() failed";
    k2 = k1 + n;
    k3 = k2 + n;
    k4 = k3 + n;
    k_2 = k4 + n;
    k_3 = k_2 + n;
    k_4 = k_3 + n;
}

template <class X, class T, class STEP>
RK4<X, T, STEP>::~RK4(void)
{
    if(k1)
    {
        free(k1);
        k1 = 0;
    }
}

// TODO: This interface needs work.  Sometimes we want an array of values
// out.
//
template <class X, class T, class STEP>
void RK4<X, T, STEP>::go(X *x, T t, T to)
{
    // This may not be the most efficient code but it may be easy to
    // follow.

    STEP dt = tStep; // The largest time step to start with.
    bool running = true;

    while(running)
    {
        if(t + dt > to)
        {
            // shorter time step, so we do not over shoot.
            dt = to - t;
            running = false;
        }

        uint32_t i;

        derivatives(t, x, k1);
        for(i=0; i<n; ++i)
        {
            k1[i] *= dt;
            k_2[i] = x[i] + k1[i]/2;
        }

        derivatives(t + dt/2, k_2, k2);
        for(i=0; i<n; ++i)
        {
            k2[i] *= dt;
            k_3[i] = x[i] + k2[i]/2;
        }

        derivatives(t + dt/2, k_3, k3);
        for(i=0; i<n; ++i)
        {
            k3[i] *= dt;
            k_4[i] = x[i] + k3[i];
        }

        derivatives(t + dt, k_4, k4);
        for(i=0; i<n; ++i)
        {
            k4[i] *= dt;
            x[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
        }

        t += dt;
    }

    t = to; // update current time, should be there but
    // what about round off.
}

#endif //#ifndef __RK4__HPP__
