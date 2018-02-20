#include <stdio.h>
#include <math.h>

#include "../include/rk4.hpp"


class Sin : public RK4<float, float, float>
{
    public:

        Sin(float period);

        void setPeriod(float period_in)
        {
            period = period_in;
            angularFreq_2 = 2.0F*M_PI/period;
            angularFreq_2 *= angularFreq_2;
            tStep = period/6.0F;
        };

        float getPeriod(void) { return period; };

    protected:

        void derivatives(float t, const float *x, float *xDot);

    private:

        float angularFreq_2, period;
};

Sin::Sin(float period):
    RK4(2/*num degree of freedom*/,
            0.0F/*time step gets set in setPeriod()*/)
{
    setPeriod(period);
}

void Sin::derivatives(float t, const float *x, float *xDot)
{
    xDot[0] = x[1]; // x_dot = v
    xDot[1] = - angularFreq_2 * x[0];
}


int main(void)
{
    float x[2] = { 1.0, 0.0 };
    float t, fromT;
    Sin sin(1.0F);

    for(t=0.0F; t<10.0F;)
    {
        fromT = t;
        sin.go(x, fromT, t += 0.01);
        printf("%g %g %g %g\n", t, x[0], x[1], cosf(2.0F*M_PI*t));
#if 1
        if(sin.getPeriod() > 5.0F)
            sin.setPeriod(1.0F);
        else
            sin.setPeriod(sin.getPeriod() * 1.003);
#endif
    }

    return 0;
}
