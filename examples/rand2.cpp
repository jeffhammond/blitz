#include <blitz/blitz.h>
#include <random/exponential.h>
#include <random/discrete-uniform.h>
#include <random/F.h>

using namespace blitz;
using namespace ranlib;

#include <time.h>

int main2()
{
    DiscreteUniform<int> rng(100);
    rng.seed((unsigned int)time(0));

    for (int i=0; i < 100; ++i)
        cout << rng.random() << "  ";
    cout << endl;
    return 0;
}

int main3()
{
    MersenneTwister x;

    for (int j=0; j<1000; j++) {
        printf("%10u ", x.random());
        if (j%8==7) printf("\n");
    }

    return 0;
}

int main()
{
    F<long double> rng(2.0,3.0);
    rng.seed((unsigned int)time(0)); 

    long double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;

    const int N = 10000;

    for (int i=0; i < N; ++i)
    {
        long double x = rng.random();
        sum1 += x;
        sum2 += x*x;
        sum3 += x*x*x;
        sum4 += x*x*x*x;
    }

    cout << "k1 = " << (sum1/N) << endl
         << "k2 = " << (sum2/N) << endl
         << "k3 = " << (sum3/N) << endl
         << "k4 = " << (sum4/N) << endl;

    return 0;
}

