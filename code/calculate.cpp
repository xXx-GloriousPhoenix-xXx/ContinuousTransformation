#include <cmath>
#include <iostream>

int main() {
    const double T1 = 3.8,
                 T2 = 4.1,
                 t  = 2.4;

    const int n_min = 6,
              n_max = 8;

    const double lim_1_min = 2.1,
                 lim_1_max = 4.3,
                 lim_2_min = 1.7,
                 lim_2_max = 5.6;

    //

    const int n_count = n_max - n_min + 1;
    auto x = new int[n_count];
    for (auto i = 0; i < n_count; i++) {
        // (2n + 3) / sin(n + 2)
        auto n = n_min + i;
        x[i] = (2 * n + 3) / Math
    }



    return 0;
}