#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <functional>
using namespace std;

struct HandlerData {
    double T, lim_start, lim_end;
};

struct CalculatableData {
    double t, j;
    int n_min, n_max;
    function<double(int)> func;
};

class FourierHandler {
    double T, lim_start, lim_end;
    public:
        FourierHandler(const HandlerData& handler_data) 
            : T(handler_data.T), lim_start(handler_data.lim_start), lim_end(handler_data.lim_end) {}
        void calculate(const CalculatableData& calculatable_data) {
            const auto t = calculatable_data.t;
            const auto j = calculatable_data.j;
            const auto n_min = calculatable_data.n_min;
            const auto n_max = calculatable_data.n_max;
            const auto func = calculatable_data.func;

            const auto n_count = n_max - n_min + 1;
            const auto lim_diff = lim_end - lim_start;

            const auto w = (2 *  M_PI) / T;
            const auto w0 = M_PI_2 / T;
            cout << "=== === w === ===" << endl;
            cout << "w = " << "(2 * n) / " << T << " = " << w << endl;
            cout << "w0 = " << "n / (2 * " << T << ") = " << w0 << endl;
            cout << "=== === = === ===" << endl;

            cout << "=== === x(t) === ===" << endl;
            vector<double> x_arr(n_count);
            for (auto i = 0; i < n_count; i++) {
                // (2n + 3) / sin(n + 2)
                const auto n = n_min + i;
                x_arr[i] = func(n);
                cout << "x(t)" << i + 1 << " = " << x_arr[i] << endl;
            }
            cout << "=== === ==== === ===\n" << endl;

            cout << "=== === a0 === ===" << endl;
            vector<double> a0_arr(n_count);
            for (auto i = 0; i < n_count; i++) {
                a0_arr[i] = 2 * x_arr[i] * lim_diff / T;
                cout << "a0" << i << " = " << a0_arr[i] << endl;
            }
            cout << "=== === == === ===\n" << endl;
            const auto abn_count = n_count * n_count;
            vector<double> an_arr(abn_count);
            vector<double> bn_arr(abn_count);
            for (auto i = 0; i < n_count; i++) {
                for (auto j = 0; j < n_count; j++) {
                    const auto index = i * n_count + j;
                    const auto n = n_min + j;

                    const auto func_arg = (2. * M_PI * t * n) / T;
                    const auto factor = (2. * x_arr[i] * lim_diff) / T;

                    an_arr[index] = factor * cos(func_arg);
                    bn_arr[index] = factor * sin(func_arg);
                }
            }

            cout << "=== === an === ===" << endl;
            for (auto i = 0; i < abn_count; i++) {
                cout << "an" << i + 1 << " = " << an_arr[i] << endl;
            }
            cout << "=== === == === ===\n" << endl;

            cout << "=== === bn === ===" << endl;
            for (auto i = 0; i < abn_count; i++) {
                cout << "bn" << i + 1 << " = " << bn_arr[i] << endl;
            }
            cout << "=== === == === ===\n" << endl;
            
            cout << "=== === === === === === === === X(t) === === === === === === === ===" << endl;
            const auto X_count = n_count * n_count * abn_count;
            vector<double> X_arr(X_count);

            const auto n_count_2 = n_count * n_count;
            const auto n_count_3 = n_count_2 * n_count;
            
            // i for a_0_i
            for (auto i = 0; i < n_count; i++) {
                const auto free = a0_arr[i] / 2;

                // j for n_j
                for (auto __j = 0; __j < n_count; __j++) {
                    const auto n = n_min + __j;
                    const auto func_arg = (2. * M_PI * t * n) / T;
                    
                    // k for a_n_k
                    for (auto k = 0; k < n_count; k++) {
                        const auto ancos = an_arr[k] * cos(func_arg);

                        // l for b_n_l
                        for (auto l = 0; l < n_count; l++) {
                            const auto bnsin = bn_arr[l] * sin(func_arg);
                            const auto index = i * n_count_3 + __j * n_count_2 + k * n_count + l;
                            X_arr[index] = free + ancos + bnsin;

                            cout << "X(t)" << index + 1 << " = ";
                            cout << an_arr[k] << " * cos(" << func_arg << ") + ";
                            cout << bn_arr[l] << " * sin(" << func_arg << ") = ";
                            cout << X_arr[index] << endl;
                        }
                    }
                }
            }
            cout << "=== === === === === === === === ==== === === === === === === === ===\n" << endl;

            cout << "=== === === === === === G === === === === === ===" << endl;
            vector<double> G_arr(X_count);
            for (auto i = 0; i < n_count; i++) {
                for (auto __j = 0; __j < n_count_3; __j++) {
                    const auto index = i * n_count_3 + __j;
                    G_arr[index] = sqrt(pow(x_arr[i], 2) + pow(X_arr[index], 2));

                    cout << "G" << index + 1 << " = sqrt(" << x_arr[i] << "^2 + " << X_arr[index] << "^2) = " << G_arr[index] << endl;
                }
            }
            cout << "=== === === === === === = === === === === === ===\n" << endl;

            cout << "=== === === === === === === === === Cw === === === === === === === === ===" << endl;
            const auto t2 = pow(t, 2);
            const auto jwt = j * w * t;

            const auto Cψ_1 = lim_diff * (pow(1 - t, 2) * pow(M_E, -t2 / 2)) / w;
            cout << "Cw_1 = int [(1 - " << t << ")^2 * e^(-" << t << "^2 / 2)] / " << w << " dw = " << Cψ_1 << endl;

            const auto Cψ_2 = lim_diff * (pow(M_E, -t2 / 2.) - pow(M_E, -t2 / 8) / 2) / w;
            cout << "Cw_2 = int [e^(" << -t << "^2 / 2) - 0.5 * e^(" << -t << "^2 / 8)] / " << w << " dw = " << Cψ_2 << endl;

            const auto Cψ_3 = lim_diff * (pow(M_E, jwt) - pow(M_E, -t2 / 2)) / w;
            cout << "Cw_3 = int [e^" << jwt << " - e^(" << -t << "^2 / 8)] / " << w << " dw = " << Cψ_3 << endl;
            cout << "=== === === === === === === === === == === === === === === === === === ===\n" << endl;

            
        }
};

int main() {
    const HandlerData handler_data_1 = { 3.8, 2.1, 4.3 };
    const HandlerData handler_data_2 = { 4.1, 1.7, 5.6 };
    const CalculatableData calculatable_data = { 2.4, -0.62, 6, 8, [](int n) -> double {
        return (2. * n + 3) / sin(n + 2.);
    }};

    auto ps1 = new FourierHandler(handler_data_1);
    ps1->calculate(calculatable_data);

    // auto ps2 = new FourierHandler(handler_data_2);
    // ps2->calculate(calculatable_data);

    return 0;
};