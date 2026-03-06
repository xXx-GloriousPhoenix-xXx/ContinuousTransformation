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

            cout << "=== === Run Data === ===" << endl;
            cout << "T = " << T << endl;
            cout << "lim = {" << lim_start << " ... " << lim_end << "}" << endl;
            cout << "=== === ======== === ===" << endl << endl;

            const auto n_count = n_max - n_min + 1;
            const auto lim_diff = lim_end - lim_start;

            const auto w = (2 *  M_PI) / T;
            const auto w0 = M_PI_2 / T;
            cout << "=== === w === ===" << endl;
            cout << "w = " << "(2 * n) / " << T << " = " << w << endl;
            cout << "w0 = " << "n / (2 * " << T << ") = " << w0 << endl;
            cout << "=== === = === ===" << endl << endl;

            cout << "=== === x(t) === ===" << endl;
            vector<double> x_arr(n_count);
            for (auto i = 0; i < n_count; i++) {
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

            const auto Cψ_1 = lim_diff * (pow(1 - t, 2) * exp(-t2 / 2)) / w;
            cout << "Cw_1 = int [(1 - " << t << ")^2 * e^(-" << t << "^2 / 2)] / " << w << " dw = " << Cψ_1 << endl;

            const auto Cψ_2 = lim_diff * (exp(-t2 / 2.) - exp(-t2 / 8) / 2) / w;
            cout << "Cw_2 = int [e^(" << -t << "^2 / 2) - 0.5 * e^(" << -t << "^2 / 8)] / " << w << " dw = " << Cψ_2 << endl;

            const auto Cψ_3 = lim_diff * (exp(jwt) - exp(-t2 / 2)) / w;
            cout << "Cw_3 = int [e^" << jwt << " - e^(" << -t << "^2 / 8)] / " << w << " dw = " << Cψ_3 << endl;
            cout << "=== === === === === === === === === == === === === === === === === === ===\n" << endl;

            cout << "=== === === === === === Wx(a,b) === === === === === ===" << endl;
            const auto wxab_factor = X_arr[0] * M_SQRT2 * M_PI;
            const auto jwe2 = pow(j * w, 2) * exp(-pow(w, 2) / 2);
            const auto ew2ew3 = exp(-pow(w, 2) / 2) - exp(-pow(w, 2) / 3);
            const auto weww02 = w * exp(-pow(w - w0, 2) / 2);
            const auto wxabweww02 = wxab_factor * weww02;

            cout << "X1 * sqrt(2) * n = " << X_arr[0] << " * " << M_SQRT2 << " * " << M_PI << " = " << wxab_factor << endl;
            cout << "(" << j << " * " << w << ")^2 * e^(-" << w << "^2 / 2) = " << jwe2 << endl;
            cout << "e^(-" << w << "^2 / 2) - e^(-" << w << "^2 / 3) = " << ew2ew3 << endl;
            cout << "w * e^(-(" << w << " - " << w0 << ")^2 / 2) = " << weww02 << endl;
            cout << wxabweww02 << " * " << weww02 << " = " << wxabweww02 << endl << endl;

            const auto wxab1 = wxab_factor * jwe2 * lim_diff;
            cout << "Wx(a,b)1 = int " << wxab_factor << " * " << jwe2 << " dt = " << wxab1 << endl; 

            const auto wxab2 = wxab_factor * ew2ew3 * lim_diff;
            cout << "Wx(a,b)2 = int " << wxab_factor << " * " << ew2ew3 << " dt = " << wxab2 << endl; 

            vector<double> wxab(X_count);
            for (auto i = 0; i < X_count; i++) {
                wxab[i] = wxab_factor * G_arr[i] * weww02 * lim_diff;
                cout << "Wx" << i + 1 << "(a,b) = int " << wxabweww02 << " * " << G_arr[i] << " dt = " << wxab[i] << endl;
            }
            cout << "=== === === === === === ======= === === === === === ===" << endl << endl;

            cout << "=== === === === === === X(t,a,b) === === === === === ===" << endl;
            const auto s2pjw2ew2 = sqrt(2) * M_PI * pow(j * w, 2) * exp(-pow(w, 2) / 2);
            const auto s2pjw2ew2Cw1 = s2pjw2ew2 / Cψ_1;
            const auto s2pjw2ew2Cw2 = s2pjw2ew2 / Cψ_2;
            const auto s2pjw2ew2Cw3 = s2pjw2ew2 / Cψ_3;

            cout << "sqrt(2) * n * (" << j << " * " << w << ")^2 * e^(-" << w << "^2 / 2) = " << s2pjw2ew2 << endl;
            cout << s2pjw2ew2 << " / " << Cψ_1 << " = " << s2pjw2ew2Cw1 << endl;
            cout << s2pjw2ew2 << " / " << Cψ_2 << " = " << s2pjw2ew2Cw2 << endl;
            cout << s2pjw2ew2 << " / " << Cψ_3 << " = " << s2pjw2ew2Cw3 << endl << endl;

            const auto Xtab1 = s2pjw2ew2Cw1 * wxab1 * pow(lim_diff, 2) / 2;
            cout << "X(t,a,b)1 = int int " << s2pjw2ew2Cw1 << " * " << wxab1 << " dadb = " << Xtab1 << endl;

            const auto Xtab2 = s2pjw2ew2Cw2 * wxab2 * pow(lim_diff, 2) / 2;
            cout << "X(t,a,b)2 = int int " << s2pjw2ew2Cw2 << " * " << wxab2 << " dadb = " << Xtab2 << endl;

            vector<double> Xtab(X_count);
            for (auto i = 0; i < X_count; i++) {
                Xtab[i] = s2pjw2ew2Cw3 * wxab[i] * pow(lim_diff, 2) / 2;
                cout << "X" << i + 1 << "(t,a,b) = int int " << s2pjw2ew2Cw3 << " * " << wxab[i] << " dadb = " << Xtab[i] << endl;
            }
            cout << "=== === === === === === ======== === === === === === ===" << endl;
        }
};

int main() {
    const HandlerData handler_data_1 = { 3.8, 2.1, 4.3 };
    const HandlerData handler_data_2 = { 4.1, 1.7, 5.6 };
    const CalculatableData calculatable_data = { 2.4, -0.62, 6, 8, [](int n) -> double {
        // return (2. * n + 3) / sin(n + 2.);
        return cos(2 * n + pow(n, 2));
    }};

    cout << "--- --- --- --- --- --- --- --- --- Run 1 --- --- --- --- --- --- --- --- ---" << endl;
    auto ps1 = new FourierHandler(handler_data_1);
    ps1->calculate(calculatable_data);
    cout << "--- --- --- --- --- --- --- --- --- ----- --- --- --- --- --- --- --- --- ---" << endl;

    cout << "--- --- --- --- --- --- --- --- --- Run 2 --- --- --- --- --- --- --- --- ---" << endl;
    auto ps2 = new FourierHandler(handler_data_2);
    ps2->calculate(calculatable_data);
    cout << "--- --- --- --- --- --- --- --- --- ----- --- --- --- --- --- --- --- --- ---" << endl;

    return 0;
};