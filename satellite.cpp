#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

constexpr double G = 6.67e-11;
constexpr double M_EARTH = 6e24;
constexpr double L_PARAM = 1.0;
constexpr double TIME_STEP = 0.01;
constexpr double MIN_MASS = 5000.0;
constexpr size_t RESERVE_CAPACITY = 1000000;

struct SystemState {
    double t, r, v, m, E;
    bool engine_active;
};

enum class SimMode {
    CONSTANT_MASS_RK4 = 1,
    VARIABLE_MASS_RK4 = 2,
    VARIABLE_MASS_EULER = 3
};

std::vector<SystemState> run_simulation(SimMode mode) {
    std::vector<SystemState> history;
    // Pre-allocation mandatory to prevent heap fragmentation in the fast-path
    history.reserve(RESERVE_CAPACITY);

    SystemState state = { 0.0, 4e5, std::sqrt(G * M_EARTH / 4e5), 6000.0, 1.0, false };
    double E0 = G * M_EARTH / state.r + 0.5 * state.m * state.v * state.v;
    
    if (mode != SimMode::CONSTANT_MASS_RK4) {
        state.E = E0;
    }

    const double energy_cutoff = (mode == SimMode::CONSTANT_MASS_RK4) ? 1e-9 : 0.01;

    while (state.E > energy_cutoff) {
        state.engine_active = (mode != SimMode::CONSTANT_MASS_RK4 && state.m > MIN_MASS);
        history.push_back(state);

        // Local aliasing preserves legacy RK4 formula syntax without abstraction overhead
        double& t = state.t;
        double& r = state.r;
        double& v = state.v;
        double& m = state.m;
        double& E = state.E;

        double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, dv, dr, dm;

        if (mode == SimMode::CONSTANT_MASS_RK4) {
            a1 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * r / (m * m)) * TIME_STEP;
            b1 = -2 * ((L_PARAM * L_PARAM * r) / m) * TIME_STEP;

            a2 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b1) / (m * m)) * TIME_STEP;
            b2 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b1)) / m) * TIME_STEP;

            a3 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b2) / (m * m)) * TIME_STEP;
            b3 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b2)) / m) * TIME_STEP;

            a4 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + b3) / (m * m)) * TIME_STEP;
            b4 = -2 * ((L_PARAM * L_PARAM * (r + b3)) / m) * TIME_STEP;

            dv = (a1 + 2 * a2 + 2 * a3 + a4) / 6.0;
            dr = (b1 + 2 * b2 + 2 * b3 + b4) / 6.0;

            v += dv;
            r += dr;
            t += TIME_STEP;
            E += -(L_PARAM * L_PARAM * G * M_EARTH / (m * r)) / E0;

        } else if (mode == SimMode::VARIABLE_MASS_RK4) {
            if (state.engine_active) {
                c1 = -(TIME_STEP * G * M_EARTH * m / (v * r * r)) - (L_PARAM * L_PARAM * TIME_STEP);
                a1 = -TIME_STEP * G * M_EARTH / (r * r) - (v / m) * (L_PARAM * L_PARAM * TIME_STEP + c1);
                b1 = -2 * ((L_PARAM * L_PARAM * r) / m) * TIME_STEP - 2 * (r / m) * c1;

                c2 = -(TIME_STEP * G * M_EARTH * (m + 0.5 * c1) / ((v + 0.5 * a1) * (r + 0.5 * b1) * (r + 0.5 * b1))) - (L_PARAM * L_PARAM * TIME_STEP);
                a2 = -TIME_STEP * G * M_EARTH / ((r + 0.5 * b1) * (r + 0.5 * b1)) - ((v + 0.5 * a1) / (m + 0.5 * c1)) * (L_PARAM * L_PARAM * TIME_STEP + c2);
                b2 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b1)) / (m + 0.5 * c1)) * TIME_STEP - 2 * ((r + 0.5 * b1) / (m + 0.5 * c1)) * c2;

                c3 = -(TIME_STEP * G * M_EARTH * (m + 0.5 * c2) / ((v + 0.5 * a2) * (r + 0.5 * b2) * (r + 0.5 * b2))) - (L_PARAM * L_PARAM * TIME_STEP);
                a3 = -TIME_STEP * G * M_EARTH / ((r + 0.5 * b2) * (r + 0.5 * b2)) - ((v + 0.5 * a2) / (m + 0.5 * c2)) * (L_PARAM * L_PARAM * TIME_STEP + c3);
                b3 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b2)) / (m + 0.5 * c2)) * TIME_STEP - 2 * ((r + 0.5 * b2) / (m + 0.5 * c2)) * c3;

                c4 = -(TIME_STEP * G * M_EARTH * (m + c3) / ((v + a3) * (r + b3) * (r + b3))) - (L_PARAM * L_PARAM * TIME_STEP);
                a4 = -TIME_STEP * G * M_EARTH / ((r + b3) * (r + b3)) - ((v + a3) / (m + c3)) * (L_PARAM * L_PARAM * TIME_STEP + c4);
                b4 = -2 * ((L_PARAM * L_PARAM * (r + b3)) / (m + c3)) * TIME_STEP - 2 * ((r + b3) / (m + c3)) * c4;

                dv = (a1 + 2 * a2 + 2 * a3 + a4) / 6.0;
                dr = (b1 + 2 * b2 + 2 * b3 + b4) / 6.0;
                dm = (c1 + 2 * c2 + 2 * c3 + c4) / 6.0;

                v += dv;
                r += dr;
                m += dm;
                t += TIME_STEP;
                E += -v * v * (L_PARAM * L_PARAM * TIME_STEP + dm);

            } else {
                a1 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * r / (m * m)) * TIME_STEP;
                b1 = -2 * ((L_PARAM * L_PARAM * r) / m) * TIME_STEP;

                a2 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b1) / (m * m)) * TIME_STEP;
                b2 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b1)) / m) * TIME_STEP;

                a3 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b2) / (m * m)) * TIME_STEP;
                b3 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b2)) / m) * TIME_STEP;

                a4 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + b3) / (m * m)) * TIME_STEP;
                b4 = -2 * ((L_PARAM * L_PARAM * (r + b3)) / m) * TIME_STEP;

                dv = (a1 + 2 * a2 + 2 * a3 + a4) / 6.0;
                dr = (b1 + 2 * b2 + 2 * b3 + b4) / 6.0;

                v += dv;
                r += dr;
                t += TIME_STEP;
                E += -(L_PARAM * L_PARAM * G * M_EARTH / (m * r));
            }
        } else if (mode == SimMode::VARIABLE_MASS_EULER) {
            if (state.engine_active) {
                m += -L_PARAM * L_PARAM * TIME_STEP;
                v -= -G * M_EARTH * TIME_STEP / (r * r * 10000.0);
                t += TIME_STEP;
            } else {
                a1 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * r / (m * m)) * TIME_STEP;
                b1 = -2 * ((L_PARAM * L_PARAM * r) / m) * TIME_STEP;

                a2 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b1) / (m * m)) * TIME_STEP;
                b2 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b1)) / m) * TIME_STEP;

                a3 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + 0.5 * b2) / (m * m)) * TIME_STEP;
                b3 = -2 * ((L_PARAM * L_PARAM * (r + 0.5 * b2)) / m) * TIME_STEP;

                a4 = 4 * (L_PARAM * L_PARAM * L_PARAM * L_PARAM * (r + b3) / (m * m)) * TIME_STEP;
                b4 = -2 * ((L_PARAM * L_PARAM * (r + b3)) / m) * TIME_STEP;

                dv = (a1 + 2 * a2 + 2 * a3 + a4) / 6.0;
                dr = (b1 + 2 * b2 + 2 * b3 + b4) / 6.0;

                v += dv;
                r += dr;
                t += TIME_STEP;
                E += -(L_PARAM * L_PARAM * G * M_EARTH / (m * r));
            }
        }
    }
    
    history.push_back(state);
    return history;
}

void save_data(const std::vector<SystemState>& data, SimMode mode) {
    std::ofstream out_file("sat.dat");
    if (!out_file) {
        std::cerr << "Fatal: Failed to acquire file lock for sat.dat\n";
        std::exit(EXIT_FAILURE);
    }

    // Output format strictly matches legacy behavior for downstream pipeline compatibility
    for (const auto& s : data) {
        if (mode == SimMode::CONSTANT_MASS_RK4) {
            out_file << s.t << "\t" << s.r << "\t" << s.v << "\t" << s.E << "\n";
        } else {
            out_file << s.t << "\t" << s.r << "\t" << s.v << "\t" << s.m << "\t" << s.E;
            if (s.engine_active) {
                out_file << "\tZ";
            }
            out_file << "\n";
        }
    }
}

int main() {
    std::cout << "1\t2\t3\n";
    
    int user_input;
    if (!(std::cin >> user_input) || user_input < 1 || user_input > 3) {
        std::cerr << "Invalid parameter. Terminating.\n";
        return EXIT_FAILURE;
    }

    SimMode mode = static_cast<SimMode>(user_input);
    
    std::vector<SystemState> results = run_simulation(mode);
    save_data(results, mode);
    
    const auto& last = results.back();
    std::cout << "\n" << last.t << "\t" << last.r << "\t" << last.v << "\t" << last.E << "\n\n";

    return 0;
}
