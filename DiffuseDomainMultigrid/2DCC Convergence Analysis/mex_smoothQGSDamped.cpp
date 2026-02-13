/*
 * mex_smoothQGSDamped.cpp
 *
 * Modern C++ MEX implementation (Fixed for compilation errors).
 */

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstring>
#include <vector>
#include <algorithm> // Required for std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    // Helper: Column-major index for U (which has ghost layers)
    inline size_t IDU(size_t i, size_t j, size_t rows) {
        return i + j * rows;
    }

    // Apply Homogeneous Neumann BCs
    void applyBCs(double* u, size_t nx, size_t ny) {
        size_t rows = nx + 2;
        size_t cols = ny + 2;

        // Left (j=0) and Right (j=ny+1) Ghost Columns
        for (size_t i = 1; i <= nx; ++i) {
            u[IDU(i, 0, rows)]      = u[IDU(i, 1, rows)];
            u[IDU(i, ny + 1, rows)] = u[IDU(i, ny, rows)];
        }

        // Bottom (i=0) and Top (i=nx+1) Ghost Rows
        for (size_t j = 0; j < cols; ++j) {
            u[IDU(0, j, rows)]      = u[IDU(1, j, rows)];
            u[IDU(nx + 1, j, rows)] = u[IDU(nx, j, rows)];
        }
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // 1. Initialize Factory first (fix for 'undeclared identifier')
        ArrayFactory factory;

        // Check inputs
        if (inputs.size() != 9) {
            // 2. Use '->' pointer access (fix for 'no member named feval')
            getEngine()->feval(u"error", 
                0, std::vector<Array>({ factory.createScalar("9 inputs required.") }));
        }

        // --- Unpack Inputs ---
        TypedArray<double> fArr   = inputs[0];
        TypedArray<double> DeWArr = inputs[1];
        TypedArray<double> DnSArr = inputs[2];
        TypedArray<double> CArr   = inputs[3];
        TypedArray<double> uInArr = inputs[4];

        // Scalars
        double hf    = inputs[5][0];
        int m        = (int)inputs[6][0];
        double omega = inputs[7][0];

        // String
        CharArray dirChar = inputs[8];
        std::string direction = dirChar.toAscii();
        bool isFwd = (direction == "fwd");

        // --- Dimensions ---
        auto dims = fArr.getDimensions();
        size_t nx = dims[0];
        size_t ny = dims[1];
        size_t uRows = nx + 2;

        // --- Create Output U ---
        // 3. Fix for 'createArray' ambiguity: Create uninitialized, then copy.
        TypedArray<double> uOut = factory.createArray<double>(uInArr.getDimensions());
        std::copy(uInArr.begin(), uInArr.end(), uOut.begin());

        // --- Get Raw Pointers ---
        double* f   = &(*fArr.begin());
        double* DeW = &(*DeWArr.begin());
        double* DnS = &(*DnSArr.begin());
        double* C   = &(*CArr.begin());
        double* u   = &(*uOut.begin());

        // Precompute constants
        double hf2 = hf * hf;
        double omegaPrime = 1.0 - omega;

        // Apply Initial BCs
        applyBCs(u, nx, ny);

        // --- Smoothing Loop ---
        for (int sweep = 0; sweep < m; ++sweep) {
            if (isFwd) {
                // Forward Sweep
                for (size_t j = 0; j < ny; ++j) {
                    for (size_t i = 0; i < nx; ++i) {
                        size_t idx_int = i + j * nx; 

                        // Indices in U
                        size_t idx_u_C = IDU(i+1, j+1, uRows);
                        size_t idx_u_E = IDU(i+2, j+1, uRows);
                        size_t idx_u_W = IDU(i,   j+1, uRows);
                        size_t idx_u_N = IDU(i+1, j+2, uRows);
                        size_t idx_u_S = IDU(i+1, j,   uRows);

                        // Coefficients
                        size_t rowsDeW = nx + 1;
                        double D_E = DeW[(i+1) + j * rowsDeW];
                        double D_W = DeW[i     + j * rowsDeW];

                        double D_N = DnS[i + (j+1) * nx];
                        double D_S = DnS[i + j      * nx];

                        // Calculation
                        double num = hf2 * f[idx_int]
                                   + u[idx_u_E] * D_E
                                   + u[idx_u_W] * D_W
                                   + u[idx_u_N] * D_N
                                   + u[idx_u_S] * D_S;

                        double den = D_E + D_W + D_N + D_S + C[idx_int] * hf2;

                        // Gauss-Seidel Update
                        u[idx_u_C] = omega * (num / den) + omegaPrime * u[idx_u_C];
                    }
                }
            } else {
                // Backward Sweep
                for (long long j = ny - 1; j >= 0; --j) {
                    for (long long i = nx - 1; i >= 0; --i) {
                        size_t idx_int = i + j * nx;

                        size_t idx_u_C = IDU(i+1, j+1, uRows);
                        size_t idx_u_E = IDU(i+2, j+1, uRows);
                        size_t idx_u_W = IDU(i,   j+1, uRows);
                        size_t idx_u_N = IDU(i+1, j+2, uRows);
                        size_t idx_u_S = IDU(i+1, j,   uRows);

                        size_t rowsDeW = nx + 1;
                        double D_E = DeW[(i+1) + j * rowsDeW];
                        double D_W = DeW[i     + j * rowsDeW];

                        double D_N = DnS[i + (j+1) * nx];
                        double D_S = DnS[i + j      * nx];

                        double num = hf2 * f[idx_int]
                                   + u[idx_u_E] * D_E
                                   + u[idx_u_W] * D_W
                                   + u[idx_u_N] * D_N
                                   + u[idx_u_S] * D_S;

                        double den = D_E + D_W + D_N + D_S + C[idx_int] * hf2;

                        u[idx_u_C] = omega * (num / den) + omegaPrime * u[idx_u_C];
                    }
                }
            }
            applyBCs(u, nx, ny);
        }

        outputs[0] = uOut;
    }
};