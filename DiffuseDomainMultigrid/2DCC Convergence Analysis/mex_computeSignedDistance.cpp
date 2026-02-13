/**
 * mex_computeSignedDistance.cpp
 * 
 * MODERN C++ MEX API VERSION
 * Fixes: Initialization Order Bug (Split Status update and Neighbor push)
 * 
 * Compile with: mex mex_computeSignedDistance.cpp
 */

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <queue>
#include <algorithm>
#include <string>

using namespace matlab::data;
using matlab::mex::ArgumentList;

const double INF = std::numeric_limits<double>::infinity();

struct Node {
    double val;
    int i, j;
    bool operator>(const Node& other) const { return val > other.val; }
};

inline size_t IDX(int i, int j, size_t rows) { return (size_t)i + (size_t)j * rows; }

double pointSegmentDistSq(double px, double py, double x1, double y1, double x2, double y2) {
    double l2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
    if (l2 == 0) return (px-x1)*(px-x1) + (py-y1)*(py-y1);
    double t = ((px-x1)*(x2-x1) + (py-y1)*(y2-y1)) / l2;
    t = std::max(0.0, std::min(1.0, t));
    double projx = x1 + t*(x2-x1);
    double projy = y1 + t*(y2-y1);
    return (px-projx)*(px-projx) + (py-projy)*(py-projy);
}

bool isInside(double px, double py, const double* polyX, const double* polyY, size_t nPoly) {
    int wn = 0;
    for (size_t i = 0; i < nPoly - 1; i++) {
        double y1 = polyY[i];
        double y2 = polyY[i+1];
        if (y1 <= py) {
            if (y2 > py) {
                if ((polyX[i+1] - polyX[i]) * (py - y1) - (px - polyX[i]) * (y2 - y1) > 0) wn++;
            }
        } else {
            if (y2 <= py) {
                if ((polyX[i+1] - polyX[i]) * (py - y1) - (px - polyX[i]) * (y2 - y1) < 0) wn--;
            }
        }
    }
    return wn != 0;
}

double solveUpwind(int i, int j, const std::vector<double>& phi, const std::vector<int>& status, double h, int nx, int ny) {
    double a = INF;
    double b = INF;
    if (i > 0 && status[IDX(i-1, j, nx)] == 2) a = std::min(a, phi[IDX(i-1, j, nx)]);
    if (i < nx-1 && status[IDX(i+1, j, nx)] == 2) a = std::min(a, std::min(a, phi[IDX(i+1, j, nx)]));
    if (j > 0 && status[IDX(i, j-1, nx)] == 2) b = std::min(b, phi[IDX(i, j-1, nx)]);
    if (j < ny-1 && status[IDX(i, j+1, nx)] == 2) b = std::min(b, std::min(b, phi[IDX(i, j+1, nx)]));
    
    if (std::isinf(a) && std::isinf(b)) return INF;
    if (std::isinf(a)) return b + h;
    if (std::isinf(b)) return a + h;
    double d = std::abs(a - b);
    if (d >= h) return std::min(a, b) + h;
    return 0.5 * (a + b + std::sqrt(2*h*h - d*d));
}

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        if (inputs.size() < 5) {
            getEngine()->feval(u"error", 0, { factory.createScalar("5 inputs required") });
            return;
        }

        TypedArray<double> arrX = std::move(inputs[0]);
        TypedArray<double> arrY = std::move(inputs[1]);
        TypedArray<double> arrPX = std::move(inputs[2]);
        TypedArray<double> arrPY = std::move(inputs[3]);
        TypedArray<double> arrH = std::move(inputs[4]);

        const double* X = &(*arrX.begin());
        const double* Y = &(*arrY.begin());
        const double* polyX = &(*arrPX.begin());
        const double* polyY = &(*arrPY.begin());
        double h = *arrH.begin(); 

        size_t nx = arrX.getDimensions()[0];
        size_t ny = arrX.getDimensions()[1];
        size_t nPoly = arrPX.getNumberOfElements();

        std::vector<double> phi(nx * ny, INF);
        std::vector<int> status(nx * ny, 0);
        std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

        double bandWidthSq = (2.5 * h) * (2.5 * h);

        // 1. Compute Exact Distance to Segments
        for (size_t k = 0; k < nPoly - 1; k++) {
            double x1 = polyX[k];     double y1 = polyY[k];
            double x2 = polyX[k+1];   double y2 = polyY[k+1];
            double minX = std::min(x1, x2) - 2.5*h;
            double maxX = std::max(x1, x2) + 2.5*h;
            double minY = std::min(y1, y2) - 2.5*h;
            double maxY = std::max(y1, y2) + 2.5*h;
            double xStart = X[0]; double yStart = Y[0];

            int i_start = std::max(0, (int)std::floor((minX - xStart)/h));
            int i_end   = std::min((int)nx-1, (int)std::ceil((maxX - xStart)/h));
            int j_start = std::max(0, (int)std::floor((minY - yStart)/h));
            int j_end   = std::min((int)ny-1, (int)std::ceil((maxY - yStart)/h));

            for (int j = j_start; j <= j_end; j++) {
                for (int i = i_start; i <= i_end; i++) {
                    int idx = IDX(i, j, nx);
                    double d2 = pointSegmentDistSq(X[idx], Y[idx], x1, y1, x2, y2);
                    if (d2 < bandWidthSq) {
                        double d = std::sqrt(d2);
                        if (d < phi[idx]) phi[idx] = d;
                    }
                }
            }
        }

        // 2. Pass 1: Mark all initialized points as KNOWN
        // We must do this BEFORE computing neighbors to ensure solveUpwind sees all boundary data.
        for (size_t idx = 0; idx < nx*ny; idx++) {
            if (phi[idx] != INF) {
                status[idx] = 2; // Known
            }
        }

        // 3. Pass 2: Initialize Heap with Neighbors of KNOWN points
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = IDX(i, j, nx);
                if (status[idx] == 2) { // If Known
                    int ni[4] = {i-1, i+1, i, i};
                    int nj[4] = {j, j, j-1, j+1};
                    for (int k=0; k<4; k++) {
                        int ii = ni[k]; int jj = nj[k];
                        if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
                            int nidx = IDX(ii, jj, nx);
                            if (status[nidx] == 0) {
                                status[nidx] = 1; // Trial
                                phi[nidx] = solveUpwind(ii, jj, phi, status, h, nx, ny);
                                pq.push({phi[nidx], ii, jj});
                            }
                        }
                    }
                }
            }
        }

        // 4. Fast Marching
        while (!pq.empty()) {
            Node top = pq.top();
            pq.pop();
            int i = top.i; int j = top.j;
            int idx = IDX(i, j, nx);

            if (status[idx] == 2) continue;
            status[idx] = 2; 

            int ni[4] = {i-1, i+1, i, i};
            int nj[4] = {j, j, j-1, j+1};

            for (int k=0; k<4; k++) {
                int ii = ni[k]; int jj = nj[k];
                if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
                    int nidx = IDX(ii, jj, nx);
                    if (status[nidx] != 2) {
                        double newVal = solveUpwind(ii, jj, phi, status, h, nx, ny);
                        if (newVal < phi[nidx]) {
                            phi[nidx] = newVal;
                            status[nidx] = 1; 
                            pq.push({newVal, ii, jj});
                        }
                    }
                }
            }
        }

        // 5. Output
        TypedArray<double> outPhi = factory.createArray<double>({nx, ny});
        double* outPtr = &(*outPhi.begin());

        for (size_t idx = 0; idx < nx*ny; idx++) {
            double val = phi[idx];
            if (!isInside(X[idx], Y[idx], polyX, polyY, nPoly)) val = -val;
            outPtr[idx] = val;
        }

        outputs[0] = outPhi;
    }
    ArrayFactory factory;
};