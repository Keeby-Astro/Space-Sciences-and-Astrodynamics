#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
/*
    RSmooth - Robust smoothing of gridded data in one and higher dimensions.

    This function implements the algorithm described by:
    Garcia D. Robust smoothing of gridded data in one and higher dimensions with missing values.
    Comput Stat Data Anal. 2010 Apr 1;54(4):1167-1178. doi: 10.1016/j.csda.2009.09.020.
    PMID: 24795488; PMCID: PMC4008475.

    y - Input data array to be smoothed.
    z - Smoothed output array.

    Note: Please ensure the proper acknowledgment of the original author
          and source when using this code. This is a direct translation
          of the MATLAB code provided in the paper to C++.
*/

// Compute the median of a 1D vector
double median(std::vector<double> data) 
{
    size_t n = data.size();
    if (n == 0) return 0.0;
    std::sort(data.begin(), data.end());
    if (n % 2 == 1) 
    {
        return data[n / 2];
    } 
    else 
    {
        return 0.5 * (data[n / 2 - 1] + data[n / 2]);
    }
}

// 1D DCT (Type II) for a vector of length N
std::vector<double> dct1d(const std::vector<double> &in) 
{
    size_t N = in.size();
    if (N == 0) return {};

    std::vector<double> out(N, 0.0);
    const double factor = M_PI / (2.0 * N);

    for (size_t k = 0; k < N; ++k) 
    {
        double sumVal = 0.0;
        for (size_t n = 0; n < N; ++n) 
        {
            sumVal += in[n] * std::cos((2.0 * n + 1.0) * k * factor);
        }
        double alpha = (k == 0) ? std::sqrt(1.0 / N) : std::sqrt(2.0 / N);
        out[k] = alpha * sumVal;
    }
    return out;
}

// 1D IDCT (Inverse DCT-II) for a vector of length N
std::vector<double> idct1d(const std::vector<double> &in) 
{
    size_t N = in.size();
    if (N == 0) return {};

    std::vector<double> out(N, 0.0);
    const double factor = M_PI / (2.0 * N);

    for (size_t n = 0; n < N; ++n) 
    {
        double sumVal = 0.0;
        for (size_t k = 0; k < N; ++k) 
        {
            double alpha = (k == 0) ? std::sqrt(1.0 / N) : std::sqrt(2.0 / N);
            sumVal += alpha * in[k] * std::cos((2.0 * n + 1.0) * k * factor);
        }
        out[n] = sumVal;
    }
    return out;
}

// 2D DCT: apply 1D DCT row-wise, then column-wise
std::vector<std::vector<double>> dct2(const std::vector<std::vector<double>> &in) 
{
    size_t N1 = in.size();
    if (N1 == 0) return {};
    size_t N2 = in[0].size();
    if (N2 == 0) return {};

    // Copy input
    std::vector<std::vector<double>> temp = in;

    for (size_t i = 0; i < N1; ++i) 
    {
        temp[i] = dct1d(temp[i]);
    }

    std::vector<std::vector<double>> out(N1, std::vector<double>(N2, 0.0));
    for (size_t j = 0; j < N2; ++j) 
    {
        // Extract column
        std::vector<double> col(N1);
        for (size_t i = 0; i < N1; ++i) 
        {
            col[i] = temp[i][j];
        }
        col = dct1d(col);
        for (size_t i = 0; i < N1; ++i) 
        {
            out[i][j] = col[i];
        }
    }
    return out;
}

// 2D IDCT: Apply 1D IDCT column-wise, then row-wise
std::vector<std::vector<double>> idct2(const std::vector<std::vector<double>> &in) 
{
    size_t N1 = in.size();
    if (N1 == 0) return {};
    size_t N2 = in[0].size();
    if (N2 == 0) return {};

    // Copy input
    std::vector<std::vector<double>> temp = in;

    for (size_t j = 0; j < N2; ++j) 
    {
        // Extract column
        std::vector<double> col(N1);
        for (size_t i = 0; i < N1; ++i) 
        {
            col[i] = temp[i][j];
        }
        col = idct1d(col);
        for (size_t i = 0; i < N1; ++i) 
        {
            temp[i][j] = col[i];
        }
    }

    std::vector<std::vector<double>> out(N1, std::vector<double>(N2, 0.0));
    for (size_t i = 0; i < N1; ++i) 
    {
        out[i] = idct1d(temp[i]);
    }
    return out;
}

// Bisquare function for robust weighting
std::vector<std::vector<double>> bisquare(const std::vector<std::vector<double>> &r, double h) 
{
    // Flatten r to get median
    std::vector<double> flat;
    flat.reserve(r.size() * (r.empty()?0:r[0].size()));
    for (const auto &row : r) 
    {
        for (auto val : row) 
        {
            flat.push_back(val);
        }
    }

    double med_all = median(flat);
    // Compute absolute deviations from median
    for (auto &val : flat) 
    {
        val = std::fabs(val - med_all);
    }
    double MAD = median(flat);
    if (MAD < 1e-12) {
        MAD = 1e-12;
    }

    double denom = 1.4826 * MAD * std::sqrt(std::max(1e-12, 1.0 - h));
    std::vector<std::vector<double>> W(r.size(), std::vector<double>(r[0].size(), 0.0));

    for (size_t i = 0; i < r.size(); ++i) 
    {
        for (size_t j = 0; j < r[i].size(); ++j) 
        {
            double u = std::fabs(r[i][j]) / denom;
            double t = u / 4.685;
            if (t < 1.0) 
            {
                double tmp = (1.0 - t * t);
                W[i][j] = tmp * tmp; 
            } 
            else 
            {
                W[i][j] = 0.0;
            }
        }
    }
    return W;
}

// "fminbnd"-like search on p in [pMin, pMax]
struct GCVParams 
{
    const std::vector<std::vector<double>> &DCTy;
    const std::vector<std::vector<double>> &y;
    const std::vector<std::vector<double>> &W;
    const std::vector<std::vector<double>> &Lambda;
    double &bestS;       
    double &bestTrH;  
    std::vector<std::vector<double>> &bestZ;
    size_t n;
};

// Compute the GCV score for a given p
double GCVscore(double p, const GCVParams &params) 
{
    double s = std::pow(10.0, p);
    size_t N1 = params.Lambda.size();
    size_t N2 = (N1 == 0)? 0 : params.Lambda[0].size();

    std::vector<std::vector<double>> Gamma(N1, std::vector<double>(N2, 0.0));
    for (size_t i = 0; i < N1; ++i) 
    {
        for (size_t j = 0; j < N2; ++j) 
        {
            double val = params.Lambda[i][j];
            double denom = 1.0 + s * val * val;
            Gamma[i][j] = 1.0 / denom;
        }
    }

    std::vector<std::vector<double>> Gmul(N1, std::vector<double>(N2, 0.0));
    for (size_t i = 0; i < N1; ++i) 
    {
        for (size_t j = 0; j < N2; ++j) 
        {
            Gmul[i][j] = Gamma[i][j] * params.DCTy[i][j];
        }
    }
    std::vector<std::vector<double>> z = idct2(Gmul);

    double RSS = 0.0;
    for (size_t i = 0; i < N1; ++i) 
    {
        for (size_t j = 0; j < N2; ++j) 
        {
            double diff = params.y[i][j] - z[i][j];
            RSS += params.W[i][j] * diff * diff;
        }
    }

    double TrH = 0.0;
    for (size_t i = 0; i < N1; ++i) 
    {
        for (size_t j = 0; j < N2; ++j) 
        {
            TrH += Gamma[i][j];
        }
    }

    double n = static_cast<double>(params.n);
    double tmp = 1.0 - TrH / n;
    if (std::fabs(tmp) < 1e-15) {
        tmp = (tmp >= 0)? 1e-15 : -1e-15;
    }
    double GCVs = (RSS / n) / (tmp * tmp);
    return GCVs;
}

// Find the best p in [pMin, pMax] for GCV score
double findBestP(double pMin, double pMax, const GCVParams &params) 
{
    const int NUM_SAMPLES = 50;
    double bestP = pMin;
    double bestVal = std::numeric_limits<double>::infinity();
    double step = (pMax - pMin) / (NUM_SAMPLES - 1);

    for (int i = 0; i < NUM_SAMPLES; ++i) 
    {
        double p = pMin + i * step;
        double val = GCVscore(p, params);
        if (val < bestVal) 
        {
            bestVal = val;
            bestP = p;
        }
    }

    return bestP;
}

std::vector<std::vector<double>> rsmooth(const std::vector<std::vector<double>> &y) 
{
    // Check input size
    size_t n1 = y.size();
    if (n1 == 0) return {};
    size_t n2 = y[0].size();
    if (n2 == 0) return {};

    size_t n = n1 * n2;
    int N = 0;
    if (n1 > 1) ++N;
    if (n2 > 1) ++N;

    // Build the Lambda matrix
    std::vector<double> rowLambda(n2), colLambda(n1);
    for (size_t j = 0; j < n2; ++j) 
    {
        rowLambda[j] = -2.0 + 2.0 * std::cos((double)j * M_PI / (double)n2);
    }
    for (size_t i = 0; i < n1; ++i) 
    {
        colLambda[i] = -2.0 + 2.0 * std::cos((double)i * M_PI / (double)n1);
    }
    std::vector<std::vector<double>> Lambda(n1, std::vector<double>(n2, 0.0));
    for (size_t i = 0; i < n1; ++i) 
    {
        for (size_t j = 0; j < n2; ++j) 
        {
            Lambda[i][j] = rowLambda[j] + colLambda[i];
        }
    }

    std::vector<std::vector<double>> W(n1, std::vector<double>(n2, 1.0));
    std::vector<std::vector<double>> z = y;   // Final estimate
    std::vector<std::vector<double>> zz = y;  // Temp old estimate

    double sVal = 0.0; 
    for (int iter = 0; iter < 6; ++iter) 
    {
        double tol = std::numeric_limits<double>::infinity();

        while (tol > 1e-5) 
        {
            std::vector<std::vector<double>> toDCT(n1, std::vector<double>(n2, 0.0));
            for (size_t i = 0; i < n1; ++i) {
                for (size_t j = 0; j < n2; ++j) 
                {
                    toDCT[i][j] = W[i][j] * (y[i][j] - zz[i][j]) + zz[i][j];
                }
            }
            auto DCTy = dct2(toDCT);
            double bestTrH = 0.0;          
            std::vector<std::vector<double>> bestZ;
            
            GCVParams params{DCTy, y, W, Lambda, sVal, bestTrH, bestZ, n};
            double pMin = -15.0;
            double pMax = 38.0;
            double bestP = findBestP(pMin, pMax, params);

            double s = std::pow(10.0, bestP);
            std::vector<std::vector<double>> Gamma(n1, std::vector<double>(n2, 0.0));
            for (size_t i = 0; i < n1; ++i) 
            {
                for (size_t j = 0; j < n2; ++j) 
                {
                    double val = Lambda[i][j];
                    double denom = 1.0 + s * val * val;
                    Gamma[i][j] = 1.0 / denom;
                }
            }

            std::vector<std::vector<double>> Gmul(n1, std::vector<double>(n2, 0.0));
            for (size_t i = 0; i < n1; ++i)
             {
                for (size_t j = 0; j < n2; ++j) 
                {
                    Gmul[i][j] = Gamma[i][j] * DCTy[i][j];
                }
            }
            auto newZ = idct2(Gmul);

            // Check tolerance
            double normDiff = 0.0;
            double normNewZ = 0.0;
            for (size_t i = 0; i < n1; ++i) 
            {
                for (size_t j = 0; j < n2; ++j) 
                {
                    double d = zz[i][j] - newZ[i][j];
                    normDiff += d * d;
                    normNewZ += newZ[i][j] * newZ[i][j];
                }
            }
            normDiff = std::sqrt(normDiff);
            normNewZ = (normNewZ < 1e-15)? 1e-15 : std::sqrt(normNewZ);
            tol = normDiff / normNewZ;

            // Update z
            z = newZ;
            zz = newZ;
            sVal = s;
        }

        double tmp = std::sqrt(1.0 + 16.0 * sVal);
        double h = (std::sqrt(1.0 + tmp) / std::sqrt(2.0) / tmp);
        h = std::pow(h, (double)N);

        std::vector<std::vector<double>> r(n1, std::vector<double>(n2, 0.0));
        for (size_t i = 0; i < n1; ++i) 
        {
            for (size_t j = 0; j < n2; ++j) 
            {
                r[i][j] = y[i][j] - z[i][j];
            }
        }
        W = bisquare(r, h);
    }
    return z;
}