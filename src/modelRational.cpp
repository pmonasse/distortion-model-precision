#include "modelRational.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>
using namespace libNumerics;

/// Constructor
ModelRational::ModelRational(int order)
: order_(order),
  a_((order+1)*(order+2)/2, 0),
  b_((order+1)*(order+2)/2, 0),
  c_((order+1)*(order+2)/2, 0),
  xi_(order+1), yi_(order+1)
{
    assert(order>=1);
    int end = static_cast<int>(a_.size());
    a_[end-3]=b_[end-2]=c_[end-1]=1; // Init to identity
}

/// Redefined from base class
void ModelRational::change_order(int order) {
    assert(order>=1 && order>=order_);
    order_ = order;
    // Preserve the coefficients with their assigned order, since they are
    // stored in decreasing degree
    std::reverse(a_.begin(), a_.end());
    std::reverse(b_.begin(), b_.end());
    std::reverse(c_.begin(), c_.end());
    a_.resize((order+1)*(order+2)/2, 0);
    b_.resize((order+1)*(order+2)/2, 0);
    c_.resize((order+1)*(order+2)/2, 0);
    std::reverse(a_.begin(), a_.end());
    std::reverse(b_.begin(), b_.end());
    std::reverse(c_.begin(), c_.end());
    xi_.resize(order+1);
    yi_.resize(order+1);
}

/// Compute powers of i, from 0 to order_ (included)
static void powers(flnum x, std::vector<flnum>& xi) {
    xi[0] = 1.0;
    for(size_t i=1; i<xi.size(); i++)
        xi[i] = x*xi[i-1];
}

/// Apply distortion model
void ModelRational::apply(flnum& x, flnum& y) const {
    // Precompute powers of x and y
    powers(x, xi_);
    powers(y, yi_);
    // Computation
    x = y = 0;
    flnum d=0;
    for(int o=order_, n=0; o>=0; o--)
        for(int i=o; i>=0; i--, n++) {
            x += a_[n]*xi_[i]*yi_[o-i];
            y += b_[n]*xi_[i]*yi_[o-i];
            d += c_[n]*xi_[i]*yi_[o-i];
        }
    x /= d;
    y /= d;
}

/// For Levenberg-Marquardt minimization of the radial model, necessary when
/// the center is variable.
class MinRational : public MinLM {
    vector<flnum> ptsU_;
    vector<flnum> ptsD_;
    std::vector<flnum>& xi_;
    std::vector<flnum>& yi_;
    int order_;
public:
    MinRational(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD,
                std::vector<flnum>& xi, std::vector<flnum>& yi, int order);

    flnum eval(std::vector<flnum>& a,
               std::vector<flnum>& b,
               std::vector<flnum>& c);

    virtual void modelData(const vector<flnum>& P,
                           vector<flnum>& ymodel) const;
    virtual void modelJacobian(const vector<flnum>& P,
                               matrix<flnum>& J) const;
};



/// Adjust the coefficients so that ptsU (undistorted) match to ptsD (distorted)
void ModelRational::evaluate(matrix<flnum>& ptsU, matrix<flnum>& ptsD) {
    assert(ptsU.nrow()==2 && ptsD.nrow()==2);
    assert(ptsU.ncol()==ptsD.ncol());
    MinRational(ptsU, ptsD, xi_, yi_, order_).eval(a_, b_, c_);
}

MinRational::MinRational(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD,
                         std::vector<flnum>& xi, std::vector<flnum>& yi,
                         int order)
: ptsU_(2*ptsU.ncol()), ptsD_(2*ptsD.ncol()), xi_(xi), yi_(yi), order_(order)
{
    const int n=ptsU.ncol();
    for(int j=0, k=0; j<n; j++, k+=2) {
        ptsU_(k+0) = ptsU(0,j);
        ptsU_(k+1) = ptsU(1,j);
        ptsD_(k+0) = ptsD(0,j);
        ptsD_(k+1) = ptsD(1,j);
    }
}

/// Evaluation
flnum MinRational::eval(std::vector<flnum>& a,
                        std::vector<flnum>& b,
                        std::vector<flnum>& c) {
    vector<flnum> P( static_cast<int>(a.size()+b.size()+c.size()) );
    std::copy(a.begin(), a.end(), P.begin());
    std::copy(b.begin(), b.end(), P.begin()+a.size());
    std::copy(c.begin(), c.end(), P.begin()+a.size()+b.size());
    float err = minimize(P, ptsD_, 1.0e-17, 1000);
    std::copy(P.begin(), P.begin()+a.size(), a.begin());
    std::copy(P.begin()+a.size(), P.begin()+a.size()+b.size(), b.begin());
    std::copy(P.begin()+a.size()+b.size(), P.end(), c.begin());
    return err/M_SQRT2;
}

/// Vectors a, b, and c are concatenated in P.
void MinRational::modelData(const vector<flnum>& P,
                            vector<flnum>& ymodel) const {
    const int sz = P.nrow()/3;
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        flnum x=ptsU_(i+0), y=ptsU_(i+1), z=0;
        powers(x, xi_);
        powers(y, yi_);
        vector<flnum>::const_iterator a(P.begin()),
                                      b(P.begin()+sz),
                                      c(P.begin()+sz*2);
        x=y=z=0;
        for(int o=order_; o>=0; o--)
            for(int k=o; k>=0; k--) {
                const flnum v = xi_[k]*yi_[o-k];
                x += *a++*v;
                y += *b++*v;
                z += *c++*v;
            }
        ymodel(i+0) = x/z;
        ymodel(i+1) = y/z;
    }
}

/// Vectors a, b, and c are concatenated in P.
void MinRational::modelJacobian(const vector<flnum>& P,
                                matrix<flnum>& J) const {
    const int sz = P.nrow()/3;
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        flnum x=ptsU_(i+0), y=ptsU_(i+1), z=0;
        powers(x, xi_);
        powers(y, yi_);
        vector<flnum>::const_iterator a(P.begin()),
                                      b(P.begin()+sz),
                                      c(P.begin()+sz*2);
        x=y=z=0;
        for(int o=order_, j=0; o>=0; o--)
            for(int k=o; k>=0; k--, j++) {
                const flnum v = xi_[k]*yi_[o-k];
                J(i+0,j) = J(i+1,j+sz) = J(i+0,j+2*sz) = J(i+1,j+2*sz) = v;
                J(i+1,j) = J(i+0,j+sz) = 0;
                x += *a++*v;
                y += *b++*v;
                z += *c++*v;
            }
        for(int j=0; j<sz; j++) {
            J(i+0,j)    /= z;
            J(i+1,j+sz) /= z;
            J(i+0,j+2*sz) *= -x/(z*z);
            J(i+1,j+2*sz) *= -y/(z*z);
        }
    }
}
