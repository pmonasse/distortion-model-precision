#include "modelPoly.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>
using namespace libNumerics;

/// Constructor
ModelPoly::ModelPoly(int order)
: order_(order),
  a_((order+1)*(order+2)/2, 0),
  b_((order+1)*(order+2)/2, 0)
{
    assert(order>=1);
    int end = static_cast<int>(a_.size());
    a_[end-3]=b_[end-2]=1; // Init to identity
}

/// Redefined from base class
void ModelPoly::change_order(int order) {
    assert(order>=1 && order>=order_);
    order_ = order;
    // Preserve the coefficients with their assigned order, since they are
    // stored in decreasing degree
    std::reverse(a_.begin(), a_.end());
    std::reverse(b_.begin(), b_.end());
    a_.resize((order+1)*(order+2)/2, 0);
    b_.resize((order+1)*(order+2)/2, 0);
    std::reverse(a_.begin(), a_.end());
    std::reverse(b_.begin(), b_.end());    
}

/// Compute powers of i, from 0 to order_ (included)
void ModelPoly::powers(flnum x, flnum* xi) const {
    xi[0] = 1.0;
    for(int i=1; i<=order_; i++)
        xi[i] = x*xi[i-1];
}

/// Apply distortion model
void ModelPoly::apply(flnum& x, flnum& y) const {
    // Precompute powers of x and y
    flnum* xi = new flnum[order_+1];
    flnum* yi = new flnum[order_+1];
    powers(x, xi);
    powers(y, yi);
    // Computation
    x = y = 0;
    for(int o=order_, n=0; o>=0; o--)
        for(int i=o; i>=0; i--) {
            x += a_[n]*xi[i]*yi[o-i];
            y += b_[n]*xi[i]*yi[o-i];
            ++n;
        }
    delete [] xi;
    delete [] yi;
}

/// Adjust the coefficients so that ptsU (undistorted) match to ptsD (distorted)
void ModelPoly::evaluate(matrix<flnum>& ptsU, matrix<flnum>& ptsD) {
    assert(ptsU.nrow()==2 && ptsD.nrow()==2);
    assert(ptsU.ncol()==ptsD.ncol());
    const int eqn=ptsU.ncol();
    const int unknowns=static_cast<int>(a_.size());
    flnum* xi = new flnum[order_+1];
    flnum* yi = new flnum[order_+1];
    matrix<flnum> A(eqn,unknowns+1);
    for(int eq=0; eq<eqn; eq++) {
        powers(ptsU(0,eq), xi);
        powers(ptsU(1,eq), yi);
        for(int o=order_, n=0; o>=0; o--)
            for(int i=o; i>=0; i--)
                A(eq,n++)=xi[i]*yi[o-i];
    }
    delete [] xi;
    delete [] yi;

    for(int eq=0; eq<eqn; eq++)
        A(eq,A.lastCol()) = -ptsD(0,eq);
    vector<flnum> N(A.ncol());
    if(! SVD::Nullspace(A, &N))
        throw std::runtime_error("Ambiguous SVD::Nullspace");
    for(int i=0; i<unknowns; i++)
        a_[i] = N(i) / N(unknowns);

    for(int eq=0; eq<eqn; eq++)
        A(eq,A.lastCol()) = -ptsD(1,eq);
    if(! SVD::Nullspace(A, &N))
        throw std::runtime_error("Ambiguous SVD::Nullspace");
    for(int i=0; i<unknowns; i++)
        b_[i] = N(i) / N(unknowns);
}
