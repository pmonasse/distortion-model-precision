#include "modelRadial.h"
#include <cassert>
#include <stdexcept>
using namespace libNumerics;

/// Constructor
ModelRadial::ModelRadial(int order, bool varCenter, bool onlyOddOrder)
: k_(order, 0), varCenter_(varCenter), odd_(onlyOddOrder) {
    assert(order>=1);
    k_[0] = 1;
    center_[0]=center_[1]=0;
}

/// Change order, normally used for expansion
void ModelRadial::change_order(int order) {
    assert(order>=1);
    k_.resize(order, 0);
}

void ModelRadial::apply(flnum& x, flnum& y) const {
    x -= center_[0];
    y -= center_[1];
    const flnum r1 = sqrt(x*x+y*y);
    flnum r2_r1=0; // r2/r1, polynomial of degree order-1 in r1
    for(std::vector<flnum>::const_reverse_iterator it=k_.rbegin();
        it!=k_.rend(); ++it)
        r2_r1 = *it + r1*r2_r1;
    x = center_[0] + x*r2_r1;
    y = center_[1] + y*r2_r1;
}

/// For Levenberg-Marquardt minimization of the radial model, necessary when
/// the center is variable.
class MinRadial : public MinLM {
    vector<flnum> ptsU_;
    vector<flnum> ptsD_;
    bool odd_;
public:
    MinRadial(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD);

    flnum eval(std::vector<flnum>& k, flnum center[2], bool odd);

    virtual void modelData(const vector<flnum>& P,
                           vector<flnum>& ymodel) const;
    virtual void modelJacobian(const vector<flnum>& P,
                               matrix<flnum>& J) const;
};

// Copy even index elements of v1 into v2
static void copy_step2_in(const std::vector<flnum>& v1, vector<flnum>& v2) {
    std::vector<flnum>::const_iterator in=v1.begin();
    vector<flnum>::iterator out=v2.begin();
    for(bool ok=true; in!=v1.end(); ++in, ok=!ok)
        if(ok) {
            assert(out != v2.end());
            *out++ = *in;
        }
}

// Copy v1 into even index elements of v2
static void copy_step2_out(const vector<flnum>& v1, std::vector<flnum>& v2) {
    vector<flnum>::const_iterator in=v1.begin();
    std::vector<flnum>::iterator out=v2.begin();
    for(; in!=v1.end(); ++in, ++out) {
        assert(out != v2.end());
        *out++ = *in;
        if(out == v2.end()) // Do not go beyond v2.end()
            break;
    }
}

/// Direct evaluation if fixed center, LM otherwise.
void ModelRadial::evaluate(matrix<flnum>& ptsU, matrix<flnum>& ptsD) {
    if(varCenter_) // Levenberg-Marquardt if variable center
        MinRadial(ptsU, ptsD).eval(k_, center_, odd_);
    else { // Direct solution
        const int eqn=ptsU.ncol();
        const int unknowns=static_cast<int>(odd_? (k_.size()+1)/2: k_.size());
        libNumerics::matrix<flnum> A(eqn,unknowns+1);
        for(int eq=0; eq<eqn; eq++) {
            flnum x=ptsD(0,eq)-center_[0], y=ptsD(1,eq)-center_[1];
            A(eq,0) = -sqrt(x*x+y*y);
            x=ptsU(0,eq)-center_[0]; y=ptsU(1,eq)-center_[1];
            flnum ri = sqrt(x*x+y*y);
            const flnum s = odd_? ri*ri: ri;
            for(int i=0; i<unknowns; i++, ri*=s)
                A(eq,i+1) = ri;
        }
        libNumerics::vector<flnum> N(A.ncol());
        if(! SVD::Nullspace(A, &N))
            throw std::runtime_error("Ambiguous SVD::Nullspace");
        N = N.copy(1,N.lastRow())/N(0);
        if(odd_) {
            std::fill(k_.begin(), k_.end(), (flnum)0);
            copy_step2_out(N, k_);
        } else
            std::copy(N.begin(), N.end(), k_.begin());
    }
}

/// Constructor
MinRadial::MinRadial(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD)
: ptsU_(2*ptsU.ncol()), ptsD_(2*ptsD.ncol()), odd_(false) {
    assert(ptsU.nrow()==2 && ptsD.nrow()==2);
    assert(ptsU.ncol()==ptsD.ncol());
    const int n=ptsU.ncol();
    for(int j=0, k=0; j<n; j++, k+=2) {
        ptsU_(k+0) = ptsU(0,j);
        ptsU_(k+1) = ptsU(1,j);
        ptsD_(k+0) = ptsD(0,j);
        ptsD_(k+1) = ptsD(1,j);
    }
}

/// Evaluation
flnum MinRadial::eval(std::vector<flnum>& k, flnum center[2], bool odd) {
    odd_ = odd;
    int nCoeff = static_cast<int>(odd? k.size()/2: k.size());
    vector<flnum> P(nCoeff + 2); // +2: center of distortion
    if(odd) copy_step2_in(k, P);
    else std::copy(k.begin(), k.end(), P.begin());
    std::copy(center, center+2, P.begin()+nCoeff);
    float err = minimize(P, ptsD_, 1.0e-17, 1000);
    if(odd) copy_step2_out(P, k);
    else std::copy(P.begin(), P.begin()+nCoeff, k.begin());
    std::copy(P.begin()+nCoeff, P.end(), center);
    return err/M_SQRT2;
}

/// P contains the coefficients of the polynomial followed by the center coords.
void MinRadial::modelData(const vector<flnum>& P,
                          vector<flnum>& ymodel) const {
    const int order = P.nrow()-2;
    const float cx=P(order), cy=P(order+1);
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0)-cx, y=ptsU_(i+1)-cy;
        const flnum r = odd_? x*x+y*y: sqrt(x*x+y*y);
        flnum r2_r1=0; // r2/r1, polynomial of degree order-1 in r1
        for(int j=order-1; j>=0; j--)
            r2_r1 = P(j) + r*r2_r1;
        ymodel(i+0) = cx + x*r2_r1;
        ymodel(i+1) = cy + y*r2_r1;
    }
}

/// P contains the coefficients of the polynomial followed by the center coords.
void MinRadial::modelJacobian(const vector<flnum>& P,
                              matrix<flnum>& J) const {
    const int order = P.nrow()-2;
    const float cx=P(order), cy=P(order+1);
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0)-cx, y=ptsU_(i+1)-cy;
        const flnum r = odd_? x*x+y*y: sqrt(x*x+y*y);
        const flnum d = (x*x+y*y>1e-10? 1/(x*x+y*y): 0);
        flnum dr2x_x=0; // derivative of x*r2/r1 wrt x
        flnum dr2x_y=0; // derivative of x*r2/r1 wrt y, =deriv of y*r2/r1 wrt x
        flnum dr2y_y=0; // derivative of y*r2/r1 wrt y
        flnum rj=1;
        for(int j=0; j<order; j++, rj*=r) {
            J(i+0,j) = x*rj;
            J(i+1,j) = y*rj;
            int k=odd_? 2*j: j;
            dr2x_x += P(j)*rj*(1+k*x*x*d);
            dr2x_y += P(j)*rj*   k*x*y*d;
            dr2y_y += P(j)*rj*(1+k*y*y*d);
        }
        J(i+0,order+0)= 1-dr2x_x;
        J(i+1,order+0)=  -dr2x_y;
        J(i+0,order+1)=  -dr2x_y;
        J(i+1,order+1)= 1-dr2y_y;
    }
}
