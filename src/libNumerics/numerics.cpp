/**
 * @file numerics.cpp
 * @brief Linear algebra basics
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2010-2012 Pascal Monasse, Pierre Moulon
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "numerics.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>

namespace libNumerics {

const flnum MinLM::DEFAULT_RELATIVE_TOL = 1E-3;
const flnum MinLM::DEFAULT_LAMBDA_INIT  = 1E-3;
const flnum MinLM::DEFAULT_LAMBDA_FACT  = 10.0;
const flnum MinLM::EPSILON_KERNEL       = 1E-9;

inline flnum ABS(flnum x)
{ return (x >= 0)? x: -x; }

/// Resolution by LU decomposition with pivot.
bool solveLU(const matrix<flnum>& A, const vector<flnum>& B, vector<flnum>& X)
{
    X = B;
    return solveLU(A, X);
}

/// Replace X by A^{-1}X, by LU solver.
bool solveLU(matrix<flnum> A, vector<flnum>& X)
{
    assert(A.nrow() == A.ncol());
    int	n = A.nrow();
    vector<flnum> rowscale(n); // Implicit scaling of each row
    std::vector<int> permut(n,0); // Permutation of rows

    // Get the implicit scaling information of each row
    for(int i=0; i< n; i++) {
        flnum max = 0.0;
        for(int j=0; j< n; j++) {
            flnum tmp = ABS(A(i,j));
            if (tmp> max)
                max = tmp;
        }
        if(max == 0.0)
            return false;
        rowscale(i) = 1.0/max;
    }
    // Perform the decomposition
    for(int k=0; k < n; k++) {
        // Search for largest pivot element
        flnum max = rowscale(k)*ABS(A(k,k));
        int imax = k;
        for(int i=k+1; i < n; i++) {
            flnum tmp = rowscale(i)*ABS(A(i,k));
            if(tmp > max) {
                max = tmp;
                imax = i;
            }
        }
        if(max == 0.0)
            return false;

        // Interchange rows if needed
        if(k != imax) {
            A.swapRows(k, imax);
            rowscale(imax) = rowscale(k); // Scale of row k no longer needed
        }
        permut[k] = imax; // permut(k) was not initialized before
        flnum Akk = 1/A(k,k);
        for(int i=k+1; i < n; i++) {
            flnum tmp = A(i,k) *= Akk; // Divide by pivot
            for (int j=k+1;j < n; j++) // Reduce the row
                A(i,j) -= tmp*A(k,j);
        }
    }
    // Forward substitution
    for (int k=0; k < n; k++) {
        flnum sum = X(permut[k]);
        X(permut[k]) = X(k);
        for(int j = 0; j < k; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum;
    }
    // Backward substitution
    for(int k=n-1; k >= 0; k--) {
        flnum sum = X(k);
        for(int j=k+1; j < n; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum/A(k,k);
    }
    return true;
}

/// Decompose A into U diag(D) V^T with U(m,m) and V(n,n) orthonormal matrices.
SVD::SVD(const matrix<flnum>& A)
: U(A.nrow(),A.nrow()), V(A.ncol(),A.ncol()), D(std::min(A.nrow(),A.ncol()))
{
    D.fill(0);
    A.SVD(U, D, V);
    sort();
}

/// Recompose from SVD. This should be the initial matrix.
matrix<flnum> SVD::compose() const
{
    return U * D.diag(U.ncol(), V.nrow()) * V.t();
}

/// Return ith greathest singular value. It is safer than picking directly in D
/// as it checks against out of bounds.
flnum SVD::sv(int i) const
{
    assert(0<=i && (i<U.nrow() || i<V.ncol()));
    if(i<D.nrow())
        return D(i);
    return 0;
}

/// Functor to sort elements of a vector by decreasing value.
class SVDElement {
public:
    SVDElement(const vector<flnum>& D, int i)
    : m_val(D(i)), m_i(i) {}
    bool operator<(const SVDElement& e) const
    { return (m_val>e.m_val); }

    flnum m_val;
    int m_i;
};

/// Sort SVD by decreasing order of singular value.
void SVD::sort()
{
    std::vector<SVDElement> vec;
    for(int i=0; i < D.nrow(); i++)
        vec.push_back( SVDElement(D, i) );
    std::sort(vec.begin(), vec.end());
    // Apply permutation
    for(int i=std::min(U.ncol(),V.ncol())-1; i >=0; i--)
        if(vec[i].m_i != i) { // Find cycle of i
            const vector<flnum> colU = U.col(i);
            const vector<flnum> colV = V.col(i);
            const flnum w = D(i);
            int j = i;
            while(vec[j].m_i != i) {
                U.paste(0,j, U.col(vec[j].m_i));
                V.paste(0,j, V.col(vec[j].m_i));
                D(j) = D(vec[j].m_i);
                std::swap(j,vec[j].m_i);
            }
            vec[j].m_i = j;
            U.paste(0,j, colU);
            V.paste(0,j, colV);
            D(j) = w;
        }
}

/// Solve the linear system Ax = 0 with ||x|| = 1.0 via SVD.
/// Check if there is a single singular value such that SV/maxSV<ratioExtremes.
/// However if the two lowest SV satisfy this, accept if minSV/minSV2<ratio2Min.
/// This indicates whether the solution is unique.
bool SVD::Nullspace(const Mat& A, Vec *nullspace,
                    double ratioExtremes, double ratio2Min)
{
    if(A.nrow()+1<A.ncol()) // rank<n-1
        return false;
    libNumerics::SVD svd(A);
    int last=A.lastCol();
    assert(last>=1);
    flnum minSV1=std::max(svd.sv(last),
			  std::numeric_limits<flnum>::epsilon()*svd.sv(0));
    flnum minSV2=std::max(svd.sv(last-1),
			  std::numeric_limits<flnum>::epsilon()*svd.sv(0));    
    if(minSV1>=ratioExtremes*svd.sv(0) || // No low SV
       (minSV2<ratioExtremes*svd.sv(0) && // 2+ low SV...
        minSV1>=ratio2Min*minSV2))        // ... and rather close
        return false;
    *nullspace = svd.V.col(svd.V.lastCol());
    return true;
}

/// Inverse of norm-2 condition value (ratio of extreme singular values)
flnum SVD::InvCond(const Mat& A) {
    assert(A.nrow()==A.ncol());
    libNumerics::SVD svd(A);
    return svd.D(svd.D.lastRow())/svd.D(0);
}

/// Make square matrix of rank<=2.
void SVD::EnforceRank2_3x3(const Mat& A, Mat *ARank)
{
    assert(A.nrow()==3 && A.ncol()==3);
    libNumerics::SVD svd(A);
    svd.D(svd.D.lastRow()) = 0;
    *ARank = svd.compose();
}

/// Save the two last nullspace vectors as 3x3 matrices.
void SVD::Nullspace2_Remap33(const Mat &A, Mat& f1, Mat& f2) {
    assert(A.ncol()==9);
    assert(f1.nrow()==3 && f1.ncol()==3);
    assert(f2.nrow()==3 && f2.ncol()==3);      
    libNumerics::SVD svd(A);
    f1.read( svd.V.col(svd.V.lastCol()-1) );
    f2.read( svd.V.col(svd.V.lastCol()) );
}

/// Constructor.
MinLM::MinLM()
: iterations(0), relativeTol(DEFAULT_RELATIVE_TOL),
  lambdaInit(DEFAULT_LAMBDA_INIT), lambdaFact(DEFAULT_LAMBDA_FACT)
{}

/// In equation JtJ X = B, remove columns of J close to 0, so that JtJ can be
/// invertible
void MinLM::compress(matrix<flnum>& JtJ, vector<flnum>& B)
{
    flnum max=0;
    for(int i=0; i < JtJ.nrow(); i++)
        if(JtJ(i,i) > max)
            max = JtJ(i,i);
    max *= EPSILON_KERNEL;
    m_nullCols.clear();
    for(int i=0; i < JtJ.nrow(); i++)
        if(JtJ(i,i) <= max)
            m_nullCols.push_back(i);
    if( m_nullCols.empty() )
        return;
    int n=(int)m_nullCols.size();
    matrix<flnum> JtJ2(JtJ.nrow()-n, JtJ.ncol()-n);
    vector<flnum> B2(B.nrow()-(int)m_nullCols.size());
    for(int i=0,i2=0; i < JtJ.nrow(); i++) {
        if(i-i2 < n && m_nullCols[i-i2]==i)
            continue;
        for(int j=0,j2=0; j < JtJ.ncol(); j++) {
            if(j-j2 < n && m_nullCols[j-j2]==j)
                continue;
            JtJ2(i2,j2) = JtJ(i,j);
            j2++;
        }
        B2(i2) = B(i);
        i2++;
    }
    swap(JtJ,JtJ2);
    swap(B,B2);
}

/// Insert 0 in rows of B that were removed by \c compress()
void MinLM::uncompress(vector<flnum>& B)
{
    if(m_nullCols.empty())
        return;
    int n=(int)m_nullCols.size();
    vector<flnum> B2(B.nrow()+(int)m_nullCols.size());
    for(int i=0,i2=0; i2 < B2.nrow(); i2++)
        if(i2-i < n && m_nullCols[i2-i]==i2)
            B2(i2)=0;
        else
            B2(i2) = B(i++);
    swap(B,B2);
}


/// Perform minimization.
/// \a targetRMSE is the root mean square error aimed at.
/// Return the reached RMSE. Since the class does not know the dimension, the
/// real RMSE should be this value multiplied by sqrt(dim). For example, for 2-D
/// points this would be sqrt(2) times the returned value.
flnum MinLM::minimize(vector<flnum>& P, const vector<flnum>& yData,
                      flnum targetRMSE, int maxIters)
{
    flnum errorMax = targetRMSE*targetRMSE*yData.nrow();
    vector<flnum> yModel( yData.nrow() );
    modelData(P, yModel);
    vector<flnum> E( yData-yModel );
    flnum error = E.qnorm();
    matrix<flnum> J( yData.nrow(), P.nrow() );
    modelJacobian(P, J);
    matrix<flnum> Jt = J.t();
    matrix<flnum> JtJ = Jt*J;
    vector<flnum> B = Jt*E;
    compress(JtJ, B);

    flnum lambda = lambdaInit;
    for(iterations=0; iterations < maxIters && error > errorMax; iterations++) {
        matrix<flnum> H(JtJ);
        for(int i = 0; i < H.nrow(); i++)
            H(i,i) *= 1+lambda;
        vector<flnum> dP( P.nrow() );
        solveLU(H, B, dP);
        uncompress(dP);
        vector<flnum> tryP = P + dP;
        modelData(tryP, yModel);
        E = yData - yModel;
        flnum tryError = E.qnorm();
        if(ABS(tryError-error) <= relativeTol*error)
            break;
        if(tryError > error)
            lambda *= lambdaFact;
        else {
            lambda /= lambdaFact;
            error = tryError;
            P = tryP;
            modelJacobian(P, J);
            Jt = J.t();
            JtJ = Jt*J;
            B = Jt*E;
            compress(JtJ, B);
        }
    }
    return sqrt(error/yData.nrow());
}

} // namespace libNumerics
