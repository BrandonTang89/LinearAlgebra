#include "Matrix.h"

#include <iomanip>

#include "Vec.h"
#include "helpers.h"

// == Constructors ==
// Construct a Matrix from a Vec
Matrix::Matrix(const Vec& v) : H(v.N), W(1) {
    data = std::vector<std::vector<double>>(H, std::vector<double>(1, 0.0));
    for (int i = 0; i < H; ++i) {
        data[i][0] = v.data[i];
    }
}

// Construct a Zero matrix
Matrix::Matrix(int h, int w) : H(h), W(w), data(h, std::vector<double>(w, 0.0)) {}
Matrix::Matrix(): H(0), W(0), data(0) {}

// Matrix Move Constructor
Matrix::Matrix(Matrix&& A) noexcept : H(A.H), W(A.W), data(std::move(A.data)) {}

// Matrix copy constructor
Matrix::Matrix(const Matrix& A) : H(A.H), W(A.W) {
    data = std::vector<std::vector<double>>(H, std::vector<double>(W, 0.0));
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            data[i][j] = A.data[i][j];
        }
    }
}

// Matrix Copy Assignment Operator
Matrix& Matrix::operator=(const Matrix& A) {
    if (this == &A) return *this;
    H = A.H;
    W = A.W;
    data = std::vector<std::vector<double>>(H, std::vector<double>(W, 0.0));
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            data[i][j] = A.data[i][j];
        }
    }
    return *this;
}

// Matrix Move Assignment Operator
Matrix& Matrix::operator=(Matrix&& A) noexcept {
    if (this == &A) return *this;
    H = A.H;
    W = A.W;
    data = std::move(A.data);
    return *this;
}

// Matrix equality operator using isCloseTo
bool Matrix::operator==(const Matrix& B) const {
    if (H != B.H || W != B.W) return false;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            if (!isCloseTo(data[i][j], B.data[i][j])) return false;
        }
    }
    return true;
}

// Matrix Addition operator
Matrix Matrix::operator+(const Matrix& B) const {
    if (H != B.H || W != B.W) throw std::invalid_argument("Matrix dimensions do not match for addition.");
    Matrix result(H, W);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            result.data[i][j] = data[i][j] + B.data[i][j];
        }
    }
    return result;
}

// Matrix subtraction operator
Matrix Matrix::operator-(const Matrix& B) const {
    if (H != B.H || W != B.W) throw std::invalid_argument("Matrix dimensions do not match for subtraction.");
    Matrix result(H, W);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            result.data[i][j] = data[i][j] - B.data[i][j];
        }
    }
    return result;
}

// Matrix multiplication operator
Matrix Matrix::operator*(const Matrix& B) const {
    if (W != B.H) throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    Matrix result(H, B.W);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < B.W; ++j) {
            for (int k = 0; k < W; ++k) {
                result.data[i][j] += data[i][k] * B.data[k][j];
            }
        }
    }
    return result;
}

// Matrix-Vector multiplication operator
Vec Matrix::operator*(const Vec v) const {
    if (W != v.N) throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    Vec result(H);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            result[i] += data[i][j] * v[j];
        }
    }
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(W, H);
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            result.data[j][i] = data[i][j];
        }
    }
    return result;
}

// Get the diagonal of the matrix
std::vector<double> Matrix::diag() const {
    std::vector<double> diag(std::min(H, W));
    for (int i = 0; i < diag.size(); i++) {
        diag[i] = data[i][i];
    }
    return diag;
}

// Matrix stream insertion operator
std::ostream& operator<<(std::ostream& os, const Matrix& A) {
    os << "Matrix of size " << A.H << " x " << A.W << std::endl;
    os << std::fixed << std::setprecision(5);
    for (int i = 0; i < A.H; i++) {
        os << "[ ";
        for (int j = 0; j < A.W; j++) {
            os << A.data[i][j] << " ";
        }
        os << " ]" << std::endl;
    }
    return os;
}

// Matrix scaling operator
Matrix operator*(double c, const Matrix& A) {
    Matrix result(A.H, A.W);
    for (int i = 0; i < A.H; i++) {
        for (int j = 0; j < A.W; j++) {
            result.data[i][j] = c * A.data[i][j];
        }
    }
    return result;
}


// Transpose function
Matrix transpose(const Matrix& A) {
    Matrix result(A.W, A.H);
    for (int i = 0; i < A.H; ++i) {
        for (int j = 0; j < A.W; ++j) {
            result.data[j][i] = A.data[i][j];
        }
    }
    return result;
}

Matrix identity(int n) {  // Returns the identity matrix of size n
    Matrix I(n, n);
    for (int i = 0; i < n; i++) {
        I.data[i][i] = 1;
    }
    return I;
}

bool isOrthogonal(const Matrix& Q) {  // checks if the COLUMNS of Q are orthogonal
    Matrix Qt = transpose(Q);
    Matrix I = Qt * Q;
    if (I != identity(Q.W)) return false;

    // check if the columns of Q are of unit length
    for (int i = 0; i < Q.W; i++) {
        double sum = 0;
        for (int j = 0; j < Q.H; j++) {
            sum += Q.data[j][i] * Q.data[j][i];
        }
        if (!isCloseTo(sum, 1)) return false;
    }
    return true;
}

bool isUpperTriangular(const Matrix& A) {  // check if the matrix is upper triangular
    for (int i = 0; i < A.H; i++) {
        for (int j = 0; j < std::min(i, A.W); j++) {
            if (!isCloseTo(A.data[i][j], 0)) {
                return false;
            }
        }
    }
    return true;
}

bool isLowerTriangular(const Matrix& A) {  // check if the matrix is lower triangular
    for (int i = 0; i < A.H; i++) {
        for (int j = i + 1; j < A.W; j++) {
            if (!isCloseTo(A.data[i][j], 0)) {
                return false;
            }
        }
    }
    return true;
}
