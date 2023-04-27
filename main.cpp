#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <cmath>

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

using namespace std;


class Matrix {
public:
    int rows, cols;
    vector<vector<double> > matr;

    Matrix(int r, int c) : rows(r), cols(c), matr(r, vector<double>(c, 0)) {}

    void set_value(int row, int col, double value) {
        matr[row][col] = value;
    }

    double get_value(int row, int col) {
        return matr[row][col];
    }

    int get_rows() {
        return rows;
    }

    int get_cols() {
        return cols;
    }

    friend ostream &operator<<(ostream &os, const Matrix &rhs) {
        os.precision(2);
        for (int i = 0; i < rhs.rows; ++i) {
            for (int j = 0; j < rhs.cols; ++j) {
                printf("%.4f ", rhs.matr[i][j]);
                //os << rhs.matr[i][j] << " ";

            }
            os << "\n";
        }
        return os;
    }

    friend istream &operator>>(istream &stream, Matrix &rhs) {
        for (int i = 0; i < rhs.rows; ++i) {
            for (int j = 0; j < rhs.cols; ++j) {
                stream >> rhs.matr[i][j];
            }
        }
        return stream;
    }

    Matrix operator+(Matrix &b) {
        Matrix result(get_rows(), get_cols());
        for (int i = 0; i < result.get_rows(); ++i) {
            for (int j = 0; j < result.get_cols(); ++j) {
                result.set_value(i, j, get_value(i, j) + b.get_value(i, j));
            }
        }
        return result;
    }

    Matrix transpose() {
        Matrix transpose_matrix(cols, rows);
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                transpose_matrix.set_value(i, j, matr[j][i]);
            }
        }
        return transpose_matrix;
    }

    Matrix operator-(Matrix &b) {
        Matrix result(get_rows(), get_cols());
        for (int i = 0; i < result.get_rows(); ++i) {
            for (int j = 0; j < result.get_cols(); ++j) {
                result.set_value(i, j, get_value(i, j) - b.get_value(i, j));
            }
        }
        return result;
    }

    Matrix operator*(Matrix &b) {
        if (get_cols() != b.get_rows()) {
            cout << "Error: the dimensional problem occurred\n";
        }
        Matrix result(get_rows(), b.get_cols());
        for (int write_i = 0; write_i < result.get_rows(); ++write_i) {
            for (int write_j = 0; write_j < result.get_cols(); ++write_j) {
                for (int now = 0; now < get_cols(); ++now) {
                    result.set_value(write_i, write_j, result.get_value(write_i, write_j) +
                                                       get_value(write_i, now) * b.get_value(now, write_j));
                }
            }
        }
        return result;
    }

    void swapRows(int r1, int r2) {
        if (r1 < 0 || r1 >= rows || r2 < 0 || r2 >= rows) {
            throw std::out_of_range("Invalid row index");
        }
        if (r1 == r2) {
            return;
        }
        for (int c = 0; c < cols; ++c) {
            swap(matr[r1 * cols + c], matr[r2 * cols + c]);
        }
    }


};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int r) : Matrix(r, r) {};

    double determinant() {
        int step = 0;

        Matrix mat(*this);

        double det = 1.0;

        for (int i = 0; i < cols; i++) {
            int max_row = i;
            for (int j = i + 1; j < cols; j++) {
                if (abs(matr[j][i]) > abs(matr[max_row][i])) {
                    max_row = j;
                }
            }


            if (max_row != i) {
                step++;
                //cout << "step #" << step << ": permutation\n";
                swap(matr[i], matr[max_row]);
                //cout << (*this);
                det = -det;
            }

            if (matr[i][i] == 0) {
                //cout << "result:\n";
                return 0.0;
            }


            for (int j = i + 1; j < cols; j++) {
                double factor = matr[j][i] / matr[i][i];
                for (int k = 0; k < cols; ++k) {
                    matr[j][k] = matr[j][k] - matr[i][k] * factor;
                }
                step++;
                //cout << "step #" << step << ": elimination\n";
                //cout << (*this);
            }


            det *= matr[i][i];
        }
        //cout << "result:\n";
        return det;
    }

    SquareMatrix invert() {
        int n = get_rows();
        Matrix augment(n, n * 2);
        int step = 0;

        for (int i = 0; i < n; i++){
            for (int j = 0; j < 2 * n; j++) {
                if (j < n) {
                    augment.set_value(i, j, get_value(i, j));
                } else {
                    augment.set_value(i, j, i == j - n ? 1 : 0);
                }
            }
        }
        //cout<< "step #0: Augmented Matrix\n";
        //cout<< augment;
        //cout<< "Direct way:\n";
        for (int i = 0; i < augment.rows; i++) {
            int max_row = i;
            for (int j = i + 1; j < augment.rows; j++) {
                if (abs(augment.matr[j][i]) > abs(augment.matr[max_row][i])) {
                    max_row = j;
                }
            }


            if (max_row != i) {
                step++;
                //cout << "step #" << step << ": permutation\n";
                swap(augment.matr[i], augment.matr[max_row]);
                //cout << augment;
            }

            if (augment.matr[i][i] == 0) {
                //cout << "error\n";
                return 0.0;
            }


            for (int j = i + 1; j < augment.rows; j++) {
                if (augment.matr[j][i] != 0){
                    double factor = augment.matr[j][i] / augment.matr[i][i];
                    for (int k = 0; k < augment.cols; ++k) {
                        augment.matr[j][k] = augment.matr[j][k] - augment.matr[i][k] * factor;
                    }
                    step++;
                    //cout << "step #" << step << ": elimination\n";
                    //cout << augment;
                }
            }
        }
        //cout<< "Way back:\n";
        for (int i = augment.rows-1; i >0; i--) {


            if (augment.matr[i][i] == 0) {
                //cout << "error\n";
                return 0.0;
            }


            for (int j = i - 1; j >= 0; j--) {
                if (augment.matr[j][i] != 0){
                    double factor = augment.matr[j][i] / augment.matr[i][i];
                    for (int k = 0; k < augment.cols; ++k) {
                        augment.matr[j][k] = augment.matr[j][k] - augment.matr[i][k] * factor;
                    }
                    step++;
                    //cout << "step #" << step << ": elimination\n";
                    //cout << augment;
                }
            }
        }
        //cout<< "Diagonal normalization:\n";

        for (int i = 0; i < n; i++) {
            double val = augment.get_value(i, i);
            for (int j = 0; j < augment.get_cols(); j++) {
                augment.set_value(i, j, augment.get_value(i, j) / val);
            }
        }
        //cout<< augment;
        SquareMatrix invert(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                invert.set_value(i, j, augment.get_value(i, j+n));
            }
        }
        //cout<<"result:\n";
        return invert;
    }
};

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int r) : Matrix(r, r) {
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < r; j++) {
                if (i == j) {
                    this->matr[i][j] = 1;
                } else {
                    this->matr[i][j] = 0;
                }
            }
        }
    }
};


int main() {
    FILE *pipe = popen(R"(C:\gnuplot\bin\gnuplot -persist)", "w");
    int m;
    cin >> m;

    vector<pair<double, double>> data(m);
    for (int i = 0; i < m; i++) {
        cin >> data[i].first >> data[i].second;
    }

    int n;
    cin >> n;

    Matrix A(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            A.set_value(i, j, pow(data[i].first, j));
        }
    }

    Matrix AT = A.transpose();
    Matrix ATA = AT * A;

    SquareMatrix ATA_square(ATA.get_rows());
    for (int i = 0; i < ATA.get_rows(); ++i) {
        for (int j = 0; j < ATA.get_cols(); ++j) {
            ATA_square.set_value(i, j, ATA.get_value(i, j));
        }
    }
    SquareMatrix ATA_inv = ATA_square.invert();

    Matrix b(m, 1);
    for (int i = 0; i < m; ++i) {
        b.set_value(i, 0, data[i].second);
    }

    Matrix ATb = AT * b;

    Matrix x = ATA_inv * ATb;

    fprintf(pipe, "plot [-2 : 8] [0 : 11] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' with points\n",
            x.matr[3][0], x.matr[2][0], x.matr[1][0], x.matr[0][0]);

    for (int i = 0; i < m; ++i) {
        fprintf(pipe, "%lf %lf\n", data[i].first, data[i].second);
    }

    fprintf(pipe, "e\n");
    pclose(pipe);

    return 0;
}


