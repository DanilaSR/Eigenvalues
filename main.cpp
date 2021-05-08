#include <iostream>
#include "Matrix.h"
#include "Polynom.cpp"
using namespace :: std;

constexpr double rel_error {1.0e-14};    // smallest relative error we want
constexpr int max_pow = 32;              // max power of 10 we wish to search to
constexpr int max_it = 800;              // max number of iterations
constexpr double small_enough {1.0e-12}; // a coefficient smaller than


struct gap{
    double a;
    double b;
};

Polynom GCD(const Polynom p, const Polynom q){

}

double pow (double a, int b){      //функция возведения в степень
    for (int i = 1; i < b; i++)
        a *= a;
    return a;
}


double abs_ (double a){            //функция взятия модуля
    if (a >= 0)
        return a;
    else
        return -a;
}

double bin_search(gap g, const Polynom& p){

    double fa, fc, eps, mid;
    fa = p.value(g.a);
    eps = 1.0e-8;
    if (p.value(g.b) * p.value(g.a) > 0)
        return 99999;
    while (true) {
        mid = (g.a + g.b) / 2;
        if (g.a - g.b >= 0){
            if (g.a - g.b < eps)
                return mid;
        }
        else{
            if (g.b - g.a < eps)
                return mid;
        }
        fc = p.value(mid);
        if (fc >= 0){
            if (fc < eps)
                return mid;
        }
        else{
            if (-fc < eps)
                return mid;
        }
        if (fc * fa < 0) {
            g.b = mid;
        }
        else {
            g.a = mid;
            fa = fc;
        }
    }
}

//int num_changes(const Matrix<Polynom>& shturm_seq, const double a )
///*
//  Return the number of sign changes in the Sturm sequence in  sseq at the value
//  a.
// */
//{
//    double f, lf;
//    Polynom *s;
//    int changes = 0;
//
//    lf = eval_poly(sturm_seq[0].ord, sturm_seq[0].coef, a);
//
//    lf = shturm_seq.data[0][0].value(a);
//
//    for (s = sturm_seq + 1; s <= sturm_seq + num_poly; s++) {
//        f = eval_poly(s->ord, s->coef, a);
//        if (lf == 0.0 || lf * f < 0) ++changes;
//        lf = f;
//    }
//
//    return changes;
//}

int modp(Polynom *u, Polynom *v, Polynom *r)
/*
  Calculates the modulus of u(x) / v(x) leaving it in r, it  returns 0 if r(x)
  is a constant. Note: this function assumes the leading coefficient of v is 1
  or -1
 */
{
    int k, j;
    double *nr, *end, *uc;

    nr = r->koef;
    end = &u->koef[u->n];

    uc = u->koef;
    while (uc <= end)
        *nr++ = *uc++;

    if (v->koef[v->n] < 0.0) {
        for (k = u->n - v->n - 1; k >= 0; k -= 2)
            r->koef[k] = -r->koef[k];
        for (k = u->n - v->n; k >= 0; k--)
            for (j = v->n + k - 1; j >= k; j--)
                r->koef[j] = -r->koef[j] - r->koef[v->n + k] * v->koef[j - k];
    }
    else {
        for (k = u->n - v->n; k >= 0; k--)
            for (j = v->n + k - 1; j >= k; j--)
                r->koef[j] -= r->koef[v->n + k] * v->koef[j - k];
    }

    k = v->n - 1;
    while (k >= 0 && abs(r->koef[k]) < small_enough) {
        r->koef[k] = 0.0;
        k--;
    }

    r->n = (k < 0) ? 0 : k;

    return r->n;
}

int build_sturm_rest(Matrix<Polynom>& Shturm)
/*
  Build up a sturm sequence for a polynomial in smat, returning the number of
  polynomials in the sequence
 */
{
    double f, *fp, *fc;
    Polynom *sp;

    // Construct the rest of the Sturm sequence
    for (sp = &Shturm.data[2][0]; modp(&Shturm.data[0][0], &Shturm.data[1][0], &Shturm.data[2][0]); sp++) {
        // Reverse the sign and normalise
        f = -abs(sp->koef[sp->n]);
        for (fp = &sp->koef[sp->n]; fp >= sp->koef; fp--)
            *fp /= f;
    }

    sp->koef[0] = -sp->koef[0]; // reverse the sign

    return Shturm.getN();
}



int main() {
    Polynom x(2), q(3), p;
    x.koef[0] = 1;
    x.koef[1] = 0;                               //Задаём х

    Matrix<Polynom> M(3, 3);
    M.data[0][0] = 4; M.data[0][1] = 1; M.data[0][2] = -1;
    M.data[1][0] = 2; M.data[1][1] = 5; M.data[1][2] = -2;
    M.data[2][0] = 4; M.data[2][1] = 4; M.data[2][2] = -1;

    cout << "Исходная матрица:" << endl;
    cout << M << endl;
    for (int i = 0; i < M.getN(); i++){
        for (int j = 0; j < M.getN(); j++){
            if (i == j){
                M.data[i][j] = M.data[i][j]  - x;
            }
        }
    }
    p = M.determinant();
    cout << "Детерминант характеристического уравнения:" << endl;
    cout << p << endl;

//    int size = p.get_size();
//    gap* g = new gap [size];    //Массив для сохранения отрезков, на которых находятся только 1 корень
//    int n = 0;                         //Счётчиек для подсчёта места в массиве
//
//    while (p.value(0) == 0){    //Избавляемся от нулевых корней
//        g[n].a = 0;
//        g[n].b = 0;
//        n++;
//        p = p/x;
//    }
//
//    Matrix<Polynom> shturm_seq (size + 1 , 1);
//    shturm_seq.data[0][0] = p;
//    shturm_seq.data[1][0] = p.derivative();
//
//    cout << "Последовательность штурма:" << endl;
//    //Дихотомия
//    int N;
//    double sum = 0;
//    for (int i = 0; i < p.get_size(); i++){
//        sum += abs_(p.koef[i]);
//    }
//    N = p.degree();
//    double delta;
//    delta = 1.73/(pow(sum, N - 1)*pow(N, -(N+2)/2));
//    cout << delta << endl;
//    double f = -999.0;
//    while (f < 999.00){
//        double a, b;
//        a = f;
//        b = f + delta;
//        if (p.value(a)*p.value(b) < 0){
//            g[n].a = a;
//            g[n].b = b;
//            n++;
//        }
//        f += delta;
//    }
//    cout << "Собственные значения:"<< endl;
//    for (int i = 0; i < n; i++){
//        cout << bin_search(g[i], p) << endl;
//    }
//
//    delete[] g;

}
