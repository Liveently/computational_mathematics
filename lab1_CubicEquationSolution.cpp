#include <iostream>
#include <vector>
#include <cmath>

#define __infinity -1 // ---inf
#define infinity__ 1 // +++inf
#define closed_interval 0
#define inf_val 42

long double a, b, c;
long double delta, eps;

long double f(long double x) { // f(x)
    return (x * x * x) + (a * x * x) + (b * x) + c;
}

long double df_dx(long double x) { // f'(x)
    return 3 * x * x + 2 * a * x + b;
}

long double d2f_dx2(long double x) { // f''(x)
    return 6 * x + 2 * a;
}

static long double D_dfdx() { // D(f')
    return 4 * a * a - 12 * b;
}

struct interval {
    long double l;
    long double r;
    bool left_Inf;
    bool right_Inf;

    interval(long double l, long double r, int limitIndicator) {
        this->l = l;
        this->r = r;

        switch (limitIndicator) {
            case (__infinity) : {
                left_Inf = true;
                right_Inf = false;
                break;
            }
            case (infinity__) : {
                left_Inf = false;
                right_Inf = true;
                break;
            }
            case (closed_interval) : {
                left_Inf = false;
                right_Inf = false;
            }
        }
    }

    bool isUnlimited() const {
        return right_Inf || left_Inf;
    }

    bool isPoint() const {
        return l == r;
    }
};

void limitInterval(interval &intervalLocal) { // ограничиваем интервал
    if (intervalLocal.left_Inf) { // (-inf, t)
        while (intervalLocal.isUnlimited()) {
            long double tmp = intervalLocal.r - delta;
            if (f(tmp) <= -eps) {
                intervalLocal.l = tmp;
                intervalLocal.left_Inf = false;
            } else if (f(tmp) >= eps) {
                intervalLocal.r = tmp;
            } else {
                intervalLocal.l = tmp;
                intervalLocal.r = tmp;
                intervalLocal.left_Inf = false;
            }
        }
    } else { // (t, inf+)
        while (intervalLocal.isUnlimited()) {
            long double tmp = intervalLocal.l + delta;
            if (f(tmp) <= -eps) {
                intervalLocal.l = tmp;
            } else if (f(tmp) >= eps) {
                intervalLocal.r = tmp;
                intervalLocal.right_Inf = false;
            } else {
                intervalLocal.l = tmp;
                intervalLocal.r = tmp;
                intervalLocal.right_Inf = false;
            }
        }
    }
}

long double findRoot(interval intervalLocal) {
    if (intervalLocal.isPoint()) { //если интервал - точка
        return intervalLocal.r;
    }

    if (intervalLocal.isUnlimited()) { //если интервал с бесконечностью - другая функция
        limitInterval(intervalLocal);
    }

    long double tmp = (intervalLocal.l + intervalLocal.r) / 2;

    if (f(intervalLocal.l) <= -eps && f(intervalLocal.r) >= eps) {
        while (std::abs(f(tmp)) >= eps) {
            if (f(tmp) <= -eps) {
                intervalLocal.l = tmp;
            } else {
                intervalLocal.r = tmp;
            }
            tmp = (intervalLocal.l + intervalLocal.r) / 2;
        }
    } else if (f(intervalLocal.l) >= eps && f(intervalLocal.r) <= -eps) {
        while (std::abs(f(tmp)) >= eps) {
            if (f(tmp) <= -eps) {
                intervalLocal.r = tmp;
            } else {
                intervalLocal.l = tmp;
            }
            tmp = (intervalLocal.l + intervalLocal.r) / 2;
        }
    } else {
        perror("\nERROR interval\n");
        printf("\n l = %Lg, r = %Lg\n", intervalLocal.l, intervalLocal.r);
        exit(1);
    }

    return tmp;
}



std::vector<interval> localizeRoots() {

    if (D_dfdx() <= -eps) {   // одно пересечение (функция или только возрастает или только убывает)
        if (f(0) >= eps) {
            return {{inf_val, 0, __infinity}};
        } else if (f(0) <= -eps) {
            return {{0, inf_val, infinity__}};
        } else {
            return {{0, 0, closed_interval}};
        }
    } else if (D_dfdx() < eps) { // дискриминант бизок к 0 - одно пересечение

        long double x = (-2 * a) / 6; // единственный экстремум производной
        if (f(x) >= eps) {
            return {{inf_val, x, __infinity}};
        } else if (f(x) <= -eps) {
            return {{x, inf_val, infinity__}};
        } else {
            return {{x, x, closed_interval}};
        }
    } else { // 2 корня у функции производной
        long double x1 = ((-2 * a) + std::sqrt(D_dfdx())) / 6;
        long double x2 = ((-2 * a) - std::sqrt(D_dfdx())) / 6;
        if (x1 > x2) std::swap(x1, x2);

        long double f_x1 = f(x1);
        long double f_x2 = f(x2);

        if (f_x1 >= eps && f_x2 <= -eps) { // случай 3х корней
            return {{inf_val, x1, __infinity},
                    {x1,      x2, closed_interval},
                    {x2, inf_val, infinity__}};
        } else if (std::abs(f_x1) < eps && f_x2 <= -eps) { // 2 корня сл 1
            return {{x1, x1,      closed_interval},
                    {x2, inf_val, infinity__}};
        } else if (f_x1 >= eps && std::abs(f_x2) < eps) { // 2 корня сл 2
            return {{inf_val, x1, __infinity},
                    {x2,      x2, closed_interval}};
        } else if (f_x1 <= -eps && f_x2 <= -eps) { // случай 1-го корня
            return {{x2, inf_val, infinity__}};
        } else if (f_x1 >= eps && f_x2 >= eps) {
            return {{inf_val, x1, __infinity}};
        } else {
            perror("\nERROR\n");
            exit(1);
        }
    }
}


int main() {

    std::cin>>a>>b>>c;
    std::cin>>delta>>eps;

    std::vector<interval> intervalsLocalization = localizeRoots();


    for (auto &intervalLocal: intervalsLocalization) {
        long double root = findRoot(intervalLocal);
        std::cout<<"Root: "<<root;
    }
}
