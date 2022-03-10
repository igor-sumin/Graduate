#include "main/interfaces/parallel/ParallelInstrumental.h"

#include <ctime>


bool ParallelInstrumental::isPrime(int num) {
    bool ok = true;
    for (int i = 2; i <= num / 2; ++i) {
        if (num % i == 0) {
            ok = false;
            break;
        }
    }

    return ok;
}

vec ParallelInstrumental::findDivisors(int num) {
    vec res;

    for (int i = 2; i <= sqrt(num); i++) {
        if (num % i == 0) {
            if (num / i != i) {
                res.push_back((double)num / i);
            }

            res.push_back(i);
        }
    }

    return res;
}

void ParallelInstrumental::prepareData() {
    size_t n, threadNums;

    do {
        std::cout << "Enter the dimension (N) = ";
        std::cin >> n;

        std::cout << "Enter the number of compute nodes (threadNum) = ";
        std::cin >> threadNums;
    } while(!checkData());

    this->prepareData(threadNums);
}

void ParallelInstrumental::prepareData(size_t tN) {
    threadNum = tN;

    if (!checkData()) {
        throw std::invalid_argument("n = " + std::to_string(N) + ", threadNum = " + std::to_string(threadNum));
    }

    blockSize = (int)N / threadNum;
    interSize = (threadNum - 1) * 2;

    this->setParallelOptions();
}

bool ParallelInstrumental::checkData() const {
    if (N < 7) {
        std::cerr << "Dimension (N = " << N << ") is too small\n"
                  << "Parallel computing is not effective for it\n"
                  << "Enter a larger dimension (N)\n\n";

        return false;
    }

    if ((N / threadNum) < 3) {
        std::cerr << "The algorithm is not working for the proportions of the dimension (N) "
                  << "with the number of computing nodes (threadNum)\n"
                  << "Enter a larger dimension (N), or reduce the number of computing nodes (threadNum)\n\n";

        return false;
    }

    vec div = findDivisors((int)N);
    if (isPrime((int)N) || std::find(div.begin(), div.end(), threadNum) == div.end()) {
        std::cerr << "It is impossible to split the dimension (N = "<< N << ") into the same blocks (blockSize)\n"
                  << "Enter the correct values of dimension (N) and number of computing nodes (threadNum)\n\n";

        return false;
    }

    return true;
}

void ParallelInstrumental::setParallelOptions() const {
    // routine enables or disables dynamic adjustment of the number of threads available
    omp_set_dynamic(0);
    // sets the number of threads to be used by all subsequent parallel regions
    omp_set_num_threads((int)threadNum);
}

void ParallelInstrumental::defineDataByTask7() {
    double total = 12. / (h * h);

    this->A = createThirdDiagMatr(-(2. * total + 5.), total, total, 1., 0., 0., 1.);
    this->b = createVecByTask7();
    this->y = createResByTask7();
}

void ParallelInstrumental::defineDataByNonTask() {
    this->A = createThirdDiagMatr(3., 1., 2., 1, -0.5, -0.5, 1);
    this->b = createVecN();
    this->y.assign(N, 0.);
}

vec ParallelInstrumental::createVecN() {
    vec a(N);
    std::iota(a.begin(), a.end(), 0);

    return a;
}

vec ParallelInstrumental::createVecByTask7() {
    vec res(N);
    size_t i;

    res[0] = 10.; res[N - 1] = 100.;
    #pragma omp parallel for private(i) shared(N, res, x) default(none) if (N > 500)
    for (i = 1; i < N - 1; i++) {
        res[i] = 2110. - 450. * x[i] * x[i];
    }

    return res;
}

vec ParallelInstrumental::createResByTask7() {
    vec res(N);
    size_t i;

    #pragma omp parallel for private(i) shared(N, res, x) default(none) if (N > 500)
    for (i = 0; i < N; i++) {
        res[i] = 10. + 90. * x[i] * x[i];
    }

    return res;
}

matr ParallelInstrumental::createThirdDiagMatr(double diag, double upDiag, double downDiag,
                                               double first, double second,
                                               double preLast, double last) {
    matr res(N, vec(N));

    #pragma omp parallel for shared(N, res, diag, upDiag, downDiag) default(none) if (N > 500)
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i][i] = diag;
            res[i][i - 1] = upDiag;
            res[i - 1][i] = downDiag;
        }
    }

    res[0][0] = first;  res[N - 1][N - 1] = last;
    res[0][1] = second; res[N - 1][N - 2] = preLast;

    return res;
}