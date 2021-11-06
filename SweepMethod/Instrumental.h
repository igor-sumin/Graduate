#pragma once

std::tuple<int, int, int, int> prepareData();
bool checkData(int N, int THREADNUM);
bool isPrime(int num);
vec findDivisors(int num);
void printVec(const vec& a, const std::string& name);
void printMatr(const matr& a, const std::string& name);
vec createVecN(int n);
vec createVecRand(int n);
matr createThirdDiagMatrI(int n);
matr createMatr(int n, double a, double b, double c);
vec calcMatrVecMult(const matr& A, const vec& b);

// подготовка пользовательских данных к параллельным вычислениям
std::tuple<int, int, int, int> prepareData() {
	int N, THREADNUM, blockSize, classicSize;
	
	do {
		// N = 8;
		// THREADNUM = 4;

		std::cout << "Введите размерность (N) = ";
		std::cin >> N;

		std::cout << "Введите количество вычислительных узлов (THREADNUM) = ";
		std::cin >> THREADNUM;

	} while (!checkData(N, THREADNUM));

	blockSize = N / THREADNUM;
	classicSize = (THREADNUM - 1) * 2;

	return std::make_tuple(N, THREADNUM, blockSize, classicSize);
}

// проверка на кратность N и THREADNUM
bool checkData(int N, int THREADNUM) {
	if (N < 7) {
		std::cout << "Размерность " << N << " слишком мала\n";
		std::cout << "Параллельные вычисления для нее не эффективны\n";
		std::cout << "Введите бОльшую размерность N\n\n";

		return false;
	}

	if ((N / THREADNUM) < 3) {
		std::cout << "Алгоритм нерабочий для пропорций размерности с количеством вычислительных узлов\n";
		std::cout << "Введите бОльшую размерность, либо уменьшите количество вычислительных узлов\n\n";

		return false;
	}

	vec div = findDivisors(N);
	if (isPrime(N) || std::find(div.begin(), div.end(), THREADNUM) == div.end()) {
		std::cout << "Невозможно разбить на одинаковые блоки размерность " << N << "\n";
		std::cout << "Введите корректные значения N, THREADNUM\n\n";

		return false;
	}
	
	return true;
}

// является ли число простым
bool isPrime(int num) {
	bool ok = true;
	for (int i = 2; i <= num / 2; ++i) {
		if (num % i == 0) {
			ok = false;
			break;
		}
	}

	return ok;
}

// найти все делители числа
vec findDivisors(int num) {
	vec res;

	for (int i = 2; i <= sqrt(num); i++) {
		if (num % i == 0) {
			if (num / i != i) {
				res.push_back(num / i);
			}

			res.push_back(i);
		}
	}
	
	return res;
}

void printVec(const vec& a, const std::string& name) {
	if (a.size() < 30) {
		std::cout << name << ":\n";
		bool flag = true;
		for (int i = 0; i < a.size(); i++) {
			if (!flag) {
				std::cout << ", ";
			} flag = false;
			std::cout << a[i];
		} std::cout << std::endl;
	}
}

void printMatr(const matr& a, const std::string& name) {
	if (a.size() < 30) {
		bool first = true;

		std::cout << name << ":\n";
		for (int i = 0; i < a.size(); i++) {
			first = true;
			for (int j = 0; j < a[0].size(); j++) {
				if (!first) {
					std::cout << ", ";
				} first = false;

				printf("%8.3f", a[i][j]);
			} printf("\n");
		} std::cout << "\n";
	}
}

vec createVecN(int n) {
	vec a(n);
	std::iota(a.begin(), a.end(), 0);

	return a;
}

vec createVecRand(int n) {
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	vec a(n);

	#pragma omp parallel for
	for (int i = 0; i < n; i++)
		a[i] = gen() % 100;

	return a;
}

matr createThirdDiagMatrI(int n) {
	matr a(n, vec(n));

	#pragma omp parallel for if (n > 500)
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i][i] = 3.;
			a[i][i - 1] = 1.;
			a[i - 1][i] = 2.;
		}
	}

	a[0][0] = 1.; a[n - 1][n - 1] = 1.;
	a[0][1] = -0.5; a[n - 1][n - 2] = -0.5;

	return a;
}

matr createMatr(int n, double a, double b, double c) {
	matr res(n, vec(n));

	#pragma omp parallel for if (n > 500)
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res[i][i] = c;
			res[i][i - 1] = b;
			res[i - 1][i] = a;
		}
	}

	res[0][0] = 1.; res[n - 1][n - 1] = 1.;
	res[0][1] = 0.; res[n - 1][n - 2] = 0.;

	return res;
}

matr createThirdDiagMatrRand(int n) {
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	matr a(n, vec(n));

	a[0][0] = gen() % 100;
#pragma omp parallel for if (n > 500)
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i][i] = gen() % 100;
			a[i][i - 1] = gen() % 100;
			a[i - 1][i] = gen() % 100;
		}
	}

	return a;
}



// Матрично-векторное умножение
vec calcMatrVecMult(const matr& A, const vec& b) {
	int n = A.size();
	vec res(n);

	#pragma omp parallel shared(res, n)
	{
		#pragma omp for
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				res[i] += A[i][j] * b[j];
	}

	return res;
}