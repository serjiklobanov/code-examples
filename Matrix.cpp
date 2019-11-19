#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;

unsigned int gcd(uint32_t a, uint32_t b) {
	while (b) {
		uint32_t r = a % b;
		a = b;
		b = r;
	}
	return a;
}

class RationalDivisionByZero : public std::exception {
};

class Rational {
private:
	constexpr static size_t MAX_LENGTH = 1000;
	int32_t p;
	uint32_t q; // is maintained to be positive

	void reduce() {
		int32_t t = gcd((uint32_t)abs(p), q);
		p /= t;
		q /= t;// fraction reducing
	}

public:
	Rational() : p(0), q(1) {}
	Rational(int P, int Q) {
		if (Q < 0) {
			this->q = -Q;
			this->p = -P;
		}
		else {
			this->q = Q;
			this->p = P;
		}
		reduce();
	}
	Rational(int p) : p(p), q(1) {}

	friend std::ostream& operator << (std::ostream& out, const Rational& r) {
		if (r.q != 1) {
			out << r.p << '/' << r.q;
		}
		else {
			out << r.p;
		}
		return out;
	}

	friend std::istream& operator >> (std::istream& in, Rational& r) {
		string buffer;
		in >> buffer;
		int length = buffer.length();
		int SlashPlace = length;
		for (int i = 0; i < length; i++) {
			if (buffer[i] == '/') {
				SlashPlace = i;
			}
		}
		if (SlashPlace < length) {
			string PString, QString;
			for (int i = 0; i < SlashPlace; i++) {
				PString.push_back(buffer[i]);
			}
			for (int i = SlashPlace + 1; i < length; i++) {
				QString.push_back(buffer[i]);
			}
			r.p = atoi(PString.c_str());
			r.q = atoi(QString.c_str());
		}
		else {
			r.p = atoi(buffer.c_str());
			r.q = 1;
		}
		r.reduce();
		return in;
	}

	int32_t getNumerator() const {
		return p;
	}

	uint32_t getDenominator() const {
		return q;
	}

	friend Rational operator + (const Rational& first, const Rational& second) {
		Rational res;
		res.p = first.p * second.q + first.q * second.p;
		res.q = first.q * second.q;
		res.reduce();
		return res;
	}

	friend Rational operator * (const Rational& first, const Rational& second) {
		Rational res;
		res.p = first.p * second.p;
		res.q = first.q * second.q;
		res.reduce();
		return res;
	}

	friend Rational operator - (const Rational& first, const Rational& second) {
		Rational res;
		res.p = first.p * second.q - first.q * second.p;
		res.q = first.q * second.q;
		res.reduce();
		return res;
	}


	Rational& operator -= (const Rational& other) {
		*this = *this - other;
		return *this;
	}

	Rational& operator += (const Rational& other) {
		*this = *this + other;
		return *this;
	}

	Rational& operator *= (const Rational& other) {
		*this = *this * other;
		return *this;
	}

	Rational& operator /= (const Rational& other) {
		*this = *this / other;
		return *this;
	}

	Rational& operator ++ () {
		return *this += 1;
	}

	Rational operator ++ (int dummy) {
		Rational res = *this;
		*this += 1;
		return res;
	}

	Rational& operator -- () {
		return *this -= 1;
	}

	Rational operator -- (int dummy) {
		Rational res = *this;
		*this -= 1;
		return res;
	}

	Rational operator - () const {
		Rational res;
		res.p = -(this->p);
		res.q = this->q;
		return res;
	}

	friend Rational operator / (const Rational& first, const Rational& second) {
		if (second == 0) {
			throw RationalDivisionByZero();
		}
		Rational res;
		if (second.p < 0) {
			res.p = -(first.p * (int32_t)second.q);
		}
		else {
			res.p = first.p * second.q;
		}
		res.q = first.q * abs(second.p);
		res.reduce();
		return res;
	}


	friend bool operator < (const Rational& first, const Rational& second) {
		return (first.p * (int32_t)second.q < (int32_t)first.q * second.p);
	}

	friend bool operator == (const Rational& first, const Rational& second) {
		return !(second < first || second > first);
	}

	friend bool operator != (const Rational& first, const Rational& second) {
		return !(second == first);
	}

	friend bool operator > (const Rational& first, const Rational& second) {
		return (second < first);
	}

	friend bool operator >= (const Rational& first, const Rational& second) {
		return !(first < second);
	}

	friend bool operator <= (const Rational& first, const Rational& second) {
		return !(first > second);
	}
};

class MatrixAllocationError {
};

class MatrixWrongSizeError {
};

class MatrixIndexError {
};

class MatrixIsDegenerateError {
};

// non-specified functions to get "zero" and "one" of type T

template<typename T> T getZero() {
	return T(0);
}

template<typename T> T getOne() {
	return T(1);
}

template<typename T>
class Matrix {
private:
	size_t rowsCnt, colsCnt;
	T **data;
public:
	explicit Matrix(const std::vector<std::vector<T>>& vectorMatrix) {
		rowsCnt = vectorMatrix.size();
		colsCnt = vectorMatrix[0].size();
		data = new T*[rowsCnt];
		for (size_t i = 0; i < rowsCnt; ++i) {
			data[i] = new T[colsCnt];
			for (size_t j = 0; j < colsCnt; ++j) {
				data[i][j] = vectorMatrix[i][j];
			}
		}
	}

	Matrix(size_t n, size_t m) {
		this->rowsCnt = n;
		this->colsCnt = m;
		this->data = new T*[n];
		for (size_t i = 0; i < n; ++i) {
			this->data[i] = new T[m];
		}
	}

	Matrix(const Matrix& other) : data(nullptr) {
		*this = other;
	}

	Matrix& operator = (const Matrix& other) {
		if (this->data != other.data) {
			if (this->data != nullptr) {
				this->~Matrix();
			}
			this->rowsCnt = other.rowsCnt;
			this->colsCnt = other.colsCnt;
			this->data = new T*[rowsCnt];
			for (size_t i = 0; i < rowsCnt; ++i) {
				this->data[i] = new T[colsCnt];
				for (size_t j = 0; j < colsCnt; ++j) {
					this->data[i][j] = other.data[i][j];
				}
			}
		}
		return *this;
	}

	~Matrix() {
		for (size_t i = 0; i < rowsCnt; ++i) {
			delete[] data[i];
		}
		delete[] data;
	}

	size_t getRowsNumber() const {
		return rowsCnt;
	}

	size_t getColumnsNumber() const {
		return colsCnt;
	}

	Matrix operator * (T mult) const {
		Matrix res(*this);
		for (size_t i = 0; i < res.getRowsNumber(); ++i) {
			for (size_t j = 0; j < res.getColumnsNumber(); ++j) {
				res.data[i][j] *= mult;
			}
		}
		return res;
	}

	Matrix& operator *= (T mult) {
		*this = *this * mult;
		return *this;
	}

	friend Matrix operator * (T mult, const Matrix& matrix) {
		return matrix * mult;
	}

	Matrix operator * (const Matrix& other) const {
		if (this->colsCnt != other.rowsCnt) {
			throw MatrixWrongSizeError();
		}
		Matrix res(this->rowsCnt, other.colsCnt);
		for (size_t i = 0; i < res.rowsCnt; ++i) {
			for (size_t j = 0; j < res.colsCnt; ++j) {
				res.data[i][j] = getZero<T>();
				for (size_t k = 0; k < colsCnt; ++k) {
					res.data[i][j] += this->data[i][k] * other.data[k][j];
				}
			}
		}
		return res;
	}

	Matrix& operator *= (const Matrix& other) {
		*this = *this * other;
		return *this;
	}

	Matrix& operator += (const Matrix& other) {
		if (this->colsCnt != other.colsCnt || this->rowsCnt != other.rowsCnt) {
			throw MatrixWrongSizeError();
		}
		for (size_t i = 0; i < rowsCnt; ++i) {
			for (size_t j = 0; j < colsCnt; ++j) {
				data[i][j] += other.data[i][j];
			}
		}
		return *this;
	}

	Matrix operator + (const Matrix& other) const {
		Matrix res(*this);
		res += other;
		return res;
	}

	Matrix& operator -= (const Matrix& other) {
		*this += other * (-1);
		return *this;
	}

	Matrix operator - (const Matrix& other) const {
		Matrix res(*this);
		res -= other;
		return res;
	}

	friend std::ostream& operator << (std::ostream& out, const Matrix& matr) {
		for (size_t i = 0; i < matr.rowsCnt; ++i) {
			for (size_t j = 0; j < matr.colsCnt; ++j) {
				out << matr.data[i][j] << ' ';
			}
			out << '\n';
		}
		return out;
	}

	friend std::istream& operator >> (std::istream& in, Matrix& matr) {
		for (size_t i = 0; i < matr.rowsCnt; ++i) {
			for (size_t j = 0; j < matr.colsCnt; ++j) {
				in >> matr.data[i][j];
			}
		}
		return in;
	}

	T operator () (int i, int j) const {
		if (i < 0 || i >= this->getRowsNumber() || j < 0 || j >= this->getColumnsNumber()) {
			throw MatrixIndexError();
		}
		return data[i][j];
	}

	T& operator () (int i, int j) {
		if (i < 0 || i >= this->getRowsNumber() || j < 0 || j >= this->getColumnsNumber()) {
			throw MatrixIndexError();
		}
		return data[i][j];
	}

	Matrix getTransposed() const {
		Matrix res(this->colsCnt, this->rowsCnt);
		for (size_t i = 0; i < res.rowsCnt; ++i) {
			for (size_t j = 0; j < res.colsCnt; ++j) {
				res.data[i][j] = this->data[j][i];
			}
		}
		return res;
	}

	Matrix& transpose() {
		*this = this->getTransposed();
		return *this;
	}

};

template<typename T>
class SquareMatrix : public Matrix<T> {
public:
	SquareMatrix(size_t size) : Matrix<T>::Matrix(size, size) {}

	SquareMatrix(const SquareMatrix& other) : SquareMatrix(other.getSize()){
		*this = other;
	}

	SquareMatrix(const Matrix<T>& other) : SquareMatrix(other.getRowsNumber()) {
		if (other.getColumnsNumber() != other.getRowsNumber()) {
			throw MatrixWrongSizeError();
		}
		Matrix<T>::operator = (other);
	}

	SquareMatrix& operator = (const SquareMatrix& other) {
		Matrix<T>::operator = (other); // ¬ызов метода базового класса
		return *this; // изменилс€ базовый класс, который тут не заметен
	}

	size_t getSize() const {
		return Matrix<T>::getRowsNumber();
	}

	SquareMatrix operator * (T mult) const {
		return Matrix<T>::operator * (mult);
	}

	SquareMatrix& operator *= (T mult) {
		*this = *this * mult;
		return *this;
	}

	friend SquareMatrix operator * (T mult, const SquareMatrix& matrix) {
		return matrix * mult;
	}

	SquareMatrix operator * (const SquareMatrix& other) const {
		return Matrix<T>::operator * (other);
	}

	SquareMatrix& operator *= (const SquareMatrix& other) {
		*this = *this * other;
		return *this;
	}

	Matrix<T> operator * (const Matrix<T>& other) const {
		return Matrix<T>::operator * (other);
	}

	Matrix<T>& operator *= (const Matrix<T>& other) {
		*this = *this * other;
		return *this;
	}

	SquareMatrix& operator += (const SquareMatrix& other) {
		*this = *this + other;
		return *this;
	}

	SquareMatrix operator + (const SquareMatrix& other) const {
		return Matrix<T>::operator + (other);
	}

	SquareMatrix& operator -= (const SquareMatrix& other) {
		*this = *this - other;
		return *this;
	}

	SquareMatrix operator - (const SquareMatrix& other) const {
		return Matrix<T>::operator - (other);
	}

	SquareMatrix getTransposed() const {
		return this->Matrix<T>::getTransposed();
	}

	SquareMatrix& transpose() {
		*this = this->getTransposed();
		return *this;
	}

	T getDeterminant() const {
		/*if (this->getSize() > 1) {
			T res = getZero<T>();
			for (size_t k = 0; k < this->getSize(); ++k) {
				SquareMatrix<T> minor(this->getSize() - 1);
				for (size_t i = 0; i < this->getSize() - 1; ++i) {
					for (size_t j = 0; j < this->getSize() - 1; ++j) {
						minor.operator()(i, j) = this->operator()(i + 1, (j < k) ? j : j + 1);
					}
				}
				if (k % 2 == 0) {
					res += minor.getDeterminant() * this->operator() (0, k);
				}
				else {
					res -= minor.getDeterminant() * this->operator() (0, k);
				}
			}
			return res;
		}
		else {
			return this->operator() (0, 0);
		}*/
		SquareMatrix<T> res(*this);
		int count = 0;
		for (size_t i = 0; i < res.getSize(); ++i) {
			size_t k = i;
			while (k < res.getSize() && res.operator() (k, i) == 0) {
				k++;
			}
			if (k == res.getSize()) {
				return getZero<T>();
			}
			T denum = res.operator() (k, i);
			if (k != i) {
				for (size_t j = 0; j < res.getSize(); ++j) {
					swap(res.operator() (k, j), res.operator() (i, j));
				}
				count++;
			}
			for (size_t j = i + 1; j < res.getSize(); ++j) {
				T num = res.operator() (j, i);
				for (size_t h = 0; h < res.getSize(); ++h) {
					res.operator() (j, h) -= res.operator() (i, h) * num / denum;
				}
			}
		}
		T result = getOne<T>();
		for (size_t i = 0; i < res.getSize(); ++i) {
			result *= res.operator() (i, i);
		}
		if (count % 2 != 0) {
			result = -result;
		}
		return result;
	}

	SquareMatrix& invert() {
		if (this->getDeterminant() == 0) {
			throw MatrixIsDegenerateError();
		}
		SquareMatrix res(this->getSize());
		for (size_t i = 0; i < this->getSize(); ++i) {
			for (size_t j = 0; j < this->getSize(); ++j) {
				res.operator() (i, j) = (i == j) ? getOne<T>() : getZero<T>();
			}
		}
		for (size_t i = 0; i < this->getSize(); ++i) {
			size_t k = i;
			while (this->operator() (k, i) == 0) {
				k++;
			}
			T denum = this->operator() (k, i);
			for (size_t j = 0; j < this->getSize(); ++j) {
				swap(this->operator() (k, j), this->operator() (i, j));
				swap(res.operator() (k, j), res.operator() (i, j));

				this->operator() (i, j) /= denum;
				res.operator() (i, j) /= denum;
			}
			for (size_t j = 0; j < this->getSize(); ++j) {
				if (j != i) {
					T num = this->operator() (j, i);
					for (size_t h = 0; h < this->getSize(); ++h) {
						this->operator() (j, h) -= this->operator() (i, h) * num;
						res.operator() (j, h) -= res.operator() (i, h) * num;
					}
				}
			}
		}
		*this = res;
		return *this;
	}

	SquareMatrix getInverse() const {
		SquareMatrix res(*this);
		return res.invert();
	}

	T getTrace() const {
		T res = getZero<T>();
		for (size_t i = 0; i < this->getSize(); ++i) {
			res += this->operator() (i, i);
		}
		return res;
	}

};

int main() {
	/*Matrix<Rational> a(
	{
		{ 1, 2, 3},
		{ 4, 5, 6},
		{ 7, 8, 9}
	});
	SquareMatrix<Rational> b(a);
	cout << b.getSize();
	cout << b.getDeterminant();

	SquareMatrix<Rational> b(a);
	std::cout << a * Rational(1, 2) << '\n';
	std::cout << b.invert() << '\n';
	std::cout << b;

	SquareMatrix<Rational> c(b);
	b.transpose();
	std::cout << a << b << c << '\n';
	c = b;
	std::cout << a << b << c << '\n';
	b += a;*/

	int m, n, p, q;
	cin >> m >> n >> p >> q;

	Matrix<int> A(m, n), B(p, q);
	cin >> A >> B;

	A = A;
	try {
		cout << A + B * 2 - m * A << endl;
		cout << (A -= B += A *= 2) << endl;
		cout << (((A -= B) += A) *= 2) << endl;
	}
	catch (const MatrixWrongSizeError&) {
		cout << "A and B are of different size." << endl;
	}
	B = A;
	cout << B << endl;

	Rational r;
	cin >> r;
	Matrix<Rational> C(m, n), D(p, q);
	cin >> C >> D;
	try {
		cout << C * D << endl;
		cout << (C *= D) << endl;
		cout << C << endl;
	}
	catch (const MatrixWrongSizeError&) {
		cout << "C and D have not appropriate sizes for multiplication." << endl;
	}
	cout << C.getTransposed() * (r * C) << endl;
	cout << C.transpose() << endl;
	cout << C << endl;

	SquareMatrix<Rational> S(m);
	cin >> S;
	SquareMatrix<Rational> P(S);
	const SquareMatrix<Rational>& rS = S;
	cout << rS.getSize() << ' ' << rS.getDeterminant() << ' ' << rS.getTrace() << endl;
	cout << (S = S) * (S + rS) << endl;
	cout << (S *= S) << endl;
	C.transpose();
	cout << rS * C << endl;
	cout << S << endl;
	S = P;
	cout << (Rational(1, 2) * S).getDeterminant() << endl;
	try {
		cout << rS(0, 0) << endl;
		(S(0, 0) *= 2) /= 2;
		cout << rS(0, m) << endl;
	}
	catch (const MatrixIndexError&) {
		cout << "Index out of range." << endl;
	}
	cout << rS.getTransposed() << endl;
	try {
		cout << rS.getInverse() << endl;
		cout << S.invert().getTransposed().getDeterminant() << endl;
		cout << S << endl;
	}
	catch (const MatrixIsDegenerateError&) {
		cout << "Cannot inverse S." << endl;
	}

	//system("pause");
	return 0;
}