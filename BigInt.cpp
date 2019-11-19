#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>

class BigInt {
private:
	constexpr static size_t MAX_LENGTH = 2e4;
	constexpr static size_t MAX_SIZE = 3000;
	constexpr static int BASE = 1000000000;
	constexpr static int RANK = 9;
	constexpr static int MagicDevisionConst = 10;

	typedef uint32_t digit_type;
	typedef uint64_t long_digit_type;

	digit_type data[MAX_SIZE];
	size_t size;
	bool is_negative;

public:
	BigInt() : size(1), is_negative(false) {
		for (size_t i = 0; i < MAX_SIZE; i++) {
			data[i] = 0;
		}
	} // приравнять число 0

	BigInt(int Int) {
		if (Int == 0) {
			*this = BigInt();
		}
		else {
			if (Int < 0) {
				is_negative = true;
				Int = -Int;
			}
			else {
				is_negative = false;
			}
			size = 0;
			while (Int > 0) {
				data[size++] = Int % BASE;
				Int /= BASE;
			}
			for (size_t i = size; i < MAX_SIZE; i++) {
				data[i] = 0;
			}
		}
	}// присвоить int

	BigInt(const char* Str) {
		size_t shift = 0;
		if (Str[0] == '-') {
			is_negative = true;
			shift = 1;
		}
		else {
			is_negative = false;
		}
		size_t length = strlen(Str);
		size = (length - shift - 1) / RANK + 1;
		for (size_t i = 0; i < MAX_SIZE; i++) {
			data[i] = 0;
		}
		for (size_t i = 0; i < length - shift; i++) {
			data[i / RANK] += (digit_type)(Str[length - 1 - i] - '0') * (digit_type)pow(10, i % RANK);
		}
	} // обратите внимание на отрицательные числа! Проверку на корректность делать не надо

	BigInt operator + (const BigInt& other) const {
		if (this->is_negative == other.is_negative) {
			BigInt res;
			res.size = 0;
			long_digit_type remainder = 0;
			while (res.size < std::max(this->size, other.size)) {
				res.data[res.size] = ((long_digit_type)this->data[res.size] + (long_digit_type)other.data[res.size] + remainder) % BASE;
				remainder = ((long_digit_type)this->data[res.size] + (long_digit_type)other.data[res.size] + remainder) / BASE;
				res.size++;
			}
			if (remainder != 0) {
				res.data[res.size] = remainder;
				res.size++;
			}
			res.is_negative = this->is_negative;
			return res;
		}
		else {
			if (this->is_negative) {
				return (other - this->abs());
			}
			else {
				return (*this - other.abs());
			}
		}
	}

	BigInt operator * (const BigInt& other) const {
		BigInt res;
		for (size_t i = 0; i < other.size; i++) {
			BigInt mid;
			mid.size = i;
			long_digit_type remainder = 0;
			for (size_t j = 0; j < this->size; j++) {
				long_digit_type temp = (long_digit_type)this->data[j] * (long_digit_type)other.data[i] + remainder;
				mid.data[i + j] = temp % BASE;
				remainder = temp / BASE;
				mid.size++;
			}
			if (remainder != 0) {
				mid.data[mid.size] = remainder;
				mid.size++;
			}
			res += mid;
		}
		res.is_negative = this->is_negative ^ other.is_negative;
		return res;
	}

	BigInt operator - (const BigInt& other) const {
		if (this->is_negative == other.is_negative) {
			BigInt res;
			res.size = 1;
			int flag = 1;
			if (this->abs() < other.abs()) {
				flag = -1;
			}
			int64_t remainder = 0;
			for (size_t i = 0; i < std::max(this->size, other.size); i++) {
				int temp = 0;
				if ((flag * ((int64_t)this->data[i] - (int64_t)other.data[i]) + remainder) < 0) {
					temp = 1;
				}
				res.data[i] = flag * ((int64_t)this->data[i] - (int64_t)other.data[i]) + remainder + BASE * temp;
				remainder = -temp;
				if (res.data[i] != 0) {
					res.size = i + 1;
				}
			}
			if (this->is_negative) {
				if (flag == 1) {
					res.is_negative = true;
				}
				else {
					res.is_negative = false;
				}
			}
			else {
				if (flag == 1) {
					res.is_negative = false;
				}
				else {
					res.is_negative = true;
				}
			}
			return res;
		}
		else {
			if (this->is_negative) {
				return -(other + this->abs());
			}
			else {
				return (*this + other.abs());
			}
		}
	}

	BigInt operator / (int Int) const {
		BigInt res = this->abs();
		digit_type divisor = std::abs(Int);
		digit_type remainder = 0;
		for (int i = res.size - 1; i >= 0; i--) {
			long_digit_type currentDigit = (long_digit_type)res.data[i] + (long_digit_type)remainder * BASE;
			res.data[i] = currentDigit / divisor;
			remainder = currentDigit % divisor;
		}
		if (this->is_negative) {
			res++;
		}
		if (this->is_negative == (Int < 0)) {
			res.is_negative = false;
		}
		else {
			res.is_negative = true;
		}
		while (res.size > 1 && res.data[res.size - 1] == 0) {
			res.size--;
		}
		return res;
	}

	BigInt operator /= (int Int) {
		return *this = (*this / Int);
	}

	BigInt operator % (int Int) const {
		return (*this - (*this / Int) * Int);
	}

	BigInt operator %= (int Int) {
		return *this = (*this % Int);
	}

	BigInt& operator /= (const BigInt& other) {
		BigInt divisor = other;
		BigInt res = 0, shift = 1;
		bool negative = true;
		if (this->is_negative == other.is_negative) {
			negative = false;
		}
		this->is_negative = false;
		divisor.is_negative = false;
		while (divisor.size < this->size) {
			divisor *= BASE;
			shift *= BASE;
		}
		while (divisor < *this) {
			divisor *= MagicDevisionConst;
			shift *= MagicDevisionConst;
		}
		while (divisor > 0) {
			while (*this >= divisor) {
				*this -= divisor;
				res += shift;
			}
			divisor /= MagicDevisionConst;
			shift /= MagicDevisionConst;
		}
		res.is_negative = negative;
		*this = res;
		return *this;
	}

	BigInt operator / (const BigInt& other) const {
		BigInt res = *this;
		res /= other;
		return res;
	}

	BigInt operator % (const BigInt& other) const {
		return (*this - (*this / other) * other);
	}

	BigInt& operator -= (const BigInt& other) {
		*this = *this - other;
		return *this;
	}

	BigInt& operator += (const BigInt& other) {
		*this = *this + other;
		return *this;
	}

	BigInt& operator *= (const BigInt& other) {
		*this = *this * other;
		return *this;
	}

	BigInt& operator ++ () {
		return *this += 1;
	} // Преинкремент. Пишется: ++a. По сути это += 1

	BigInt operator ++ (int dummy) {
		BigInt res = *this;
		*this += 1;
		return res;
	}// Постинкремент. Пишется: a++. По сути это возвращение старого объекта и изменение текущего. (как с int). //Переменная dummy не используется и не передается. Так вышло по историческим причинам, чтобы отличать пре- от пост- инкремента.

	BigInt& operator -- () {
		return *this -= 1;
	}// Предекремент. Пишется: --a. По сути это -= 1

	BigInt operator -- (int dummy) {
		BigInt res = *this;
		*this -= 1;
		return res;
	} // Постдекремент.

	BigInt operator - () const {
		BigInt res = *this;
		if (res != 0) {
			res.is_negative = !res.is_negative;
		}
		return res;
	}// унарный минус пишется: -a

	bool operator < (const BigInt& other) const {
		if (this->is_negative == other.is_negative) {
			bool temp = false;
			if (this->size < other.size) {
				temp = true;
			}
			if (this->size > other.size) {
				temp = false;
			}
			if (this->size == other.size) {
				bool stop = false;
				for (int i = this->size - 1; i >= 0 && !stop; i--) {
					if (this->data[i] < other.data[i]) {
						temp = true;
						stop = true;
					}
					if (this->data[i] > other.data[i]) {
						temp = false;
						stop = true;
					}
				}
				if (!stop) {
					return false;
				}
			}
			if (this->is_negative) {
				return !temp;
			}
			else {
				return temp;
			}
		}
		else {
			if (this->is_negative) {
				return true;
			}
			else {
				return false;
			}
		}
	}

	bool operator == (const BigInt& other) const {
		return !(other < *this || other > *this);
	}

	bool operator != (const BigInt& other) const {
		return !(other == *this);
	}

	bool operator > (const BigInt& other) const {
		return (other < *this);
	}

	bool operator >= (const BigInt& other) const {
		return !(*this < other);
	}

	bool operator <= (const BigInt& other) const {
		return !(*this > other);
	}

	friend std::ostream& operator << (std::ostream& out, const BigInt& bigint) {
		if (bigint.is_negative) {
			out << '-';
		}
		out << bigint.data[bigint.size - 1];
		for (int i = bigint.size - 2; i >= 0; i--) {
			int count = 0;
			int num = bigint.data[i] + 1;
			while (num > 0) {
				num /= 10;
				count++;
			}
			for (int j = 0; j < RANK - count; j++) {
				out << '0';
			}
			out << bigint.data[i];
		}
		return out;
	}

	friend std::istream& operator >> (std::istream& in, BigInt& bigint) {
		char buffer[MAX_LENGTH];
		in >> buffer;
		bigint = BigInt(buffer);
		return in;
	}

	BigInt abs() const {
		BigInt res = *this;
		res.is_negative = false;
		return res;
	}

	friend BigInt operator + (int a, const BigInt& b) {
		return b + a;
	}

};

BigInt gcd(BigInt first, BigInt second) {
	while (first != 0 && second != 0) {
		if (first > second) {
			first = first % second;
		}
		else {
			second = second % first;
		}
	}
	return (first + second);
}

int main() {
	return 0;
}