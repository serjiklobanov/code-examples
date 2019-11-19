#include <iostream>
#include <math.h>
#include <algorithm>

class GeometryObject;
class Point;
class Segment;
class Line;
class Ray;
class Polygon;

struct DoublePoint {
	double x, y;
};

int sign(long long int i) {
	int sign;
	if (i > 0) {
		sign = 1;
	}
	if (i < 0) {
		sign = -1;
	}
	if (i == 0) {
		sign = 0;
	}
	return sign;
}

class Vector {
private:
	int x, y;
public:
	friend Point;
	friend Line;
	friend Ray;
	Vector(int x, int y) : x(x), y(y) {}
	Vector(const Point& A, const Point& B);
	Vector operator + (const Vector& other) const {
		Vector res(this->x + other.x, this->y + other.y);
		return res;
	}
	Vector& operator += (const Vector& other) {
		*this = *this + other;
		return *this;
	}
	Vector operator - () const {
		Vector res(-this->x, -this->y);
		return res;
	}
	Vector operator - (const Vector& other) const {
		Vector res = *this + (-other);
		return res;
	}
	double length() const {
		return sqrt(x * x + y * y);
	}
	long long int operator * (const  Vector& other) const {
		return ((long long)this->x * other.y - (long long)this->y * other.x);
	}
	long long int operator ^ (const  Vector& other) const {
		return ((long long)this->x * other.x + (long long)this->y * other.y);
	}
	friend std::ostream& operator << (std::ostream& out, const Vector& vec) {
		out << vec.x << ' ' << vec.y;
		return out;
	}
};

class GeometryObject {
public:
	virtual void shift(const Vector &) = 0;
	virtual bool contains_point(const Point &) const = 0;
	virtual bool cross_segment(const Segment &) const = 0;
};

class Point : public GeometryObject {
private:
	int x, y;
public:
	friend Vector;
	friend Segment;
	friend Line;
	friend Ray;
	Point(int x, int y) : x(x), y(y) {};
	Point() {};
	bool contains_point(const Point &point) const override {
		if (this->x == point.x && this->y == point.y) {
			return true;
		}
		else {
			return false;
		}
	}
	double point_distance(const Point& point) const {
		Vector dist(*this, point);
		return dist.length();
	}
	void shift(const Vector &vector) override {
		this->x += vector.x;
		this->y += vector.y;
	}
	bool cross_segment(const Segment &segment) const override;
	friend std::istream& operator >> (std::istream& in, Point& point) {
		in >> point.x >> point.y;
		return in;
	}
	friend std::ostream& operator << (std::ostream& out, const Point& point) {
		out << point.x << ' ' << point.y;
		return out;
	}
	friend int CompareX(const Point& point1, const Point& point2) {
		if (point1.x < point2.x) {
			return -1;
		}
		if (point1.x == point2.x) {
			return 0;
		}
		if (point1.x > point2.x) {
			return 1;
		}
	}
	friend int CompareY(const Point& point1, const Point& point2) {
		if (point1.y < point2.y) {
			return -1;
		}
		if (point1.y == point2.y) {
			return 0;
		}
		if (point1.y > point2.y) {
			return 1;
		}
	}
	friend int ComparePolarAngle(const Point& point1, const Point& point2, const Point& point) {
		Vector vector1(point, point1), vector2(point, point2);
		if (vector1 * vector2 < 0) {
			return 1;
		}
		if (vector1 * vector2 == 0) {
			return 0;
		}
		if (vector1 * vector2 > 0) {
			return -1;
		}
	}
};

class Segment : public GeometryObject {
private:
	Point A, B;
public:
	friend Line;
	friend Ray;
	Segment(Point A, Point B) : A(A), B(B) {};
	bool contains_point(const Point &point) const override {
		Vector AB(A, B), AP(A, point), BP(B, point);
		return (AB * AP == 0 && (AB ^ AP) >= 0 && (AB ^ BP) <= 0);
	}
	double point_distance(const Point& point) const;
	void shift(const Vector &vector) override {
		A.shift(vector);
		B.shift(vector);
	}
	bool cross_segment(const Segment &segment) const override {
		Vector A1B1(A, B), A2B2(segment.A, segment.B), A1A2(A, segment.A), A1B2(A, segment.B), A2A1(segment.A, A), A2B1(segment.A, B);
		return (sign(A1B1 * A1A2) * sign(A1B1 * A1B2) <= 0 &&
			sign(A2B2 * A2A1) * sign(A2B2 * A2B1) <= 0 &&
			std::min(segment.A.x, segment.B.x) <= std::max(A.x, B.x) &&
			std::min(A.x, B.x) <= std::max(segment.A.x, segment.B.x) &&
			std::min(segment.A.y, segment.B.y) <= std::max(A.y, B.y) &&
			std::min(A.y, B.y) <= std::max(segment.A.y, segment.B.y));
	}
	double segment_distance(const Segment &segment) const {
		if (cross_segment(segment)) {
			return 0;
		}
		else {
			return std::min(std::min(point_distance(segment.A), point_distance(segment.B)), std::min(segment.point_distance(A), segment.point_distance(B)));
		}
	}
};

class Line : public GeometryObject {
private:
	int A, B, C;
public:
	Line(int A, int B, int C) : A(A), B(B), C(C) {}
	Line(const Point& a, const Point& b) {
		A = a.y - b.y;
		B = b.x - a.x;
		C = -A * a.x - B * a.y;
	}
	Line(const Ray& ray);
	Line(const Segment& segment) {
		Line line(segment.A, segment.B);
	}
	long long int substitute_point(const Point& point) const {
		return (long long)A * point.x + (long long)B * point.y + (long long)C;
	}
	bool contains_point(const Point &point) const override {
		if (sign(this->substitute_point(point)) == 0) {
			return true;
		}
		else {
			return false;
		}
	}
	void shift(const Vector &vector) override {
		C -= (A * vector.x + B * vector.y);
	}
	bool cross_segment(const Segment &segment) const override {
		if (sign(this->substitute_point(segment.A)) * sign(this->substitute_point(segment.B)) <= 0) {
			return true;
		}
		else {
			return false;
		}
	}
	double point_distance(const Point& point) const {
		return (double)(abs(this->substitute_point(point)) / sqrt(A * A + B * B));
	}
	friend DoublePoint line_intersection(const Line& line1, const Line& line2) {
		double det = line1.A * line2.B - line1.B * line2.A;
		double x = (-line1.C * line2.B + line1.B * line2.C) / det;
		double y = (-line1.A * line2.C + line1.C * line2.A) / det;
		DoublePoint res;
		res.x = x;
		res.y = y;
		return res;
	}
	double line_distance(const Line& line) const {
		double k;
		if (this->A != 0) {
			k = line.A / this->A;
		}
		else {
			k = line.B / this->B;
		}
		return abs(line.C - this->C * k) / sqrt(line.A * line.A + line.B * line.B);
	}
	bool paralel_line(const Line& line) const {
		Vector normal1(A, B), normal2(line.A, line.B);
		if (normal1 * normal2 == 0) {
			return true;
		}
		else {
			return false;
		}
	}
};

class Ray : public GeometryObject {
private:
	Point start;
	Vector direct;
public:
	friend Line;
	Ray(const Point& point, const Vector& vec) : start(point), direct(vec) {}
	bool contains_point(const Point &point) const override {
		Vector vec(start, point);
		if (vec * direct == 0 && (vec ^ direct) >= 0) {
			return true;
		}
		else {
			return false;
		}
	}
	double point_distance(const Point& point) const {
		Vector vectorCA(point, start);
		Line lineAB(*this);
		if (sign(vectorCA ^ direct) >= 0) {
			return start.point_distance(point);
		}
		else {
			return lineAB.point_distance(point);
		}
	}
	void shift(const Vector &vector) override {
		start.shift(vector);
	}

	bool cross_segment(const Segment &segment) const override {
		Vector AB(segment.A, segment.B), SA(start, segment.A), SB(start, segment.B);
		if (sign(direct * SA) * sign(direct *SB) <= 0
			&& sign(AB * SA) * sign(AB * direct) >= 0
			&& sign(AB * SB) * sign(AB * direct) >= 0
			&& ((SA ^ direct) >= 0 || (SB ^ direct) >= 0)) {
			return true;
		}
		else {
			return false;
		}
	}
};

class Polygon : public GeometryObject {
private:
	int size;
	Point *points;
public:
	Polygon(int size) :size(size) {
		points = new Point[size];
	};
	Polygon(int size, const Point* ptr) :size(size) {
		points = new Point[size];
		for (int i = 0; i < size; i++) {
			points[i] = ptr[i];
		}
	};
	Polygon(const Polygon& other) {
		Polygon(other.size);
		*this = other;
	}
	~Polygon() {
		delete[] points;
	}
	Polygon& operator = (const Polygon& other) {
		this->~Polygon();
		this->size = other.size;
		this->points = new Point[size];
		for (int i = 0; i < size; i++) {
			this->points[i] = other.points[i];
		}
		return *this;
	}

	friend std::istream& operator >> (std::istream& in, Polygon& polygon) {
		Point InputPoint;
		for (int i = 0; i < polygon.size; i++) {
			in >> InputPoint;
			polygon.points[i] = InputPoint;
		}
		return in;
	}

	friend std::ostream& operator << (std::ostream& out, const Polygon& polygon) {
		for (int i = 0; i < polygon.size; i++) {
			out << polygon.points[i] << '\n';
		}
		return out;
	}

	void shift(const Vector & vector) override {
		for (int i = 0; i < size; i++) {
			points[i].shift(vector);
		}
	}
	bool contains_point(const Point & point) const override {
		Point LastPoint = points[size - 1];
		for (int i = 0; i < size; i++) {
			Segment segment(LastPoint, points[i]);
			if (segment.contains_point(point)) {
				return true;
			}
			LastPoint = points[i];
		}
		Vector vector(1, 0);
		Ray ray(point, vector);
		int CrossPoint = 0, CrossSegment = 0, ExeptSign = 0;
		for (int i = 0; i < size; i++) {
			if (ray.contains_point(points[i])) {
				Vector PrevVector(point, points[(i - 1 + size) % size]), NextVector(point, points[(i + 1) % size]);
				if (sign(vector * PrevVector) * sign(vector * NextVector) > 0) {
					CrossPoint--;
				}
				if (sign(vector * NextVector) == 0) {
					ExeptSign = sign(vector * PrevVector);
				}
				if (sign(vector * PrevVector) == 0) {
					if (ExeptSign * sign(vector * NextVector) > 0) {
						CrossPoint--;
					}
				}
				CrossPoint++;
			}
		}
		LastPoint = points[size - 1];
		for (int i = 0; i < size; i++) {
			Segment segment(LastPoint, points[i]);
			if (ray.cross_segment(segment)) {
				CrossSegment++;
			}
			LastPoint = points[i];
		}
		if ((CrossSegment - CrossPoint) % 2 == 0) {
			return false;
		}
		else {
			return true;
		}
	}

	bool is_convex() const {
		int CountPositive = 0, CountNegative = 0;
		for (int i = 0; i < size; i++) {
			Vector LastVector(points[(i - 1 + size) % size], points[i]), NewVector(points[i], points[(i + 1) % size]);
			if (sign(LastVector * NewVector) >= 0) {
				CountPositive++;
			}
			if (sign(LastVector * NewVector) <= 0) {
				CountNegative++;
			}
		}
		if (CountNegative == size || CountPositive == size) {
			return true;
		}
		else {
			return false;
		}
	}

	double square() const {
		long long square = 0;
		for (int i = 0; i < size; i++) {
			Vector LastVector(points[0], points[i]), NewVector(points[0], points[(i + 1) % size]);
			square += (LastVector * NewVector);
		}
		return abs(square) * 0.5;
	}

	bool cross_segment(const Segment& segment) const override {
		Point LastPoint = points[size - 1];
		for (int i = 0; i < size; i++) {
			Segment NewSegment(LastPoint, points[i]);
			if (NewSegment.cross_segment(segment)) {
				return true;
			}
			LastPoint = points[i];
		}
		return false;
	}
};

bool Point::cross_segment(const Segment& segment) const {
	return segment.contains_point(*this);
}

double Segment::point_distance(const Point& point) const {
	Line lineAB(A, B);
	Vector vectorAB(A, B);
	Vector vectorCA(point, A);
	Vector vectorCB(point, B);
	if (sign(vectorAB ^ vectorCA) * sign(vectorAB ^ vectorCB) <= 0) {
		return lineAB.point_distance(point);
	}
	else {
		return std::min(vectorCA.length(), vectorCB.length());
	}
}

Vector::Vector(const Point& A, const Point& B) {
	x = B.x - A.x;
	y = B.y - A.y;
}

Line::Line(const Ray& ray) {
	Point A = ray.start, B = ray.start;
	A.shift(ray.direct);
	Line line(A, B);
	*this = line;
}

int main() {
	return 0;
}
