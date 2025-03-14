#include <iostream>
#include <cmath>
#include <algorithm>

class Vector3D {
public:
    double x, y, z;

    Vector3D(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {
    }

    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3D cross(const Vector3D& other) const {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    double length() const {
        return std::sqrt(dot(*this));
    }

    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }
};

class Point3D {
public:
    double x, y, z;

    Point3D(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {
    }

    Vector3D operator-(const Point3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    Point3D operator+(const Vector3D& vec) const {
        return Point3D(x + vec.x, y + vec.y, z + vec.z);
    }
};

class Segment3D {
public:
    Point3D start, end;

    Segment3D(const Point3D& start_, const Point3D& end_)
        : start(start_), end(end_) {
    }

    Vector3D direction() const {
        return end - start;
    }

    double length() const {
        return direction().length();
    }
};

double distanceBetweenSegments(const Segment3D& seg1, const Segment3D& seg2) {
    Vector3D u = seg1.direction();
    Vector3D v = seg2.direction();
    Vector3D w = seg1.start - seg2.start;

    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    double denom = a * c - b * b;

    if (std::abs(denom) < 1e-10) {
        double t = (d < 0) ? 0.0 : (d > a ? 1.0 : d / a);
        Vector3D closestOnSeg1 = u * t;
        Vector3D diff = w + closestOnSeg1;
        double dist = diff.length();

        double s = v.dot(diff) / c;
        s = std::clamp(s, 0.0, 1.0);
        Vector3D closestOnSeg2 = v * s;
        diff = w + closestOnSeg1 - closestOnSeg2;
        return diff.length();
    }

    double s = (b * e - c * d) / denom;
    double t = (a * e - b * d) / denom;

    s = std::clamp(s, 0.0, 1.0);
    t = std::clamp(t, 0.0, 1.0);

    Point3D closest1 = seg1.start + (u * s);
    Point3D closest2 = seg2.start + (v * t);

    Vector3D diff = closest1 - closest2;
    return diff.length();
}

int main() {
    Segment3D seg1(Point3D(0, 0, 0), Point3D(1, 1, 1));
    Segment3D seg2(Point3D(0, 1, 0), Point3D(1, 2, 0));

    double dist = distanceBetweenSegments(seg1, seg2);
    std::cout << "Distance between segments: " << dist << std::endl;

    return 0;
}