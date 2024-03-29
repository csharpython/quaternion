#include <cmath>
namespace qtnion{
    template <typename TYPE>
	struct quaternion {
		TYPE one;
		TYPE i;
		TYPE j;
		TYPE k;
		quaternion() {
			one = 0;
			i = 0;
			j = 0;
			k = 0;
		}
		//Copy constructor omitted
		quaternion operator+(const quaternion& other) {
			quaternion ret;
			ret.one = one + other.one;
			ret.i = i + other.i;
			ret.j = j + other.j;
			ret.k = k + other.k;
			return ret;
		}
		quaternion operator-(const quaternion& other) {
			quaternion ret;
			ret.one = one - other.one;
			ret.i = i - other.i;
			ret.j = j - other.j;
			ret.k = k - other.k;
			return ret;
		}
		quaternion operator*(const quaternion& other){
			quaternion ret;
			ret.one = one * other.one - i * other.i - j * other.j - k*other.k;
			ret.i = one * other.i + i * other.one + j * other.k - k * other.j;
			ret.j = one * other.j - i * other.k + j * other.one + k * other.i;
			ret.k = one * other.k + i * other.j - j * other.i - k * other.one;
			return ret;
		}
	};
	template <typename TYPE>
	quaternion<TYPE> conjugate(const quaternion<TYPE> target){
		quaternion<TYPE> buf;
		buf.i*=-1;
		buf.j*=-1;
		buf.k*=-1;
		return buf;
	}
	template <typename TYPE>
	TYPE norm(const quaternion<TYPE> target){return sqrt(pow(target.one,2)+pow(target.i,2)+pow(target.j,2)+pow(target.k,2));}
	template <typename TYPE>
	quaternion<TYPE> inverse(const quaternion<TYPE> target){
		quaternion<TYPE> buf=qtnion::conjugate(target);
		const TYPE TARGET_NORM=pow(qtnion::norm(target),2);
		buf.one/=TARGET_NORM;
		buf.i/=TARGET_NORM;
		buf.j/=TARGET_NORM;
		buf.k/=TARGET_NORM;
		return buf;
	}
}
