/**
 * @file quaternion.h
 * @brief 四元数がらみのライブラリ
 * @author csharpython
 * @version v0.0.0
 */
#pragma once
#include <cmath>
namespace qtnion{
	template <typename T>
	struct quaternion {
		T one;
		T i;
		T j;
		T k;
		/** 
		 * @brief 四元数を表す構造体。
		 */
		quaternion() { one = i = j = k = 0; }
		quaternion(T one_,T i_ , T j_,T k_) {
			one = one_;
			i=i_;
			j=j_;
			k=k_;
		}
		/**
		 * @brief 四元数に対する通常の加算
		 * @param rhs 加数
		 */
		quaternion operator+(const quaternion& rhs) const {
			quaternion ret;
			ret.one = one + rhs.one;
			ret.i = i + rhs.i;
			ret.j = j + rhs.j;
			ret.k = k + rhs.k;
			return ret;
		}
		/**
		 * @brief 四元数に対する通常の減算
		 * @param rhs 減数
		 */
		quaternion operator-(const quaternion& rhs) const {
			quaternion ret;
			ret.one = one - rhs.one;
			ret.i = i - rhs.i;
			ret.j = j - rhs.j;
			ret.k = k - rhs.k;
			return ret;
		}
		/**
		 * @brief 四元数に対する通常の乗算
		 * @remarks 非可換です。つまり、左辺と右辺を入れ替えると結果が変わ(ることがあ)ります。
		 * @param rhs 乗数
		 */
		quaternion operator*(const quaternion& rhs) const {
			quaternion ret;
			ret.one = one * rhs.one - i * rhs.i - j * rhs.j - k*rhs.k;
			ret.i = one * rhs.i + i * rhs.one + j * rhs.k - k * rhs.j;
			ret.j = one * rhs.j - i * rhs.k + j * rhs.one + k * rhs.i;
			ret.k = one * rhs.k + i * rhs.j - j * rhs.i + k * rhs.one;
			return ret;
		}
		/**
		 * @brief 四元数の実数による除算
		 * @param rhs 除数
		 */
		quaternion operator/(const T rhs) const {
			quaternion ret;
			ret.one = one/rhs;
			ret.i = i/rhs;
			ret.j = j/rhs;
			ret.k = k/rhs;
			return ret;
		}
	};
	/**
	 * @brief 四元数に対する共役
	* @return 四元数の共役
	*/
	template <typename T>
	quaternion<T> conjugate(quaternion<T> val) {
		quaternion<T> ret;
		ret.one=val.one;
		ret.i=-val.i;
		ret.j=-val.j;
		ret.k=-val.k;
		return ret;
	}
	/**
	 * @brief 四元数に対するノルムの二乗を計算
	 * @return 引数のノルムの二乗
	 */
	template <typename T>
	
	T squ_norm(quaternion<T> val) {return val.one*val.one+val.i*val.i+val.j*val.j+val.k*val.k;}
	/**
	 * @brief 四元数の逆数を計算
	 * @return 引数の逆数
	 */
	template <typename T>
	quaternion<T> inverse(quaternion<T> val) {return conjugate(val)/squ_norm(val);}
	/**
	 * @brief 四元数に対するノルムを計算
	 * @return 自身のノルム
	 * @remark 二乗を求めたい場合、squ_normの方がよいと思われます。
	 * @see squ_norm()
	 */
	template <typename T>
	T norm(const quaternion<T> arg){return sqrt(squ_norm(arg));}
	/**
	 * @brief 四元数の標準化をする
	 * @return もともとの値/それのノルム
	 */
	template <typename T>
	quaternion<T> normalize(quaternion<T> val) {return val/norm(val);}
}