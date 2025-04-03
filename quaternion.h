/**
 * @file quaternion.h
 * @brief 四元数がらみのライブラリ
 * @author csharpython
 * @version v0.0.0
 */
#pragma once
#include <cmath>
#include <type_traits>
#include <stdexcept>
namespace qtnion{
	template <typename T>
	struct quaternion {
		T one;
		T i;
		T j;
		T k;
		static_assert(std::is_arithmetic<T>::value,"quaternion requires arithmetic type.");
		/** 
		 * @brief 四元数を表す構造体。
		 */
		quaternion() {one = i = j = k = 0;}//0
		quaternion(T one_,T i_=0) { one = one_; i=i_; j = k = 0;}//From ℝ and ℂ
		quaternion(T i_,T j_,T k_) {one = 0; i = i_; j=j_; k=k_;}//From 3D vector
		quaternion(T one_,T i_ , T j_,T k_) { one = one_; i=i_; j=j_; k=k_;}//From ℍ
		/**
		 * @brief 四元数に対する通常の加算
		 * @param rhs 加数
		 */
		inline quaternion operator+(const quaternion& rhs) const
			{ return {one+rhs.one,i+rhs.i,j+rhs.j,k+rhs.k};}
		/**
		 * @brief 四元数に対する通常の減算
		 * @param rhs 減数
		 */
		inline quaternion operator-(const quaternion& rhs) const
			{return {one-rhs.one,i-rhs.i,j-rhs.j,k-rhs.k};}
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
		inline quaternion operator/(const T rhs) const { 
			if(rhs == 0) throw std::invalid_argument("Cannot divide by 0");
			return {one/rhs,i/rhs,j/rhs,k/rhs};
		}
		inline bool operator==(const quaternion& rhs) const 
			{ return (one==rhs.one&&i==rhs.i&&j==rhs.j&&k==rhs.k);}
	};
	/**
	 * @brief 四元数に対する共役
	* @return 四元数の共役
	*/
	template <typename T>
	inline quaternion<T> conjugate(quaternion<T> val)
		{ return {val.one,-val.i,-val.j,-val.k};}
	/**
	 * @brief 四元数に対するノルムの二乗を計算
	 * @return 引数のノルムの二乗
	 */
	template <typename T>
	inline T squ_norm(quaternion<T> val)
		{return val.one*val.one+val.i*val.i+val.j*val.j+val.k*val.k;}
	/**
	 * @brief 四元数の逆数を計算
	 * @return 引数の逆数
	 */
	template <typename T>
	inline quaternion<T> inverse(quaternion<T> val) {return conjugate(val)/squ_norm(val);}
	/**
	 * @brief 四元数に対するノルムを計算
	 * @return 自身のノルム
	 * @remark 二乗を求めたい場合、squ_normの方がよいと思われます。
	 * @see squ_norm()
	 */
	template <typename T>
	inline T norm(const quaternion<T> arg){return sqrt(squ_norm(arg));}
	/**
	 * @brief 四元数の標準化をする
	 * @return もともとの値をそれのノルムで割ったもの
	 */
	template <typename T>
	inline quaternion<T> normalize(quaternion<T> val) {return val/norm(val);}
	/**
	 * @brief 四元数を3次元空間での回転から作る
	 * @attention 与えられる空間ベクトルは正規化済みであることを前提とする。
	 * @return 回転軸と回転角度から求められる四元数
	 */
	template <typename T>
	inline quaternion<T> polarturn(T x,T y,T z,T theta){
		return {cos(theta/2),x*sin(theta/2),y*sin(theta/2),z*sin(theta/2)};}
	/**
	 * @brief 空間ベクトルを四元数に基づいて回転させる
	 * @return 空間ベクトルを回転させたものを表す四元数
	 */
	template <typename T>
	inline quaternion<T> turn3Dvec(quaternion<T> v,quaternion<T> q){return q*v*conjugate(q);}
}