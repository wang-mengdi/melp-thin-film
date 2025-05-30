//////////////////////////////////////////////////////////////////////////
// Common header
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Common_h__
#define __Common_h__
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <vector>
#include <list>
#include <queue>
#include <array>
#include <memory>
#include <iostream>
#include <cmath>

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <filesystem>



////Eigen vector math type alias
using Quaternionf=Eigen::Quaternionf;
using Quaterniond=Eigen::Quaterniond;
using Rotation2f=Eigen::Rotation2Df;
using Rotation2d=Eigen::Rotation2Dd;
using Transform2f=Eigen::Transform<float,2,Eigen::Affine>;
using Transform3f=Eigen::Transform<float,3,Eigen::Affine>;
using Transform2d=Eigen::Transform<double,2,Eigen::Affine>;
using Transform3d=Eigen::Transform<double,3,Eigen::Affine>;

#define Declare_Eigen_Types(my_type,t)  \
using real=my_type;                     \
using Vector1=Eigen::Matrix<real,1,1>;  \
using Vector2=Eigen::Vector2##t;        \
using Vector3=Eigen::Vector3##t;        \
using Vector4=Eigen::Vector4##t;        \
using Vector5=Eigen::Matrix<real,5,1>;; \
using Vector6=Eigen::Matrix<real,6,1>;; \
using Vector7=Eigen::Matrix<real,7,1>;; \
using Vector8=Eigen::Matrix<real,8,1>;; \
using Vector9=Eigen::Matrix<real,9,1>;; \
using VectorX=Eigen::VectorX##t;        \
using Matrix2=Eigen::Matrix2##t;        \
using Matrix3=Eigen::Matrix3##t;        \
using Matrix4=Eigen::Matrix4##t;        \
using MatrixX=Eigen::MatrixX##t;        \
using C=std::complex<real>;             \
using Vector1c=Eigen::Matrix<C,1,1>;    \
using VectorXc=Eigen::VectorXc##t;      \
using Vector2c=Eigen::Vector2c##t;      \
using Vector3c=Eigen::Vector3c##t;      \
using Vector4c=Eigen::Vector4c##t;      \
using Vector5c=Eigen::Matrix<C,5,1>;;	\
using Vector6c=Eigen::Matrix<C,6,1>;;	\
using Quaternion=Eigen::Quaternion##t;  \
using AngleAxis=Eigen::AngleAxis##t;	\
using Rotation2=Eigen::Rotation2D##t;	\
using Transform2=Transform2##t;	\
using Transform3=Transform3##t; \
using ArrayX=Eigen::ArrayX##t;

#define Declare_Eigen_Vector_Types(type,t)		\
using Vector1##t=Eigen::Matrix<type,1,1>;       \
using Vector2##t=Eigen::Vector2##t;             \
using Vector3##t=Eigen::Vector3##t;             \
using Vector4##t=Eigen::Vector4##t;             \
using VectorX##t=Eigen::VectorX##t;				\
using Vector5##t=Eigen::Matrix<type,5,1>;;		\
using Vector6##t=Eigen::Matrix<type,6,1>;;		\
using Vector7##t=Eigen::Matrix<type,7,1>;;		\
using Vector8##t=Eigen::Matrix<type,8,1>;;		\
using Vector9##t=Eigen::Matrix<type,9,1>;;		

#define Declare_Eigen_Matrix_Types(type,t)		\
using Matrix1##t=Eigen::Matrix<type,1,1>;       \
using Matrix2##t=Eigen::Matrix2##t;             \
using Matrix3##t=Eigen::Matrix3##t;             \
using Matrix4##t=Eigen::Matrix4##t;             \
using MatrixX##t=Eigen::MatrixX##t;             

#ifdef USE_FLOAT
Declare_Eigen_Types(float,f)
#else
Declare_Eigen_Types(double,d)
#endif
Declare_Eigen_Vector_Types(int,i)
Declare_Eigen_Vector_Types(float,f)
Declare_Eigen_Vector_Types(double,d)
Declare_Eigen_Matrix_Types(int,i)
Declare_Eigen_Matrix_Types(float,f)
Declare_Eigen_Matrix_Types(double,d)

using TI=int;
using uchar=unsigned char;
using ushort=unsigned short;
template<class T,int d> using Vector=Eigen::Matrix<T,d,1>;
template<class T,int d> using Matrix=Eigen::Matrix<T,d,d>;

////Container alias
////Array
template<class T> using Array=std::vector<T>;
template<class T> using ArrayPtr=std::shared_ptr<Array<T> >;
using size_type=Array<int>::size_type;
////Array with fixed size
template<class T,int n> using ArrayF=std::array<T,n>;
template<class T, int n, int m> using Array2DF = std::array<std::array<T, n>,m>;
constexpr int Pow(int x,int p){return p==1?x:x*Pow(x,p-1);}
constexpr int Factorial(int n){return n<=1?1:(n*Factorial(n-1));}
template<class T,int d> using ArrayF2P=ArrayF<T,Pow(2,d) >;
template<class T,int d> using ArrayF3P=ArrayF<T,Pow(3,d) >;
////Other containers
template<class T> using List=std::list<T>;
template<class T,class CMP=std::less<T> > using Heap=std::priority_queue<T,std::vector<T>,CMP>;
template<class T1,class T2> using Pair=std::pair<T1,T2>;

////Eigen alias macros
#define Typedef_VectorD(d) \
using VectorD=Vector<real,d>
#define Typedef_VectorDi(d) \
using VectorDi=Vector<int,d>
#define Typedef_VectorDii(d) \
using VectorD=Vector<real,d>; \
using VectorDi=Vector<int,d>
#define Typedef_VectorEi(d) \
using VectorEi=Vector<int,d>;
#define Typedef_VectorDc(d) \
using VectorDc=Vector<C,d>
#define Typedef_MatrixD(d) \
using MatrixD=Matrix<real,d>
#define Typedef_VectorD_Alias(d,VectorD_Alias) \
using VectorD_Alias=Vector<real,d>
#define Typedef_VectorDi_Alias(d,VectorDi_Alias) \
using VectorDi_Alias=Vector<int,d>
#define Typedef_MatrixD_Alias(d,MatrixD_Alias) \
using MatrixD_Alias=Matrix<real,d>
#define Typedef_TransformD(d) \
using TransformD=Eigen::Transform<real,d,Eigen::Affine>; \
using TransformDf=Eigen::Transform<float,d,Eigen::Affine>;

////Eigen Vector1 compiler fix
Vector1 Vec1(const real s); 
Vector1i Vec1i(const int s);
Vector1d Vec1d(const double s);
Vector1f Vec1f(const float s);
Vector1c Vec1c(const C s);
////Eigen zero compiler fix
template<class T> T Zero(){return (T)0;}
#define Define_Zero(VecT) template<> inline VecT Zero<VecT>(){return VecT::Zero();}
Define_Zero(Matrix2);Define_Zero(Matrix3);Define_Zero(Matrix4);
Define_Zero(Vector1);Define_Zero(Vector2);Define_Zero(Vector3);
Define_Zero(Vector4);Define_Zero(Vector5);Define_Zero(Vector6);
Define_Zero(Vector7);Define_Zero(Vector8);Define_Zero(Vector9);

template <typename... Args>
void Assert(const bool flg, const char* fmt="", const Args &...args) {
    if (!flg) {
        std::cerr << "[Error]" << fmt::format(fmt, args...);
        exit(-1);
    }
}

template<typename ...Args>
void Error(const char* fmt_str, const Args&...args) {
    std::string msg = "#     " + fmt::format(fmt_str, args...) + "\n";
    fmt::print(fg(fmt::color::red), msg);
    throw msg;
}

template<typename ...Args>
void Info(const char* fmt, const Args&...args) {
    std::cout << "#     ";
    fmt::print(fmt, args...);
    std::cout << "\n";
}

void Info(const std::string& str);

//fmt adaptor for eigen vector
template <class T, int d> struct fmt::formatter<Vector<T, d> > {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        //https://fmt.dev/latest/api.html#udt
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && *it != '}') throw format_error("invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }

    // Formats the point p using the parsed format specification (presentation)
    // stored in this formatter.
    template <typename FormatContext>
    auto format(const Vector<T, d>& vec, FormatContext& ctx) -> decltype(ctx.out()) {
        std::stringstream ss;
        ss << vec.transpose();
        // ctx.out() is an output iterator to write to.
        return format_to(
            ctx.out(),
            "{}",
            ss.str());
    }
};

//fmt adaptor for eigen vector
template <class T> struct fmt::formatter<Eigen::Matrix<T, 3, 3> > {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        //https://fmt.dev/latest/api.html#udt
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && *it != '}') throw format_error("invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }

    // Formats the point p using the parsed format specification (presentation)
    // stored in this formatter.
    template <typename FormatContext>
    auto format(const Eigen::Matrix<T, 3, 3>& mat, FormatContext& ctx) -> decltype(ctx.out()) {
        std::stringstream ss;
        ss << mat;
        // ctx.out() is an output iterator to write to.
        return format_to(
            ctx.out(),
            "{}",
            ss.str());
    }
};

template<class T> struct fmt::formatter<Eigen::Quaternion<T>> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && *it != '}') throw format_error("invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }
    template <typename FormatContext>
    auto format(const Eigen::Quaternion<T>& quaternion, FormatContext& ctx) -> decltype(ctx.out()) {
        std::stringstream ss;
        ss << quaternion.coeffs().transpose();
        // ctx.out() is an output iterator to write to.
        return format_to(
            ctx.out(),
            "{}",
            ss.str());
    }
};

// Specialization for Eigen::Matrix<double, 3, 1>
// You can make this more generic for other Eigen types if needed.
template <>
struct fmt::formatter<Eigen::Matrix<double, 3, 1>> {
    // You can add support for parsing format specifiers if needed, e.g., precision.
    // For simplicity, this example uses a fixed format.
    constexpr auto parse(fmt::format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        // If you have custom presentation options, parse them here.
        // For now, we expect no custom options, just the closing '}'.
        if (it != end && *it != '}') {
            throw fmt::format_error("invalid format");
        }
        return it;
    }

    // This function defines how the Eigen::Matrix<double, 3, 1> is formatted.
    template <typename FormatContext>
    auto format(const Eigen::Matrix<double, 3, 1>& vec, FormatContext& ctx) const -> decltype(ctx.out()) {
        // Example format: "[x, y, z]"
        return fmt::format_to(ctx.out(), "[{}, {}, {}]", vec(0, 0), vec(1, 0), vec(2, 0));
        // For a column vector, vec(i) is often equivalent to vec(i,0).
        // Or, more compactly:
        // return fmt::format_to(ctx.out(), "[{}, {}, {}]", vec(0), vec(1), vec(2));
    }
};

template <>
struct fmt::formatter<Eigen::Matrix<double, 2, 1>> {
    // Parses the format string. Only accepts default {}
    constexpr auto parse(fmt::format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && *it != '}') {
            throw fmt::format_error("invalid format specifier for Eigen::Matrix<double, 2, 1> (only default {} is supported)");
        }
        return it;
    }

    // Formats the Eigen::Matrix<double, 2, 1> as "[x, y]"
    template <typename FormatContext>
    auto format(const Eigen::Matrix<double, 2, 1>& vec, FormatContext& ctx) const -> decltype(ctx.out()) {
        return fmt::format_to(ctx.out(), "[{}, {}]", vec(0), vec(1));
    }
};

template <>
struct fmt::formatter<std::filesystem::path> : fmt::formatter<std::string> {
    template <typename FormatContext>
    auto format(const std::filesystem::path& p, FormatContext& ctx) {
        return fmt::formatter<std::string>::format(p.string(), ctx);
    }
};

//if a is not nan or inf. Supports floating point types and Vector1~Vector3.
template<class T> bool Is_Valid_Number(const T& a);

void Print(const std::string& msg);

// CUDA programming
enum DataHolder { UNKNOWN=0, HOST, DEVICE };

// CPX programming
enum class SolverType { AUTO = 0, KRYLOV_CPU, MULTIGRID_AUTO, CPX_GPU };

#endif
