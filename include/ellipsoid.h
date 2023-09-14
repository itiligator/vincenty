//
// Created by andrei on 12.09.23.
//

#ifndef VINCENTY_ELLIPSOID_H
#define VINCENTY_ELLIPSOID_H

#include <type_traits>
#include <cstdint>

/// \brief "Static" класс для хранения параметров эллипсоида
/// Помимо хранения параметров позволяет различать точки на различных эллипсоидах
/// \tparam _a Большая полуось
/// \tparam _f_inv величина, обратная сжатию _f_inv = 1e9 * 1/f
template<uint32_t _a, uint64_t _f_inv>
class Ellipsoid {
public:
    /// \brief большая полуось
    constexpr static double a{_a};
    /// \brief сжатие
    constexpr static double f{1e9 / _f_inv};
    /// \brief малая полуось
    constexpr static double b{(1 - f) * a};

    Ellipsoid() = delete;

};

/// \brief Реализует из WSG84
using WSG84 = Ellipsoid<6378137, 298257223563>;
/// \brief Реализует эллипсоид Красовского
using PZ90 = Ellipsoid<6378136, 298257840000>;

/// @private
template<class T>
struct is_ellipsoid_helper : std::false_type {
};
/// @private
template<uint32_t a, uint64_t f_inv>
struct is_ellipsoid_helper<Ellipsoid<a, f_inv>> : std::true_type {
};
/// \brief Концепт для проверки является ли переданный тип корректным эллипсоидом
/// \tparam T Проверяемый тип
template<class T>
concept IsEllipsoid = is_ellipsoid_helper<T>::value && (T::a > 6378132) && (T::a<
        6378142) && (T::f > 1e9 / 298259000000) && (T::f < 1e9 / 298256000000);


#endif //VINCENTY_ELLIPSOID_H