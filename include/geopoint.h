//
// Created by andrei on 28.03.23.
//

#ifndef TEMPLATE_GEOPOINT_H
#define TEMPLATE_GEOPOINT_H

#include <cmath>
#include <stdexcept>

#include "ellipsoid.h"


using std::atan, std::sin, std::cos, std::tan, std::abs, std::sqrt, std::pow;

/// \brief Реализация точки на эллипсоиде
/// \tparam T Эллипсоид
template<typename T = WSG84> requires IsEllipsoid<T>
class Geopoint {
public:
    using ellipsoid = T;

    /// \param lat широта в градусах. Допустимые значения от -90 до 90 градусов
    /// \param lon долгота в градусах. Допустимые значения от -360 до 360 градусов
    Geopoint(double lat, double lon) : m_lat(lat * M_PI / 180),
                                       m_lon(lon * M_PI / 180),
                                       m_lat_exact(lat * 1e9),
                                       m_lon_exact(lon * 1e9) {
        if (lat < -90 || lat > 90 || lon < -360 || lon > 360) {
            throw std::invalid_argument("Bad latitude or longitude provided");
        }
    };

    /// \brief Кратчайшая дистанция до точки. Реализовано по формулам Vincenty
    /// \param other Целевая точка
    /// \param N предельное количество итераций нахождения расстояния по долготе на вспомогательной сфере
    /// \param dL условие схождения итеративного поиска расстояния по долготе на вспомогательной сфере
    /// \return Дистанция в метрах
    double distanceTo(Geopoint<T> &other, uint N = 30, double dL = 1e-9) {
        if (operator==(other)) {
            return 0.0;
        } else {
            try {
                return vincentyInverseProblem(other, N, dL);
            }
            catch (const std::exception &) {
                return NAN;
            }
        }
    };

    /// \brief Дистанция до точки "по прямой", сквозь землю
    /// \param other Целевая точка
    /// \return Дистанция в метрах
    double straightDistanceTo(Geopoint<T> &other) {
        if (operator==(other)) {
            return 0.0;
        } else {
            double X1, Y1, Z1, X2, Y2, Z2;
            toECEF(X1, Y1, Z1);
            other.toECEF(X2, Y2, Z2);
            return sqrt(pow(X1 - X2, 2) + pow(Y1 - Y2, 2) + pow(Z1 - Z2, 2));
        }
    };

    bool operator==(const Geopoint<T> &other) {
        return (m_lat_exact == other.m_lat_exact && m_lon_exact == other.m_lon_exact);
    }


private:
    double m_lat;
    double m_lon;
    long long m_lat_exact; // (широта * 1e9), используется для различения точек
    long long m_lon_exact; // (долгота * 1e9) используется для различения точек


    // реализация метода Винценти для решения обратной проблемы
    double vincentyInverseProblem(Geopoint<T> &other, uint N = 30, double d_lambda = 1e-9) {
        // https://en.wikipedia.org/wiki/Vincenty%27s_formulae
        uint cnt{0}; // счетчик итераций

        // константы для расчета
        auto const U1 = atan((1 - ellipsoid::f) * tan(m_lat));
        auto const sin_U1 = sin(U1);
        auto const cos_U1 = cos(U1);
        auto const U2 = atan((1 - ellipsoid::f) * tan(other.m_lat));
        auto const sin_U2 = sin(U2);
        auto const cos_U2 = cos(U2);
        auto const L = other.m_lon - m_lon;

        // расстояние по широте между точками на вспомогательной сфере
        double lambda{L};

        // вспомогательная переменная для определения выполнения условия сходимости
        double lambda_prev{lambda + 2 * d_lambda};

        // вспомогательные и искомые переменные, уточняемые на каждом шаге итерации
        double sin_sigma, cos_sigma, sigma, sin_alpha, cos2_alpha, C, cos_2sigma_m;

        // итеративный поиск sigma (углового расстояния между точками)
        while ((cnt <= N || N == 0) && abs(lambda - lambda_prev) > abs(d_lambda)) {
            cnt++;
            sin_sigma = pow(cos_U2 * sin(lambda), 2) + pow(cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos(lambda), 2);
            sin_sigma = sqrt(sin_sigma);
            cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos(lambda);
            sigma = atan2(sin_sigma, cos_sigma);
            sin_alpha = cos_U1 * cos_U2 * sin(lambda) / sin_sigma;
            cos2_alpha = 1 - sin_alpha * sin_alpha;
            cos_2sigma_m = cos_sigma - 2 * sin_U1 * sin_U2 / cos2_alpha;
            C = ellipsoid::f * cos2_alpha * (4 + ellipsoid::f * (4 - 3 * cos2_alpha)) / 16;
            lambda_prev = lambda;
            lambda = L + (1 - C) *
                         ellipsoid::f *
                         sin_alpha *
                         (sigma +
                          C * sin_sigma * (cos_2sigma_m + C * cos_sigma * (2 * pow(cos_2sigma_m, 2) - 1)));

        }
        // если lambda не сошлась
        [[unlikely]] if (abs(lambda - lambda_prev) > abs(d_lambda))
            throw std::runtime_error("Vincenty algorithm fails to converge");

        // расчет дистанции по найденной sigma
        auto u2 = cos2_alpha * (pow(ellipsoid::a, 2) / pow(ellipsoid::b, 2) - 1);
        auto A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (360 - 175 * u2)));
        auto B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
        auto delta_sigma = B *
                           sin_sigma *
                           (cos_2sigma_m +
                            (B / 4) * (
                                    cos_sigma *
                                    (2 * pow(cos_2sigma_m, 2) - 1) -
                                    (B / 6) * cos_2sigma_m * (4 * pow(sin_sigma, 2) - 3) *
                                    (4 * pow(cos_2sigma_m, 2) - 3)
                            )
                           );
        return ellipsoid::b * A * (sigma - delta_sigma);
    }

    void toECEF(double &X, double &Y, double &Z) noexcept {
        // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
        double N = pow(ellipsoid::a, 2) / sqrt(
                pow(ellipsoid::a * cos(m_lat), 2) + pow(ellipsoid::b * sin(m_lat), 2)
        );
        X = N * cos(m_lat) * cos(m_lon);
        Y = N * cos(m_lat) * sin(m_lon);
        Z = pow((1 - ellipsoid::f), 2) * N * sin(m_lat);
    };
};


#endif //TEMPLATE_GEOPOINT_H
