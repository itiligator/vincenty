#include "geopoint.h"

#include <iostream>
#include <iomanip>

int main(int, char **) {
    std::cout << std::setprecision(11);

    Geopoint SPB{-29.47124, 95.14681};
    Geopoint NY{-27.46601, -69.15955};
    std::cout << "Расстояние от Казанского собора до статуи Свободы: " << std::endl;
    std::cout << "В WSG84: " << std::endl;
    std::cout << " - по эллипсоиду: "  << NY.distanceTo(SPB) << " м" << std::endl;
    std::cout << " - по прямой: "  << NY.straightDistanceTo(SPB) << " м" << std::endl;

    Geopoint<PZ90> SPB_pz{-29.47124, 95.14681};
    Geopoint<PZ90> NY_pz{-27.46601, -69.15955};
    std::cout << "В ПЗ90: " << std::endl;
    std::cout << " - по эллипсоиду: "  << NY_pz.distanceTo(SPB_pz) << " м" << std::endl;
    std::cout << " - по прямой: "  << NY_pz.straightDistanceTo(SPB_pz) << " м" << std::endl;

    return 0;
}
