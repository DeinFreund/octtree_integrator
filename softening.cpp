#include "softening.hpp"

double PlummerSoftening::getSoftening(const double dist, const double softening) const{

  return -1 * pow(dist * dist + softening * softening, -0.5);
}

double PlummerSoftening::getForceSoftening(const double dist, const double softening) const{

  return dist * pow(dist * dist + softening * softening, -1.5);
}
