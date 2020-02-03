#include "softening.hpp"

float PlummerSoftening::getSoftening(const float dist, const float softening) const{

  return -1 * pow(dist * dist + softening * softening, -0.5);
}

float PlummerSoftening::getForceSoftening(const float dist, const float softening) const{

  return dist * pow(dist * dist + softening * softening, -1.5);
}
