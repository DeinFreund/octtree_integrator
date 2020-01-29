#ifndef H_SOFTENING
#define H_SOFTENING

#include <bits/stdc++.h>

using namespace std;

class PlummerSoftening{
public:
  double getSoftening(const double dist, const double softening) const;
  double getForceSoftening(const double dist, const double softening) const;

};
class FasterSoftening{
public:
  double getSoftening(const double dist, const double softening) const;
  double getForceSoftening(const double dist, const double softening) const;

  
};

#endif
