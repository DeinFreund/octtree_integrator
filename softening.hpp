#ifndef H_SOFTENING
#define H_SOFTENING

#include <bits/stdc++.h>

using namespace std;

class PlummerSoftening{
public:
  float getSoftening(const float dist, const float softening) const;
  float getForceSoftening(const float dist, const float softening) const;

};
class FasterSoftening{
public:
  float getSoftening(const float dist, const float softening) const;
  float getForceSoftening(const float dist, const float softening) const;

  
};

#endif
