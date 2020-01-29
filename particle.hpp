#ifndef H_PARTICLE
#define H_PARTICLE

#include <Eigen/Dense>
#include <vector>

using namespace std;

template<typename S>
class Particle
{
public:
  Eigen::Vector3d pos, vel;
  double mass, potential, softening;
  S softener;

  template<typename T>
  double distanceTo(const Particle<T>& other) const{
    return (pos - other.pos).norm();
  }
  double distanceTo(const Eigen::Vector3d& other) const{
    return (pos - other).norm();
  }


  template<typename T>  
  Eigen::Vector3d getForce(const Particle<T> other) const{
    return (other.pos - pos).normalized() * mass * other.mass * softener.getForceSoftening((other.pos - pos).norm(),1);//softening);
  }

  template<typename T>
  Eigen::Vector3d getForces(const vector<Particle<T>> &particles) const{
    Eigen::Vector3d force;
    for (const auto& p : particles){
      if (p.pos == pos) continue;
      force += getForce(p);
    }
    return force;
  }

};

#endif
