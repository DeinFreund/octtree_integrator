#ifndef H_PARTICLE
#define H_PARTICLE

#include "octtree.hpp"
#include <bits/stdc++.h>
#include <Eigen/Dense>

using namespace std;

template<typename S>
class Particle
{
public:
  double mass, potential, softening;
  S softener;

  double distanceTo(const Particle<S>* other) const{
    return (pos() - other->pos()).norm();
  }
  double distanceTo(const Eigen::Vector3d& other) const{
    return (pos() - other).norm();
  }

  bool inside(const Eigen::Vector3d& bottom,const Eigen::Vector3d& top) const{
    for (int i = 0; i < 3; i++) if (pos()(i) < bottom(i) || pos()(i) >= top(i)) return false;
    return true;
  }

  Eigen::Vector3d getForce(const Particle<S>* other) const{
    return (other->pos() - pos()).normalized() * mass * other->mass * softener.getForceSoftening((other->pos() - pos()).norm(),softening);
  }

  Eigen::Vector3d getForces(const vector<Particle<S>*> &particles) const{
    Eigen::Vector3d force({0,0,0});
    for (const auto* p : particles){
      if (p->pos() == pos()) continue;
      //cerr << "f:\n" << getForce(p) << endl;
      force += getForce(p);
    }
    return force;
  }
  double getPotential(const Particle<S>* other) const{
    return other->mass * softener.getSoftening((other->pos() - pos()).norm(),softening);
  }

  double getPotentials(const vector<Particle<S>*> &particles) const{
    double potential(0);
    for (const auto* p : particles){
      if (p->pos() == pos()) continue;
      potential += getPotential(p);
    }
    return potential;
  }

  void advanceTime(double dt){
    timePassed += dt;
  }

  Eigen::Vector3d updateAcceleration(const Octtree<Particle<S>>& tree){
    int cost = 0;
    newAcc = tree.getAccelerations(this, cost);
    return newAcc;
  }
  
  void integrateStep(){
    
    lastPos = pos();
    lastVel = vel_halfstep() + newAcc * timePassed / 2;
    timePassed = 0;
    lastAcc = newAcc;
  }

  double kinetic_energy() const{
    return 0.5 * mass * vel().squaredNorm();
  }
  
  Eigen::Vector3d pos() const{
    return lastPos + timePassed * vel_halfstep();
  }

  Eigen::Vector3d vel() const{
    return lastVel + timePassed * lastAcc;
  }

  
  Eigen::Vector3d acc() const{
    return lastAcc;
  }

  
  Eigen::Vector3d vel_halfstep() const{
    return lastVel + timePassed / 2 * lastAcc;
  }

  Eigen::Vector3d& setPos(){
    return lastPos;
  }
  
  Eigen::Vector3d& setVel(){
    return lastVel;
  }

  double time() const{
    return timePassed;
  }
  
  friend bool operator== (const Particle<S>& l, const Particle<S>& r){
    return l.lastPos == r.lastPos;
  }
  
  friend ostream& operator<< (ostream& stream, const Particle<S>& particle) {
    stream << "Particle of mass " << particle.mass << " at position: \n";
    stream << particle.pos() << "\nVelocity: \n";
    stream << particle.vel() << "\nTime: " << particle.time() << "\n";
    return stream;
  }
  
private:
  Eigen::Vector3d lastPos, lastVel, lastAcc = Eigen::Vector3d::Zero();
  Eigen::Vector3d newAcc = Eigen::Vector3d::Zero();
  double timePassed = 0;
  
};

#endif
