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
  float mass, potential, softening;
  S softener;

  float distanceTo(const Particle<S>* other) const{
    return (pos() - other->pos()).norm();
  }
  float distanceTo(const Eigen::Vector3f& other) const{
    return (pos() - other).norm();
  }

  bool inside(const Eigen::Vector3f& bottom,const Eigen::Vector3f& top) const{
    for (int i = 0; i < 3; i++) if (pos()(i) < bottom(i) || pos()(i) >= top(i)) return false;
    return true;
  }

  Eigen::Vector3f getForce(const Particle<S>* other) const{
    return (other->pos() - pos()).normalized() * mass * other->mass * softener.getForceSoftening((other->pos() - pos()).norm(),softening);
  }

  Eigen::Vector3f getForces(const vector<Particle<S>*> &particles) const{
    Eigen::Vector3f force({0,0,0});
    for (const auto* p : particles){
      if (p->pos() == pos()) continue;
      //cerr << "f:\n" << getForce(p) << endl;
      force += getForce(p);
    }
    return force;
  }
  float getPotential(const Particle<S>* other) const{
    return other->mass * softener.getSoftening((other->pos() - pos()).norm(),softening);
  }

  float getPotentials(const vector<Particle<S>*> &particles) const{
    float potential(0);
    for (const auto* p : particles){
      if (p->pos() == pos()) continue;
      potential += getPotential(p);
    }
    return potential;
  }

  void advanceTime(float dt){
    timePassed += dt;
  }

  Eigen::Vector3f updateAcceleration(const Octtree<Particle<S>>& tree){
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

  float kinetic_energy() const{
    return 0.5 * mass * vel().squaredNorm();
  }
  
  Eigen::Vector3f pos() const{
    return lastPos + timePassed * vel_halfstep();
  }

  Eigen::Vector3f vel() const{
    return lastVel + timePassed * lastAcc;
  }

  
  Eigen::Vector3f acc() const{
    return lastAcc;
  }

  
  Eigen::Vector3f vel_halfstep() const{
    return lastVel + timePassed / 2 * lastAcc;
  }

  Eigen::Vector3f& setPos(){
    return lastPos;
  }
  
  Eigen::Vector3f& setVel(){
    return lastVel;
  }

  float time() const{
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
  Eigen::Vector3f lastPos, lastVel, lastAcc = Eigen::Vector3f::Zero();
  Eigen::Vector3f newAcc = Eigen::Vector3f::Zero();
  float timePassed = 0;
  
};

#endif
