#ifndef H_OCTTREE
#define H_OCTTREE

#include <bits/stdc++.h>
#include <Eigen/Dense>

using namespace std;

template<typename T>
class Octtree{
public:
  Octtree(const vector<T*>& pars, const Eigen::Vector3f& bottom,const Eigen::Vector3f& top, const Octtree<T>* parent = NULL) : bottom(bottom), top(top), parent(parent)
  {
    if (!pars.empty()){
      for (T* p : pars){
	if (p->inside(bottom, top)){
	  particles.push_back(p);
	  constructionTime = p->time();
	}
      }
    }
    if (!particles.empty()){
      if (particles.size() > 1){
	Eigen::Vector3f subsize = 0.5 * (top - bottom);
	for (size_t dx = 0; dx < 2; dx++){
	  for (size_t dy = 0; dy < 2; dy++){
	    for (size_t dz = 0; dz < 2; dz++){
	      Eigen::Vector3f pos = bottom;
	      if (dx) pos(0) += subsize(0);
	      if (dy) pos(1) += subsize(1);
	      if (dz) pos(2) += subsize(2);
	      children[dx * 4 + dy * 2 + dz] = new Octtree(particles, pos, pos + subsize, this);
	    }
	  }
	}
      }
      //calculate center of mass and total mass
      center_of_mass << 0,0,0;
      velocity << 0,0,0;
      total_mass = 0;
      for (T* p : particles){
	center_of_mass += p->mass * p->cpos;
	//velocity += p->mass * p->vel();
	total_mass += p->mass;
      }
      center_of_mass *= 1.0/total_mass;
      velocity *= 1.0/total_mass;
      //calculate enclosing sphere diameter and relative positions
      diameter = 0;
      vector<Eigen::Vector3f> relative_pos(particles.size());
      for (size_t i = 0; i < particles.size(); i++){
	relative_pos[i] = (particles[i]->cpos - center_of_mass);
	diameter = max(diameter, 2 * relative_pos[i].norm());
      }
      
      //initialize quadrupole matrix
      quadrupole = Eigen::Matrix3f::Zero();
      for (size_t i = 0; i < 3; i++){
	float sum = 0;
	for (size_t k = 0; k < particles.size(); k++){
	  sum += particles[k]->mass * relative_pos[k].squaredNorm();
	}
	quadrupole(i,i) -= sum;
	for (size_t j = 0; j < 3; j++){
	  sum = 0;
	  for (size_t k = 0; k < particles.size(); k++){
	    sum += particles[k]->mass * relative_pos[k](i) * relative_pos[k](j);
	  }
	  quadrupole(i,j) += 3 * sum;
	}
      }
	
    }
    //cerr << "Initialized " << *this << endl;

  }

  ~Octtree() {
    if (particles.size() > 1){
      for (auto* o : children){
	delete o;
      }
    }
  }

  Octtree(const Octtree&) = delete;

  Eigen::Vector3f getForces(const T* particle, int& nodesEvaluated = 0, bool useQuadrupole = true) const{
    return getAccelerations(particle, nodesEvaluated, useQuadrupole) * particle->mass;
  }
  
  Eigen::Vector3f getAccelerations(const T* particle, int& nodesEvaluated = 0, bool useQuadrupole = true) const{

    if (particles.empty() || (particles.size() == 1 && particles[0]->cpos == particle->cpos)) return Eigen::Vector3f::Zero();
    const Eigen::Vector3f delta = particle->cpos - center_of_mass;
    if (particle->inside(bottom, top) || diameter / delta.norm() > ACCURACY_CRITERION){
      Eigen::Vector3f sum({0,0,0});
      for (const auto* o : children){
	sum += o->getAccelerations(particle, nodesEvaluated, useQuadrupole);
      }
      return sum;
    }

    Eigen::Vector3f curForce;
    if (diameter / delta.norm() < MONOPOLE_ACCURACY_CRITERION){
      curForce = getAcceleration(delta, particle->softening, false);
    }else{
      nodesEvaluated ++;
      curForce = getAcceleration(delta, particle->softening, useQuadrupole);
    }
    return curForce;
  }
  
  Eigen::Vector3f getAcceleration(const Eigen::Vector3f& delta, float softening, bool useQuadrupole = true) const
  {
    const float norm2 = delta.squaredNorm() + softening * softening;
    const float norm3 = pow(norm2, 3/2.);
    const Eigen::Vector3f monopoleGrad({
	-delta(0)/norm3,
	-delta(1)/norm3,
	-delta(2)/norm3,
	});
    if (useQuadrupole){
      
      const float norm5 = pow(norm2, 5/2.);
      const float norm7 = pow(norm2, 7/2.);
      const Eigen::Vector3f quadrupoleGrad({
	  (2*quadrupole(0,0)*delta(0)+quadrupole(0,1)*delta(1)+quadrupole(0,2)*delta(2)+quadrupole(1,0)*delta(1)+quadrupole(2,0)*delta(2))/norm5-(5*delta(0)*(delta(0)*(quadrupole(0,0)*delta(0)+quadrupole(1,0)*delta(1)+quadrupole(2,0)*delta(2))+delta(1)*(quadrupole(0,1)*delta(0)+quadrupole(1,1)*delta(1)+quadrupole(2,1)*delta(2))+delta(2)*(quadrupole(0,2)*delta(0)+quadrupole(1,2)*delta(1)+quadrupole(2,2)*delta(2))))/norm7,
	  (quadrupole(0,1)*delta(0)+quadrupole(1,0)*delta(0)+2*quadrupole(1,1)*delta(1)+quadrupole(1,2)*delta(2)+quadrupole(2,1)*delta(2))/norm5-(5*delta(1)*(delta(0)*(quadrupole(0,0)*delta(0)+quadrupole(1,0)*delta(1)+quadrupole(2,0)*delta(2))+delta(1)*(quadrupole(0,1)*delta(0)+quadrupole(1,1)*delta(1)+quadrupole(2,1)*delta(2))+delta(2)*(quadrupole(0,2)*delta(0)+quadrupole(1,2)*delta(1)+quadrupole(2,2)*delta(2))))/norm7,
	  (quadrupole(0,2)*delta(0)+quadrupole(1,2)*delta(1)+quadrupole(2,0)*delta(0)+quadrupole(2,1)*delta(1)+2*quadrupole(2,2)*delta(2))/norm5-(5*delta(2)*(delta(0)*(quadrupole(0,0)*delta(0)+quadrupole(1,0)*delta(1)+quadrupole(2,0)*delta(2))+delta(1)*(quadrupole(0,1)*delta(0)+quadrupole(1,1)*delta(1)+quadrupole(2,1)*delta(2))+delta(2)*(quadrupole(0,2)*delta(0)+quadrupole(1,2)*delta(1)+quadrupole(2,2)*delta(2))))/norm7
	  });

      return total_mass * monopoleGrad + 0.5 * quadrupoleGrad;
    }else{
      return total_mass * monopoleGrad;
    }
  }

  float getPotentials(const T* particle) const{

    if (particles.empty() || (particles.size() == 1 && *particles[0] == *particle)) return 0;
    const Eigen::Vector3f delta = particle->cpos - center_of_mass;
    if (particle->inside(bottom, top) || diameter / delta.norm() > ACCURACY_CRITERION){
      float sum = 0;
      for (const auto* o : children){
	sum += o->getPotentials(particle);
      }
      return sum;
    }

    return getPotential(delta, particle->softening);
  }
  
  float getPotential(const Eigen::Vector3f& delta, float softening) const
  {
    const float norm2 = delta.squaredNorm() + softening * softening;
    const float norm = pow(norm2, 1/2.);
    const float norm5 = pow(norm2, 5/2.);
    return - ( total_mass / norm + 0.5 * (delta.transpose() * quadrupole * delta)(0,0) / norm5);
  }

  
  
  Octtree<T>* getParticleNode(const T& p) const{
    if (!p.inside(bottom, top)) return NULL;
    Octtree<T>* ret = NULL;
    for (Octtree* o : children){
      ret = max(ret, o->getParticleNode(p));
    }
    return ret;
  }


  friend ostream& operator<< (ostream& stream, const Octtree<T>& octtree){
    stream << "Octtree of diameter " << octtree.diameter << " with " << octtree.particles.size() << " particles, constructed at " << octtree.constructionTime << " and center of mass: \n";
    stream << octtree.center_of_mass << "\nBottom: \n";
    stream << octtree.bottom << "\nTop: \n";
    stream << octtree.top << "\n";
    return stream;
  }

  
  const float ACCURACY_CRITERION = 0.1;//2.0;
  const float MONOPOLE_ACCURACY_CRITERION = 0.1 * ACCURACY_CRITERION;
  Eigen::Matrix3f quadrupole;
  Eigen::Vector3f center_of_mass;
  Eigen::Vector3f velocity;
  float total_mass, diameter;
  const Eigen::Vector3f bottom, top;
  vector<T*> particles;
private:
  const Octtree<T>* parent = NULL;
  array<Octtree<T>*, 8> children;
  float constructionTime = -1;
};

#endif
