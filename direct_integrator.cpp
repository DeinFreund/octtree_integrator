#include "particle.hpp"
#include "softening.hpp"
#include "octtree.hpp"
#include <bits/stdc++.h>


using namespace std;



const float distBucketSize = 0.1;
const float softening = 1;
const float timestep_mult = 0.1;

template<typename T>
void integrateTimestep(const vector<Particle<T>*>& correct, const float timestep, const vector<Particle<T>*>& all){

  cerr << "Integrating " << correct.size() << " particles over " << timestep << "s" << endl;
  vector<Particle<T>*> defer;
  vector<Particle<T>*> keep;
  for (auto* p : correct){
    float optimalTimestep = timestep_mult * sqrt(0.1 / p->acc().norm());//softening
    if (optimalTimestep < timestep){
      defer.push_back(p);
    }else{
      keep.push_back(p);
    }
  }
  if (defer.empty()){
    for (auto* p : all){
      p->advanceTime(timestep);
    }
  }else{
    integrateTimestep(defer, timestep / 2, all);
    integrateTimestep(defer, timestep / 2, all);
  }
  Octtree<Particle<T>> tree(all, Eigen::Vector3f({-1000, -1000, -1000}), Eigen::Vector3f({1000, 1000, 1000}));
  float ekin = 0;
  float epot = 0;
//  float epot2 = 0;
  for (auto* p : all){
    ekin += p->kinetic_energy();
    epot += tree.getPotentials(p);
//    epot2 += p->getPotentials(all);
  }
  cerr << "Kinetic energy: " << ekin << endl;
  cerr << "Potential energy: " << epot << endl;
//  cerr << "Potential energy: " << epot2 << endl;
  cerr << "Total energy: " << ekin + epot << endl;
  cerr << all[1]->pos() << endl << all[2]->pos() << endl;
  bool anyExcessiveStep = 0;
  for (auto* p : keep){
    float acc = p->updateAcceleration(tree).norm();
    float optimalTimestep = timestep_mult * sqrt(0.1 / acc);//softening
    if (optimalTimestep < timestep) {
      anyExcessiveStep = true;
    }
  }/*
  if (anyExcessiveStep){

  }else{*/
    for (auto* p : keep){
      p->integrateStep();
    }
  
}

template<typename T>
vector<float> getAverageForce(const vector<vector<Particle<T>*>>& buckets, const vector<Particle<T>*>& all, const Octtree<Particle<T>>& octtree, bool plot = true){
  ofstream outfile;
  outfile.open("forceFunction.dat", ios::out | ios::trunc );
  Eigen::Vector3f origin({0,0,0});
  
  vector<float> avgForce(buckets.size());
  vector<float> analyticalForce(buckets.size());
  float totalMass = 0;
  float testMass = 1*all[3]->mass;
  cerr << setprecision(10);
  int nodesEvaluated = 0;
  for (size_t i = 0; i < buckets.size(); i++){
    float r = (i+0.5) * distBucketSize;
    //cerr << "mass is " << totalMass << " divider is " << r*r << endl;
    analyticalForce[i] = -totalMass * testMass / pow(r,2);
    for (Particle<T>* p : buckets[i]){
      //auto forces = p->getForces(all);
      //nodesEvaluated += all.size();
      //auto forces2 = octtree.getForces(p, false);
      auto forces = octtree.getForces(p, nodesEvaluated);
      avgForce[i] += forces.dot(p->pos().normalized());
      //cerr << *p << endl;
      //cerr << "force is " << forces << endl;
      //cerr << "force2 is " << forces2 << endl;
      //cerr << "force3 is " << forces3 << endl;
      //cerr << "directed force is " << forces.dot(p->pos().normalized()) << endl;
      totalMass += p->mass;
    }
    avgForce[i] /= buckets[i].size();
    cerr << i << ": " << analyticalForce[i] << " vs " << avgForce[i] << endl;
  }
  cerr << nodesEvaluated << " nodes evaluated" << endl;
  
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;
  
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << avgForce[i] << " ";
  }outfile << endl;
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;

  for (size_t i = 0; i < buckets.size(); i++){
    outfile << analyticalForce[i] << " ";
  }outfile << endl;
  outfile.close();

  if (plot){
    system("python3 plot.py < forceFunction.dat &");
  }
  return avgForce;

}

template<typename T>
vector<vector<Particle<T>*>> getDistanceFunction(const vector<Particle<T>*>& particles, bool plot = true){
  ofstream outfile;
  outfile.open("distanceFunction.dat", ios::out | ios::trunc );
  Eigen::Vector3f origin({0,0,0});
  
  vector<vector<Particle<T>*>> distBuckets(100);
  vector<float> massBuckets(distBuckets.size());
  vector<float> volume(distBuckets.size(), 0);
  float totalMass = 0;

  for (size_t i = 1; i < distBuckets.size(); i++){
    float r = i * distBucketSize;
    volume[i] = 4 * 0.33333 * M_PI * pow(r,3) - volume[i - 1];
  }
  for (Particle<T>* p : particles){
    totalMass += p->mass;
    int bucket = min((int)distBuckets.size()-1, int(p->distanceTo(origin) / distBucketSize));
    distBuckets[bucket].push_back(p);
    massBuckets[bucket] += p->mass;
  }
  vector<float> analytical;
  
  for (size_t i = 0; i < distBuckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;
  
  for (size_t i = 0; i < distBuckets.size(); i++){
    outfile << massBuckets[i] / volume[i] << " ";
  }outfile << endl;
  float a = 0.4;
  for (size_t i = 0; i < distBuckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;

  for (size_t i = 0; i < distBuckets.size(); i++){
    outfile << totalMass / 2 / 3.14159265 * a / (i+0.5)/distBucketSize / pow(distBucketSize * (i+0.5) + a,3) << " ";
  }outfile << endl;
  outfile.close();

  if (plot){
    system("python3 plot.py < distanceFunction.dat &");
  }
  return distBuckets;
}

template<typename T>
void readInput(vector<T*> & particles){
  size_t N;
  cin >> N >> N >> N;
  for (size_t i = 0; i < N; i++){
    particles.push_back(new T());
  }
  for (size_t i = 0; i < N; i++){
    cin >> particles[i]->mass;
  }
  for (size_t c = 0; c < 3; c++){
    for (size_t i = 0; i < N; i++){
      cin >> particles[i]->setPos()(c);
    }
  }
  for (size_t c = 0; c < 3; c++){
    for (size_t i = 0; i < N; i++){
      cin >> particles[i]->setVel()(c);
    }
  }
  for (size_t i = 0; i < N; i++){
    cin >> particles[i]->softening;
    particles[i]->softening = softening;
  }
  for (size_t i = 0; i < N; i++){
    cin >> particles[i]->potential;
  }
}

template<typename T>
void removeDuplicates(vector<T*> & particles){
  sort(particles.begin(), particles.end(), [](const T* a, const T* b){
      return a->pos()(0) > b->pos()(0)
	|| (a->pos()(0) == b->pos()(0) && a->pos()(1) > b->pos()(1))
	|| (a->pos()(0) == b->pos()(0) && a->pos()(1) == b->pos()(1) && a->pos()(2) > b->pos()(2));
    });
  for (int i = particles.size() - 2; i >= 0; i--){
    if (*particles[i] == *particles[i+1]) particles[i]->mass += particles[i+1]->mass;
  }
  auto last = particles.end();
  size_t j = 0;
  for (size_t i = 1; i < particles.size(); i++){
    if (*particles[i] == *particles[j]) {
      last--;
    }else{
      swap(particles[j+1], particles[i]);
      j++;
    }
  }
  cerr << "Removing " << particles.end() - last << " duplicate particles..\n";
  particles.erase(last, particles.end());
}

int main(){
  vector<Particle<PlummerSoftening>*> particles;
  readInput(particles);
  removeDuplicates(particles);
  Octtree<Particle<PlummerSoftening>> tree(particles, Eigen::Vector3f({-1000, -1000, -1000}), Eigen::Vector3f({1000, 1000, 1000}));
  //auto distBuckets = getDistanceFunction(particles, false);
  //auto avgForce = getAverageForce(distBuckets, particles, tree, false);
  cerr << "Built tree "<< endl;
  for (auto* p : particles){
    p->updateAcceleration(tree);
    p->integrateStep();
  }
  cerr << "Ready for integration "<< endl;
  integrateTimestep(particles, 10, particles);
}
