#include "particle.hpp"
#include "softening.hpp"
#include "octtree.hpp"
#include <bits/stdc++.h>


using namespace std;



const float distBucketSize = 0.1;
const float softening = 0;
const float timestep_mult = 0.1;//e-12;


template<typename T>
void exportPositions(const vector<Particle<T>*>& all, int id){
  ofstream outfile;
  outfile.open(string("frame_") + to_string(id) + string(".dat"), ios::out | ios::trunc );
  Eigen::Vector3f bottom({-0.5f, -0.5f, -0.5f});
  Eigen::Vector3f top({0.5f, 0.5f, 0.5f});
  top*= 0.1;
  bottom*= 0.1;
  float scale = 100;
  for (const auto* p : all){
    if (p->inside(bottom, top)){
      outfile << scale * p->pos()(0) << "," << scale * p->pos()(1) << "," << scale * p->pos()(2) << "\n" ;
    }
  }
  outfile.close();

}



template<typename T>
void integrateTimestep(const vector<Particle<T>*>& correct, const float timestep, const vector<Particle<T>*>& all){

  cerr << "Integrating " << correct.size() << " particles over " << timestep / 6.7034898045672115724987120078450556863038856284648168e-8 << " years" << endl;
  vector<Particle<T>*> defer;
  vector<Particle<T>*> keep;
  for (auto* p : correct){
    float optimalTimestep = timestep_mult * sqrt(0.1 / p->acc().norm());//softening
    //cerr << optimalTimestep << endl;
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
  cerr << all[0]->pos() << endl;
  cerr << all[1]->pos() << endl;

    for (auto* p : keep){
      p->updateAcceleration(tree).norm();
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
  float avgError =0;
  for (size_t i = 0; i < buckets.size(); i++){
    float r = (i+0.5) * distBucketSize;
    double mass = totalMass;
    //cerr << "mass is " << totalMass << " divider is " << r*r << endl;
    for (Particle<T>* p : buckets[i]){
      //auto forcesReal = p->getForces(all);
      //nodesEvaluated += all.size();
      //auto forces2 = octtree.getForces(p, false);
      auto forces = octtree.getForces(p, nodesEvaluated);
      avgForce[i] += forces.dot(p->pos().normalized());
      //avgError += ((forces - forcesReal)/forcesReal.norm()).squaredNorm();
      //cerr << *p << endl;
      //cerr << "force is " << forces << endl;
      //cerr << "force2 is " << forcesReal << endl;
      //cerr << "force3 is " << forces3 << endl;
      //cerr << "directed force is " << forces.dot(p->pos().normalized()) << endl;
      totalMass += p->mass;
      if (p->pos().norm() < r) mass += p->mass;
    }
    analyticalForce[i] = -mass * testMass / pow(r,2);
    avgForce[i] /= buckets[i].size();
    cerr << i << ": " << analyticalForce[i] << " vs " << avgForce[i] << endl;
  }
  cerr << nodesEvaluated << " nodes evaluated" << endl;
  avgError = sqrt(avgError / all.size());
  cerr << "Mean square error: " << avgError << endl;
  /*
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize * 100 << " ";
  }outfile << endl;

  float forceMult = 2.95002701690675868785668327030469360673831619299e14;
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << avgForce[i] * forceMult  << " ";
  }outfile << endl;
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize * 100 << " ";
  }outfile << endl;

  for (size_t i = 0; i < buckets.size(); i++){
    outfile << analyticalForce[i]  * forceMult<< " ";
  }outfile << endl;
  outfile.close();
  */
  if (plot){
    system("python3 plot.py < forceFunction.dat &");
  }
  return avgForce;

}


template<typename T>
vector<float> getAverageVelocity(const vector<vector<Particle<T>*>>& buckets, const vector<Particle<T>*>& all, bool plot = true){
  ofstream outfile;
  outfile.open("velocityFunction.dat", ios::out | ios::trunc );
  Eigen::Vector3f origin({0,0,0});
  
  vector<float> avgVelocity(buckets.size());
  vector<float> analyticalVelocity(buckets.size());
  float totalMass = 0;
  cerr << setprecision(10);
  int nodesEvaluated = 0;
  for (size_t i = 0; i < buckets.size(); i++){
    float r = (i+0.5) * distBucketSize;
    double mass = totalMass;
    //cerr << "mass is " << totalMass << " divider is " << r*r << endl;
    for (Particle<T>* p : buckets[i]){
      //auto velocitys = p->getVelocitys(all);
      //nodesEvaluated += all.size();
      //auto velocitys2 = octtree.getVelocitys(p, false);
      auto velocitys = p->vel();
      //avgVelocity[i] += velocitys.cross(p->pos().normalized()).norm();
      avgVelocity[i] += velocitys.norm();
      //cerr << *p << endl;
      //cerr << "velocity is " << velocitys << endl;
      //cerr << "velocity2 is " << velocitys2 << endl;
      //cerr << "velocity3 is " << velocitys3 << endl;
      //cerr << "directed velocity is " << velocitys.dot(p->pos().normalized()) << endl;
      totalMass += p->mass;
      if (p->pos().norm() < r) mass += p->mass;
    }
    avgVelocity[i] /= buckets[i].size();
    analyticalVelocity[i] = sqrt(mass / (r));
    //cerr << buckets[i].size() << endl;
    //cerr << i << ": " << analyticalVelocity[i] << " vs " << avgVelocity[i] << endl;
  }
  cerr << nodesEvaluated << " nodes evaluated" << endl;
  
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize * 100 << " ";
  }outfile << endl;
  float velMult = 11.84562361769985255465093916276684402846336303609205718315;
  
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << avgVelocity[i] * velMult << " ";
  }outfile << endl;
  for (size_t i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize  * 100<< " ";
  }outfile << endl;

  for (size_t i = 0; i < buckets.size(); i++){
    outfile << analyticalVelocity[i] *velMult<< " ";
  }outfile << endl;
  outfile.close();

  if (plot){
    system("python3 plot.py < velocityFunction.dat &");
  }
  return avgVelocity;

}


template<typename T>
float getHalfMassRadius(const vector<Particle<T>*>& particles){
  float totalMass = 0;
  for (const auto* p : particles){
    totalMass += p->mass;
  }
  float s = 0;
  float e = 1000;
  while (e-s > 1e-5){
    float m = (e+s)/2;
    float m2 = m*m;
    float mass = 0;
    for (const auto* p : particles){
      if (p->pos().squaredNorm() < m2)
      mass += p->mass;
    }
    if (mass > 0.5 * totalMass){
      e = m;
    }else{
      s = m;
    }
  }
  return (e+s)/2;
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

  for (size_t i = 0; i < distBuckets.size(); i++){
    float r = (i+1) * distBucketSize;
    volume[i] = 4 * 0.33333 * M_PI * pow(r,3);
  }
  for (int i = distBuckets.size(); i > 0; i--){
    float r = (i+1) * distBucketSize;
    volume[i] = 4 * 0.33333 * M_PI * pow(r,3) - volume[i - 1];
  }
  for (Particle<T>* p : particles){
    totalMass += p->mass;
    int bucket = min((int)distBuckets.size()-1, int(p->distanceTo(origin) / distBucketSize));
    distBuckets[bucket].push_back(p);
    massBuckets[bucket] += p->mass;
  }
  vector<float> analytical;
  
  for (size_t i = 0; i < distBuckets.size()-1; i++){
    outfile << i * distBucketSize * 100 << " ";
  }outfile << endl;
  
  for (size_t i = 0; i < distBuckets.size()-1; i++){
    outfile << (massBuckets[i] / volume[i]) / 1e6<< " ";
  }outfile << endl;
  //float a = 0.1222;
  float a = 0.09;
  //cin >> a;
  for (size_t i = 0; i < distBuckets.size()-1; i++){
    outfile << i * distBucketSize *100<< " ";
  }outfile << endl;

  for (size_t i = 0; i < distBuckets.size()-1; i++){
    float r = (i+0.5)*distBucketSize;
    outfile << (totalMass / 2 / 3.14159265 * a / r / pow(r + a,3)) / 1e6 << " ";
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
vector<T*> reduce(const vector<T*> &particles, float fraction){
  vector<int> indexes(particles.size());
  iota(indexes.begin(), indexes.end(), 0);
  random_shuffle(indexes.begin(), indexes.end());
  int N = lround(particles.size() * fraction);
  vector<T*> ret;
  ret.reserve(N);
  for (int i = 0; i < N; i++){
    ret.push_back(particles[indexes[i]]);
  }
  for (size_t i = N; i < particles.size(); i++){
    delete particles[indexes[i]];
  }
  return ret;
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
  for (auto* p : particles){
    p->advanceTime(0);
  }
//particles = reduce(particles, 0.1);
  //exportPositions(particles, 0);
  //cerr << "Half mass radius is " << getHalfMassRadius(particles) << endl;
  Octtree<Particle<PlummerSoftening>> tree(particles, Eigen::Vector3f({-1000, -1000, -1000}), Eigen::Vector3f({1000, 1000, 1000}));
  cerr << "Built tree "<< endl;
  auto distBuckets = getDistanceFunction(particles, false);
  auto avgForce = getAverageForce(distBuckets, particles, tree, false);
  /*
  auto avgVel = getAverageVelocity(distBuckets, particles, true);
  for (auto* p : particles){
    p->updateAcceleration(tree);
    p->integrateStep();
  }
  cerr << "Ready for integration "<< endl;
  integrateTimestep(particles, 6.7034898045672115724987120078450556863038856284648168e-8, particles);
//*/
}
