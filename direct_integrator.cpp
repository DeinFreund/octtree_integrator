#include "particle.hpp"
#include "softening.hpp"
#include <bits/stdc++.h>


using namespace std;



double distBucketSize = 0.1;

template<typename T>
vector<double> getAverageForce(const vector<vector<Particle<T>>>& buckets, const vector<Particle<T>>& all, bool plot = true){
  ofstream outfile;
  outfile.open("forceFunction.dat", ios::out | ios::trunc );
  Eigen::Vector3d origin({0,0,0});
  
  vector<double> avgForce(buckets.size());
  vector<double> analyticalForce(buckets.size());
  double totalMass = 0;
  double testMass = 1*all[0].mass;
  for (int i = 0; i < buckets.size(); i++){
    double r = (i+0.5) * distBucketSize;
    cerr << "mass is " << totalMass << " divider is " << r*r << endl;
    analyticalForce[i] = -totalMass * testMass / pow(r,2);
    for (Particle<T> p : buckets[i]){
      avgForce[i] += p.getForces(all).dot(p.pos.normalized());
      totalMass += p.mass;
    }
    avgForce[i] /= buckets[i].size();
    cerr << i << ": " << analyticalForce[i] << " vs " << avgForce[i] << endl;
  }
  
  for (int i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;
  
  for (int i = 0; i < buckets.size(); i++){
    outfile << avgForce[i] << " ";
  }outfile << endl;
  for (int i = 0; i < buckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;

  for (int i = 0; i < buckets.size(); i++){
    outfile << analyticalForce[i] << " ";
  }outfile << endl;
  outfile.close();

  if (plot){
    system("python3 plot.py < forceFunction.dat &");
  }
  return avgForce;

}

template<typename T>
vector<vector<Particle<T>>> getDistanceFunction(const vector<Particle<T>>& particles, bool plot = true){
  ofstream outfile;
  outfile.open("distanceFunction.dat", ios::out | ios::trunc );
  Eigen::Vector3d origin({0,0,0});
  
  vector<vector<Particle<T>>> distBuckets(100);
  vector<double> massBuckets(distBuckets.size());
  vector<double> volume(distBuckets.size(), 0);
  double totalMass = 0;

  for (int i = 1; i < distBuckets.size(); i++){
    double r = i * distBucketSize;
    volume[i] = 4 * 0.33333 * M_PI * pow(r,3) - volume[i - 1];
  }
  for (const auto& p : particles){
    totalMass += p.mass;
    int bucket = min((int)distBuckets.size()-1, int(p.distanceTo(origin) / distBucketSize));
    distBuckets[bucket].push_back(p);
    massBuckets[bucket] += p.mass;
  }
  vector<double> analytical;
  
  for (int i = 0; i < distBuckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;
  
  for (int i = 0; i < distBuckets.size(); i++){
    outfile << massBuckets[i] / volume[i] << " ";
  }outfile << endl;
  double a = 0.4;
  for (int i = 0; i < distBuckets.size(); i++){
    outfile << i * distBucketSize << " ";
  }outfile << endl;

  for (int i = 0; i < distBuckets.size(); i++){
    outfile << totalMass / 2 / 3.14159265 * a / (i+0.5)/distBucketSize / pow(distBucketSize * (i+0.5) + a,3) << " ";
  }outfile << endl;
  outfile.close();

  if (plot){
    system("python3 plot.py < distanceFunction.dat &");
  }
  return distBuckets;
}

template<typename T>
void readInput(vector<Particle<T>> & particles){
  int N = particles.size();
  for (int i = 0; i < N; i++){
    cin >> particles[i].mass;
  }
  for (int c = 0; c < 2; c++)
    for (int i = 0; i < N; i++){
      cin >> particles[i].pos(c);
    }
  for (int c = 0; c < 2; c++)
    for (int i = 0; i < N; i++){
      cin >> particles[i].vel(c);
    }
  for (int i = 0; i < N; i++){
    cin >> particles[i].softening;
  }
  for (int i = 0; i < N; i++){
    cin >> particles[i].potential;
  }
}

int main(){
  int N;
  cin >> N >> N >> N;
  vector<Particle<PlummerSoftening>> particles(N);
  readInput(particles);
  auto distBuckets = getDistanceFunction(particles);
  auto avgForce = getAverageForce(distBuckets, particles);
}
