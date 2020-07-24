/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using std::normal_distribution;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles    

  for (int i =0; i < num_particles; i++)
  {
    auto val = gaussians(x, y, theta, std);
    particles.emplace_back(Particle{i, val[0], val[1], val[2], 1});
    weights.push_back(1);
  }
  is_initialized = true;

}

std::vector<double> ParticleFilter::gaussians(double x, double y, double theta, double std[])
{
  normal_distribution<double> dist_x(x, std[0]);
  std::default_random_engine gen;
  // TODO: Create normal distributions for y and theta
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_t(theta, std[2]);
  float x_ = dist_x(gen);
  float y_ =  dist_y(gen);
  float theta_ = dist_t(gen);

  return vector<double>{x_, y_, theta_};

}
void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;
    /* Standard deviation values for x, y and theta*/
    double std_x{}, std_y{}, std_theta{};
    /* Get standard deviation values for x, y and theta */
    std_x = std_pos[0];
    std_y = std_pos[1];
    std_theta = std_pos[2];
    /* Create a noise (Gaussian noise) distribution for x, y and theta
    around mean 0 and standard deviation std_x, std_y and std_theta */
    std::normal_distribution<double> dist_x(0, std_x);
    std::normal_distribution<double> dist_y(0, std_y);
    std::normal_distribution<double> dist_theta(0, std_theta);
    /* Every particle is moved at certain distance at a certain heading after delta t */
  for (int i=0; i < num_particles; i++)
        {
            if (fabs(yaw_rate) < 0.001) {
                particles[i].x += velocity * delta_t * cos(particles[i].theta);
                particles[i].y += velocity * delta_t * sin(particles[i].theta);
            }
            else {
                particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
                particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
                particles[i].theta += yaw_rate * delta_t;
            }
            /* Add noise to every particle after upating it with motion */
            particles[i].x += dist_x(gen);
            particles[i].y += dist_y(gen);
            particles[i].theta += dist_theta(gen);
        };
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

    for (LandmarkObs& obs: observations)
    {
      double min = std::numeric_limits<double>::max();
      int id = -1;
  
      for (const LandmarkObs& pred: predicted)
      {
        double eucli_dist = dist(obs.x, obs.y, pred.x, pred.y);
        if (eucli_dist < min)
        {
          min = eucli_dist;
          id = pred.id;
        }

      }
      obs.id = id;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) 

{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (Particle& particle : particles)
  {
    double x = particle.x;
    double y = particle.y;
    double theta = particle.theta;
    //get predictions in sensor range 
    vector<LandmarkObs> predictions;

    // we care only about landmarks in a sensor range to particle
    for (auto& obs: map_landmarks.landmark_list)
    {
         if (fabs(obs.x_f - x) <= sensor_range && fabs(obs.y_f - y) <= sensor_range) 
          predictions.emplace_back(LandmarkObs{ obs.id_i, obs.x_f, obs.y_f });
    }

    //change to coordinate frame
    vector<LandmarkObs> os;
    for (auto& obs : observations) 
    {
      auto transformation = mapTransform(x, y, obs.x, obs.y, theta);
      os.push_back(LandmarkObs{obs.id, transformation[0], transformation[1]});
    }

    dataAssociation(predictions, os);
    
    particle.weight = 1.0;

    for(LandmarkObs& obs : os)
    {
         LandmarkObs temp;
         for(LandmarkObs& pred : predictions)
         {
           if(obs.id == pred.id)
               temp=pred;
         }
         double sig_x = std_landmark[0];
         double sig_y = std_landmark[1];
         particle.weight *= multiv_prob(sig_x, sig_y, temp.x, temp.y, obs.x, obs.y);
    }
  }
}  

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   std::vector<Particle> particles_2;
   std::vector<double> weight;
   std::default_random_engine gen;

   for (int i = 0; i < num_particles; i++) {
    weight.push_back(particles[i].weight);
   }

   std::discrete_distribution<> dist(weight.begin(), weight.end());

   for(int i = 0; i < num_particles; i++) {
        particles_2.push_back(particles[dist(gen)]);
   }

   particles=particles_2;
}
double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

std::vector<double> ParticleFilter::mapTransform(double xp, double yp, double xc, double yc, double theta)
{
  double xm = xp + (cos(theta) * xc) - (sin(theta)*yc);
  double ym = yp + (sin(theta) * xc) + (cos(theta)*yc);

  return vector<double>{xm, ym};
}
void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}