/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // Create a random engine
  default_random_engine gen;

  // Set the number of particles
  num_particles = 100;

  // Creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {

    // Create a new particle
    Particle particle;

    // Sample random Gaussian noise from normal distrubtions using the random engine
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);

    // Initialize the weight with 1.0
    particle.weight = 1.0;

    // Add the new particle to the vector of particles
    particles.push_back(particle);
  }

  // Set the initialization status to true
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // Create a random engine
  default_random_engine gen;

  // Generate (Gaussian) noise for x, y and theta
  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);


  if (abs(yaw_rate) > 0.00001) {

    // Calculate theta_dot_dt as well as velocity divided by yaw_rate outside the loop to save calculation time
    double theta_dot_dt = yaw_rate * delta_t;
    double v_yr = velocity / yaw_rate;

    // Predict the new state for each particle and add the (Gaussian) noise
    for (Particle particle : particles) {
      particle.x += v_yr * (sin(particle.theta + theta_dot_dt) - sin(particle.theta)) + dist_x(gen); // Noise is added at the end
      particle.y += v_yr * (cos(particle.theta) - cos(particle.theta + theta_dot_dt)) + dist_y(gen); // Noise is added at the end
      particle.theta += theta_dot_dt; // Noise is added at the end
    }

  } else {

    // If the movement is nearly straight, avoid division by zero and use simpler formulas
    for (Particle particle : particles) {
      particle.x += velocity * cos(particle.theta) * delta_t + dist_x(gen); // Noise is added at the end
      particle.y += velocity * sin(particle.theta) * delta_t + dist_y(gen); // Noise is added at the end
      particle.theta += dist_theta(gen); // Noise is added at the end
    }
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


  double distance;
  double min_distance = numeric_limits<double>::max();
  LandmarkObs closest_observation;

  // For each predicted landmark
  for (LandmarkObs pred_landmark : predicted) {
    // Run through each observed landmark
    for (LandmarkObs ob_landmark : observations) {
      // Calculate the Euclidean distance
      distance = dist(pred_landmark.x, pred_landmark.y, ob_landmark.x, ob_landmark.y);
      // Memorize the closest observation
      if (distance < min_distance) {
        min_distance = distance;
        closest_observation = ob_landmark;
      }
    }
    // Assign the observed measurement to the predicted landmark
    pred_landmark.id = closest_observation.id;
    // Reset the minimum distance
    min_distance = numeric_limits<double>::max();
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Create a new vector for to store the survived particles
  vector<Particle> survived_particles;

  random_device rd;
  mt19937 gen(rd());
  discrete_distribution<> d(weights.begin(), weights.end());
  for(int i = 0; i < num_particles; i++) {
    // Add the survived particle to the vector
    survived_particles.push_back(particles[d(gen)]);
  }
  // Replace the particles with the resampled particles
  particles = survived_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
