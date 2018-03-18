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
#include <limits.h>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	if (is_initialized) {
    return;
  }

  // Initializing the number of particles
  num_particles = 100;

  // Extracting standard deviations
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Creating normal distributions, so that we can sample the particles
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  // Generate particles with normal distribution with mean on GPS values.
  // This is to get a good initial guess
  for (int i = 0; i < num_particles; i++) {

    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
	}

  // The filter is now initialized.
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// Extracting standard deviations

  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // Creating normal distributions
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);

  // Calculate new state.
  for (int i = 0; i < num_particles; i++) {

  	double theta = particles[i].theta;

    if ( fabs(yaw_rate) < 0.001 ) { // When yaw is not changing.
      particles[i].x += velocity * delta_t * cos( theta );
      particles[i].y += velocity * delta_t * sin( theta );
      // yaw continue to be the same.
    } else {
      particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
      particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate * delta_t;
    }

    // Adding noise.
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	// for(int i = 0 ; i < observations.size() ; i++){
	// 	double lowest_dist = INT_MAX;
	// 	int closest_landmark_id = -1;
	// 	double obs_x = observations[i].x;
	// 	double obs_y = observations[i].y;
	// 	for(int j = 0 ; j < predicted.size() ; j++){
	// 		double pred_x = predicted[j].x;
	// 		double pred_y = predicted[j].y;
	// 		int pred_id = predicted[j].id;
	// 		double current_dist = dist(obs_x,obs_y,pred_x,pred_y);
	// 		if(current_dist < lowest_dist){
	// 			lowest_dist = current_dist;
	// 			closest_landmark_id = pred_id;
	// 		}
	// 	}
	// 	observations[i].id = closest_landmark_id;
	// }
	int nObservations = observations.size();
	int nPredictions = predicted.size();
	for(int i = 0 ; i < nObservations ; i++){
		double minDistance = INT_MAX;
		int mid = -1;
		for(int j = 0 ; j < nPredictions ; j++){
			double xd = observations[i].x - predicted[j].x;
			double yd = observations[i].y - predicted[j].y;
			double distance = xd*xd + yd*yd;
			if(distance < minDistance){
				minDistance = distance;
				mid = predicted[j].id;
			}
		}
		observations[i].id = mid;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	// double weight_normalizer = 0.0;
	// for(int i = 0 ; i < num_particles ; i++){
	// 	double px = particles[i].x;
	// 	double py = particles[i].y;
	// 	double ptheta = particles[i].theta;

	// 	// transform the observations from the vehicle frame to the map frame
	// 	std::vector<LandmarkObs> transformed_observations;
	// 	for(int j = 0 ; j < observations.size() ; j++){
	// 		LandmarkObs transformed_obs;
	// 		transformed_obs.id = j;
	// 		transformed_obs.x = px + (cos(ptheta) * observations[j].x) - (sin(ptheta) * observations[j].y);
	// 		transformed_obs.y = py + (sin(ptheta) * observations[j].x) + (cos(ptheta) * observations[j].y);
	// 		transformed_observations.push_back(transformed_obs);
	// 	}

	// 	// filter to keep the data within the sensor range
	// 	vector<LandmarkObs> predicted_landmarks;
 //    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
 //      Map::single_landmark_s current_landmark = map_landmarks.landmark_list[j];
 //      if ((fabs((px - current_landmark.x_f)) <= sensor_range) && (fabs((py - current_landmark.y_f)) <= sensor_range)) {
 //        predicted_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
 //      }
 //    }

 //    // data association step
 //    dataAssociation(predicted_landmarks, transformed_observations);

 //    // calculate the weights of the particle
 //    weights[i] = 1.0;
 //    double sigma_x = std_landmark[0];
 //    double sigma_y = std_landmark[1];
 //    double sigma_x_2 = sigma_x * sigma_x;
 //    double sigma_y_2 = sigma_y * sigma_y;
 //    double normalizer = (1.0/(2.0 * M_PI * sigma_x * sigma_y));
 //    for(int a = 0 ; a < transformed_observations.size() ; a++){
 //    	double trans_obs_x = transformed_observations[a].x;
 //      double trans_obs_y = transformed_observations[a].y;
 //      double trans_obs_id = transformed_observations[a].id;
 //      double multi_prob = 1.0;
	//     for(int b = 0 ; b < predicted_landmarks.size() ; b++){
	//     	double pred_landmark_x = predicted_landmarks[b].x;
	//       double pred_landmark_y = predicted_landmarks[b].y;
	//       double pred_landmark_id = predicted_landmarks[b].id;

	//       if (trans_obs_id == pred_landmark_id) {
	//         multi_prob = normalizer * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
	//         particles[i].weight *= multi_prob;
	//       }
	//     }
 //    }
 //    weight_normalizer += particles[i].weight;
	// }
	// for (int i = 0; i < particles.size(); i++) {
 //    particles[i].weight /= weight_normalizer;
 //    weights[i] = particles[i].weight;
 //  }
	double stdLandmarkRange = std_landmark[0];
	double stdLandmarkBearing = std_landmark[1];
	for(int i = 0 ; i < num_particles ; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		// keeping only the landmarks that are in the sensor range
		double sensor_range_2 = sensor_range*sensor_range;
		vector<LandmarkObs> inRangeLandmarks;
		for(int j = 0 ; j < map_landmarks.landmark_list.size() ; j++){
			float landmarkX = map_landmarks.landmark_list[j].x_f;
			float landmarkY = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			double dX = x - landmarkX;
			double dY = y - landmarkY;
			if(dX*dX + dY*dY <= sensor_range_2){
				inRangeLandmarks.push_back(LandmarkObs{id,landmarkX,landmarkY});
			}
		}
		// transofrming to the map coordinates
		vector<LandmarkObs> mappedObservations;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mappedObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
    }
    // observation association to landmark
    dataAssociation(inRangeLandmarks, mappedObservations);
    // Reseting weight.
    particles[i].weight = 1.0;
    // Calculate weight
    for(int j = 0 ; j < mappedObservations.size() ; j++){
    	double observationX = mappedObservations[j].x;
      double observationY = mappedObservations[j].y;
      int landmarkId = mappedObservations[j].id;
      double landmarkX, landmarkY;
      unsigned int k = 0;
      unsigned int nLandmarks = inRangeLandmarks.size();
      bool found = false;
      while( !found && k < nLandmarks ) {
        if ( inRangeLandmarks[k].id == landmarkId) {
          found = true;
          landmarkX = inRangeLandmarks[k].x;
          landmarkY = inRangeLandmarks[k].y;
        }
        k++;
      }
      // Calculating weight.
      double dX = observationX - landmarkX;
      double dY = observationY - landmarkY;

      double weight = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp( -( dX*dX/(2*stdLandmarkRange*stdLandmarkRange) + (dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)) ) );
      if (weight == 0) {
        particles[i].weight *= 0.001;
      } 
      else {
        particles[i].weight *= weight;
      }
    }
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// vector<Particle> resampled_particles;
	// default_random_engine gen;
	// uniform_int_distribution<int> particle_index(0,num_particles-1);
	// int current_index = particle_index(gen);
	// double beta = 0.0;
	// double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
	// for(int i = 0 ; i < particles.size() ; i++){
	// 	uniform_int_distribution<double> random_weight(0.0,max_weight_2);
	// 	beta += random_weight(gen);
	// 	while(beta > weights[current_index]){
	// 		beta -= weights[current_index];
	// 		current_index = (current_index + 1) % num_particles;
	// 	}
	// 	resampled_particles.push_back(particles[current_index]);
	// }
	// particles = resampled_particles;
	// Get weights and max weight.
  vector<double> weights;
  double maxWeight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > maxWeight ) {
      maxWeight = particles[i].weight;
    }
  }

  // Creating distributions.
  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generating index.
  int index = distInt(gen);

  double beta = 0.0;

  // the wheel
  vector<Particle> resampledParticles;
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
