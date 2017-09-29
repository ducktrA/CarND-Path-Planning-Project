#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define LANE(l) (2. + l * 4.)
#define PREF_SIZE 20
#define MAX_LANES 3.0
#define TARGET_VELOCITY 22.5

#define PENALTY_NOT_ON_RIGHT_MOST_LANE 1e4
#define PENALTY_ON_DISTANCE 1e4
#define PENALTY_ON_VELOCITY 1e3

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

double notKeepingRightPenalty(double lane, vector<double> pack, vector<vector<double>> sensor_fusion) {
	// 	if lane right to my one is empty return penalty. else not
	//vector<double> pack = {lane, car_s, car_d};
	return (2 - lane) * 200;
}

double lowDistancePenalty(double lane, vector<double> pack, vector<vector<double>> sensor_fusion) {
	
	const double CLEAR_ROAD = 150;
	double closestCarAhead = CLEAR_ROAD;
	double car_s = pack[1];

	for(int i = 0; i < sensor_fusion.size(); i++) {
  		/*  ["sensor_fusion"] A 2d vector of cars and then that car's 
  			[0][car's unique ID, 
  			[1]car's x position in map coordinates, 
  			[2]car's y position in map coordinates, 
  			[3]car's x velocity in m/s, 
  			[4]car's y velocity in m/s, 
  			[5]car's s position in frenet coordinates, 
  			[6]car's d position in frenet coordinates. 

  			vx = sf[i][3]
  			vy = sf[i][4]
  			s = sf[i][5]
			d = sf[i][6]
  		*/
  		double d = sensor_fusion[i][6];

  		// is a vehicle in my lane? I occupy the range [(car_d) - 2, (car_d) + 2]
  		if(d > (LANE(lane)) - 2 && d < (LANE(lane)+2)){
  			//YES, how far away is it?
  			double vx = sensor_fusion[i][3], vy = sensor_fusion[i][4];
  			double other_veh_speed = sqrt(pow(vx,2)+pow(vy,2));
  			double other_veh_s = sensor_fusion[i][5];

  			if(other_veh_s > car_s){
  				double dist = other_veh_s - car_s;
  				closestCarAhead = min(closestCarAhead, dist);
  			}
  		}
  	}

	return PENALTY_ON_DISTANCE * ((CLEAR_ROAD-closestCarAhead)/CLEAR_ROAD);
} 

double adjustSpeed(double lane, vector<double> pack, vector<vector<double>> sensor_fusion) {
	const double CLEAR_ROAD = 300;
	double closestCarAhead = CLEAR_ROAD;
	double speedCarAhead = TARGET_VELOCITY;
	double ref_target = speedCarAhead;
	double car_s = pack[1];

	for(int i = 0; i < sensor_fusion.size(); i++) {
  		double d = sensor_fusion[i][6];

  		// is a vehicle in my lane? I occupy the range [(car_d) - 2, (car_d) + 2]
  		if(d > (LANE(lane)) - 2 && d < (LANE(lane)+2)){
  			//YES, how far away is it?
  			double vx = sensor_fusion[i][3], vy = sensor_fusion[i][4];
  			double other_veh_speed = sqrt(pow(vx,2)+pow(vy,2));
  			double other_veh_s = sensor_fusion[i][5];

  			if(other_veh_s > car_s){
  				double dist = other_veh_s - car_s;
  				closestCarAhead = min(closestCarAhead, dist);
  				speedCarAhead = min(speedCarAhead, other_veh_speed);
  			}
  		}
  	}

	if(closestCarAhead > 50.)
		ref_target = TARGET_VELOCITY;
	if(closestCarAhead < 30.)
		ref_target = speedCarAhead;
	if(closestCarAhead < 20.)
		ref_target = speedCarAhead - 3;

  	return ref_target;
}

double slowLanePenalty(double lane, vector<double> pack, vector<vector<double>> sensor_fusion) {
	double speedCarAhead = TARGET_VELOCITY;
	double ref_target = speedCarAhead;
	double car_s = pack[1];

	for(int i = 0; i < sensor_fusion.size(); i++) {
  		double d = sensor_fusion[i][6];

  		// is a vehicle in my lane? I occupy the range [(car_d) - 2, (car_d) + 2]
  		if(d > (LANE(lane)) - 2 && d < (LANE(lane)+2)){
  			//YES, how far away is it?
  			double vx = sensor_fusion[i][3], vy = sensor_fusion[i][4];
  			double other_veh_speed = sqrt(pow(vx,2)+pow(vy,2));
  			double other_veh_s = sensor_fusion[i][5];

  			if(other_veh_s > car_s){
  				speedCarAhead = min(speedCarAhead, other_veh_speed);
  			}
  		}
  	}
  	return (TARGET_VELOCITY - speedCarAhead) * PENALTY_ON_VELOCITY;
}

// this function estimates how much other traffic participants have when the ego performs an action. How much time is left assuming the follower vehicle does an emergency break down to the velocity of the ege vehicle minus an optimistic 0.5 second reaction time
// the higher the values the better
double reactionTime(double lane, vector<double> pack, vector<vector<double>> sensor_fusion) {

	double car_s = pack[1];
	double car_speed = pack[3];
	double ret = 5.0;

	for(int i = 0; i < sensor_fusion.size(); i++) {
  		
  		double d = sensor_fusion[i][6];

  		// is a vehicle in my lane? I occupy the range [(car_d) - 2, (car_d) + 2]
  		if(d > (LANE(lane) -2) && d < (LANE(lane) +2)){
  			//YES, how far away is it?
  			double vx = sensor_fusion[i][3], vy = sensor_fusion[i][4];
  			double other_veh_speed = sqrt(pow(vx,2)+pow(vy,2));
  			double other_veh_s = sensor_fusion[i][5];
  			
  			if((other_veh_s > car_s) && (other_veh_s - car_s) < 25.) {
  				ret = 0.5;
  			}

  			if((other_veh_s < car_s)) {
  				double distanceBackwards = car_s - other_veh_s;
  				//printf("Distance Backwards %f\n", distanceBackwards);

  				const double a = 7.0;
  				
  				double bdist = (pow(other_veh_speed,2) - pow(car_speed,2)) / (2*a);
  				double ttcoll = ((distanceBackwards - bdist) / car_speed) - 0.5; // 0.5 seconds to react

  				if(ttcoll < 0)
  					ttcoll = 0;

  				ret = min(ret, ttcoll);
			}
  		}
  	}
	return ret;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  double lane = 1.0;
  double ref_v = 0;
  double ref_target = TARGET_VELOCITY;

  vector<double> costs = {0.0, 0.0, 0.0};
  vector<int> supporters = {0,0,0};
  vector<double> ttcollision = {0.0, 0.0, 0.0};


  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  printf("map.waypoints_x.size(): %ld\n", map_waypoints_x.size());

  h.onMessage([&lane, &costs, &supporters, &ref_v, &ref_target, &ttcollision, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	vector<double> pack = {lane, car_s, car_d, car_speed};

          	json msgJson;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	
          	int prev_size = previous_path_x.size();

          	if(prev_size > 0){
          		car_s = end_path_s;
          	}
          	
          	ref_target = adjustSpeed(lane, pack, sensor_fusion);
          	
          	if(ref_v < ref_target - 0.5){
          		ref_v += .2;
          	}
          	if(ref_v > ref_target){
          		ref_v -= 0.5;
          	}

          	costs = {0.0, 0.0, 0.0};

          	//determine costs
          	for(int i = 0; i < MAX_LANES; i++) {
          		double toLane = double(i);
          		// a penalty on not keeping right when possible leads to local minima the car can only get out with a double lane change.
          		// this requires different checks if both lane is clear
          		//costs[i] += notKeepingRightPenalty(toLane, pack, sensor_fusion);
				costs[i] += lowDistancePenalty(toLane, pack, sensor_fusion);
				costs[i] += slowLanePenalty(toLane, pack, sensor_fusion);

				ttcollision[i] = reactionTime(toLane, pack, sensor_fusion);
          	}

          	printf("Lane Penalty: %f %f %f\n", costs[0], costs[1], costs[2]);
          	printf("Lane Clear: %f %f %f\n", ttcollision[0], ttcollision[1], ttcollision[2]);

          	std::vector<double>::iterator costs_iter = std::min_element(std::begin(costs), std::end(costs));
          	int min_cost_element = std::distance(std::begin(costs), costs_iter);
          	
          	//printf("min element %d\n", min_cost_element);
          	//collect supporters
          	// max and min iterators return first element. 
          	// In case the lane the ego is in has the same costs, that will cause a unnecessary lane change

          	if(costs[lane] == costs[min_cost_element]){
          		min_cost_element = lane;
          	}

         	supporters[min_cost_element]++;

          	std::vector<int>::iterator supporters_iter = std::max_element(std::begin(supporters), std::end(supporters));
          	int max_supporters_element = std::distance(std::begin(supporters), supporters_iter);

          	if(supporters[max_supporters_element] > 5) {

          		// determine intended direction of movement
          		// in this case the vehicle is right to the new target lane, check if ttcollision[lane-1] is large enough
          		if(abs(lane - max_supporters_element) > 1) {
          			printf("double lane change - requires different safety metric");
          			// here we actuall wait until a single lane change oppertunity comes up. might result getting stuck in local minima
          		}
          		else{

          		if(lane > max_supporters_element) {
          			// a lane change is only allowed if the following traffic has enough time to react on the new situation
          			if(ttcollision[lane-1] > 2.) {
          				printf("Lane Penalty: %f %f %f\n", costs[0], costs[1], costs[2]);
          				printf("Lane Clear: %f %f %f\n", ttcollision[0], ttcollision[1], ttcollision[2]);

          				printf("change lane to the left\n");
          				lane = lane - 1;
          			}
          		}

          		if(lane < max_supporters_element) {
          			if(ttcollision[lane+1] > 2.) {
          				printf("Lane Penalty: %f %f %f\n", costs[0], costs[1], costs[2]);
          				printf("Lane Clear: %f %f %f\n", ttcollision[0], ttcollision[1], ttcollision[2]);

          				printf("change lane to the right\n");
          				lane = lane + 1;
          			}
          		}
          	}
				//lane = max_supporters_element;
				supporters = {0,0,0};
          	}

          	// evaluating the sensor_fusion end

          	vector<double> splinepoints_x, splinepoints_y;
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);

          	if(prev_size < 2){
          		// get two points so the spline attaches to the current heading of the vehicle
          		splinepoints_x.push_back(car_x - cos(car_yaw));
          		splinepoints_y.push_back(car_y - sin(car_yaw));

          		splinepoints_x.push_back(car_x);
          		splinepoints_y.push_back(car_y);
          	}

          	// instead of attaching to the end of previous_path, lets cut it after a set of points to be more responsive
          	else{
          		// same here, but we can access the previous_path vector
          		ref_x = previous_path_x[prev_size-1];
          		ref_y = previous_path_y[prev_size-1];

          		const double ref_x_prev = previous_path_x[prev_size-2];
          		const double ref_y_prev = previous_path_y[prev_size-2];

          		ref_yaw = atan2(ref_y-ref_y_prev,ref_x - ref_x_prev);

          		// the "older" first, the newer second
          		splinepoints_x.push_back(ref_x_prev);
          		splinepoints_y.push_back(ref_y_prev);

          		// 
          		splinepoints_x.push_back(ref_x);
          		splinepoints_y.push_back(ref_y);
          	}

			// spline generation depends on reference velocity. If intended to got fast, stretch the anchor waypoints
          	for(double s_step = 30.; s_step < 120.0; s_step += 10 + 1.2 * ref_target){
          		vector<double> xy = getXY(car_s + s_step, LANE(lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          		splinepoints_x.push_back(xy[0]);
          		splinepoints_y.push_back(xy[1]);
          	}
          	
          	// now bring points in line with car heading. Necessary to create uniform stepping

          	//shift & rotate
			for(int i = 0; i < splinepoints_x.size(); i++) {
				double shift_x = splinepoints_x[i] - ref_x;
				double shift_y = splinepoints_y[i] - ref_y;

				splinepoints_x[i] = shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw);
				splinepoints_y[i] = shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw);
			}

			tk::spline s;
/*
			printf("splinepoints:\n");
			for(int i = 0; i < splinepoints_x.size(); i++) {
				printf("x: %f y: %f\n", splinepoints_x[i], splinepoints_y[i]);
			}
*/
			s.set_points(splinepoints_x, splinepoints_y);

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// linear approximation of point distance
          	double target_x = 30.0;
          	double target_y = s(target_x);
          	double target_dist = sqrt(pow(target_x,2) + pow(target_y,2));

          	double x_add_on = 0.0;

          	for(int i = 0; i < prev_size; i++) {
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}
    
          	for(int i = 1; i <= 75 - prev_size; i++) {
          		double N = target_dist/(0.02 * ref_v);

          		double x_point = x_add_on + (target_x / N);
          		double y_point = s(x_point);

          		x_add_on = x_point;

          		// reverse rotate & shift
          		double x_p = (x_point * cos(ref_yaw) - y_point * sin(ref_yaw)) + ref_x;
          		double y_p = (x_point * sin(ref_yaw) + y_point * cos(ref_yaw)) + ref_y;

          		next_x_vals.push_back(x_p);
          		next_y_vals.push_back(y_p);
          	}

          	// check the values and their distance
          	//double diff_avg = 0.0, diff_max = 0;
          	//int index;

            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
