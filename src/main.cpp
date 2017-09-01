#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <utility>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define POINTS_TO_SEND 50 // Numer of points to send to the simulator
#define DT 0.02 // 20 ms in s
#define THRES 7.5 // Max acc and max jerk
#define POINTS_TO_KEEP 50
#define MAX_VEL 22.15 // 50 mph is 22.352 mps
#define MAX_S 6945.554 // Length of a lap
#define STEP_SIZE 45.0 // Step for path planning in m
#define SAFETY_ZONE STEP_SIZE/2 // m between cars which we consider as safe for lane change

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
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
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

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
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


// Power function for int. https://stackoverflow.com/questions/29787310/does-pow-work-for-int-data-type-in-c
int int_pow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
           result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

class Node {
	public:
		int lane; // Frequency
		double speed; // current speed
		double c_speed; // cummulative speed
		Node *p; // Perent of the node
	Node(){
	}
	Node(int ln, double sp, double score, Node *par){ // Constructor
		lane =  ln;
	    speed = sp;
	    c_speed = speed*score + par->c_speed; // Add score to sum speed
	    p = par;
	}
	Node(int ln, double sp){ // Constructor
		lane =  ln;
	    speed = sp;
	    c_speed = speed;
	    p = NULL;
	}
};
// Compare for Node class, it needs for sorting
struct MyCompare {
    bool operator()(const Node* l, const Node* r) const {
		return l->c_speed > r->c_speed;
	}
};
// Check if it is safe to change lines 
bool unsafe(double car_s, double car_speed, double s_front, double speed_front, double s_back, double speed_back){
	double t = STEP_SIZE*2 / car_speed;
	double car_s_shifted = car_s + t * car_speed;
	if (fabs(car_s_shifted - (s_front + t * speed_front + car_s)) < SAFETY_ZONE || fabs(car_s_shifted - (car_s - s_front + t * speed_front)) < SAFETY_ZONE){
		return true;
	}
	return false;
}
// Calculate speed of a car, if we prefer to use the lane 'des_lane'
double speed_of_lane(vector<vector<double>> sens, double car_s, double dt, double cur_speed, int cur_lane, int des_lane, double tl){
	if (abs(cur_lane - des_lane) > 1){  // Prohibited maneuver, only one lane shift is allowed
		return 0;
	}
	double speed_of_lane = MAX_VEL;
	double speed_front = MAX_VEL;
	double speed_back = 0.0;
	double s_front = STEP_SIZE; 
	double s_back = STEP_SIZE;
	double car_s_shifted = car_s + (dt - tl) * cur_speed;
	for (int i =0; i < sens.size(); i++){
		double d = sens[i][6];
		if (d < 4*(des_lane+1) && d > 4*des_lane){
			double vx = sens[i][3];
			double vy = sens[i][4];
			double speed = sqrt(vx*vx+vy*vy);
			double s = sens[i][5];
			s += dt * speed;
			if (s > MAX_S){ // We have a round track!
				s -= MAX_S;
			}
			if (s > car_s_shifted && (s - car_s_shifted) < s_front){
				s_front = s - car_s_shifted; // Calculate distance for the car in front
				if (speed < speed_front){
					speed_front = speed;
				}
			}
			if (s < car_s_shifted && (car_s_shifted - s) < s_back){
				s_back = car_s_shifted - s;
				if (speed > speed_back){
					speed_back = speed;
				}
			}	
		}
	}
	if (des_lane != cur_lane && unsafe(car_s_shifted, cur_speed, s_front, speed_front, s_back, speed_back)){ // Unsafe to change lane
		speed_of_lane = 0;
	}
	else{
		speed_of_lane = max(speed_front, speed_back);
	}
	return speed_of_lane;
}

// Return the next lane and desired speed for it.
pair<int, double> path_plan(int l, double car_s, double dt, double cur_speed, vector<vector<double>> sens, int depth){
	int next_lane = l;
	double tl = dt;
	double next_speed = cur_speed;
	int tot = 0; // Total numbers of Nodes
	for (int i = 0; i < depth; i++){
		tot += int_pow(3, i);
	}
	const int elem = tot;
	vector<Node*> path_map(elem); // 3d tree
	Node *first_node = new Node(l, cur_speed);  // The top of the tree
	path_map[0] = first_node;
	int p = 1; // counter
	double est_speed = cur_speed;
	for (int i = 1; i < depth; i++){ // Build the tree
		dt += STEP_SIZE / est_speed;
		int layer_size = int_pow(3, i);
		est_speed = 0;
		for (int j = 0; j < layer_size; j++){
			Node* parent = path_map[(p-1)/3];
			int des_lane = j % 3;
			double speed = speed_of_lane(sens, car_s, dt, parent->speed, parent->lane, des_lane, tl);
			if (speed > est_speed){
				est_speed = speed;
			}
			double score = 0.9; // Punish lane changes
			if (des_lane == parent->lane){
				score = 1.0;
			}
			score *= 1.0 - 0.05*des_lane; // Left line is better because it is closer to the center line of the track, hence, we can complete a lap faster
			Node *new_node = new Node(des_lane, speed, score, parent);
			path_map[p] = new_node;
			p++;
		}
	}
	sort(path_map.begin(), path_map.end(), MyCompare()); // Sort elements
	for (int i = 0; i < elem; i++){ // Find the best path
		Node* node = path_map[i];
		if (node->speed <= 0.001){ // if the node is unreachable
			continue;
		}
		bool path_found = true;
		while ((node->p)->p != NULL){
			node = node->p;
			if (node->speed <= 0.001){ // if the node is unreachable
				path_found = false;
				break; // Go to the next option
			}
		} 
		if (path_found){
			next_lane = node->lane;
			next_speed = node->speed;
			break;
		}
	}
	cout << next_lane << " " << path_map[0]->c_speed << " " << path_map[1]->c_speed << endl;
	pair<int, double> res = {next_lane, next_speed};
	return res;
}

// Check if planned maneuver if safe and return speed to keep
double check(double car_s, int lane, int next_lane, double cur_speed, vector<vector<double>> sens, double dt){
	double next_speed = MAX_VEL;
	int minl = min(lane, next_lane);
	int maxl = max(lane, next_lane);
	double speed_front = MAX_VEL;
	double speed_back = 0.0;
	double s_front = STEP_SIZE;
	double s_back = STEP_SIZE;
	for (int i =0; i < sens.size(); i++){
		double d = sens[i][6];
		if (d < 4*(maxl+1) && d > 4*minl){
			double vx = sens[i][3];
			double vy = sens[i][4];
			double speed = sqrt(vx*vx+vy*vy);
			double s = sens[i][5];
			s += dt * speed;
			if (s > MAX_S){ // We have a round track!
				s -= MAX_S;
			}
			if (s > car_s && (s - car_s) < s_front){
				s_front = s - car_s;
				if (speed < speed_front){
					speed_front = speed;
				}
			}
			if (s < car_s && (car_s - s) < s_back){
				s_back = car_s-s;
				if (speed > speed_back){
					speed_back = speed;
				}
			}	
		}
	}
	if (unsafe(car_s, cur_speed, s_front, speed_front, s_back, speed_back)){
		next_speed = 0;
	}
	else{
		next_speed = max(speed_back, speed_front);
	}
	return next_speed;
}

int main() {
  uWS::Hub h;
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;
  double car_speed_p = 0; // car acceleration on the prev step
  double car_acc_p = 0; // car acceleration on the prev step
  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = MAX_S;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

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
  
  int lane = 1;
  double max_vel = MAX_VEL; // mps
  double cur_speed = 0; // current speed
  double cur_acc = 0;
  double c_time = 0;
  bool change = false;
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
  &car_speed_p, &car_acc_p, &lane, &max_vel, &cur_speed, &max_s, &cur_acc, &c_time, &change](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
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
          	
          	int prev_size = previous_path_x.size();
          	
          	double end_s;
          	if (prev_size>0){
				end_s = end_path_s;
			}
			const double t_step = (POINTS_TO_SEND-prev_size) * DT; 
          	double max_acc_pos = THRES+cur_acc;
          	if (max_acc_pos > THRES){
				max_acc_pos = THRES;
			}
          	double max_acc_neg = (cur_acc-THRES);
          	if (max_acc_neg < -THRES){
				max_acc_neg = -THRES;
			}
			if (c_time < 2.5 ){ // Deal with the start lag
				max_acc_pos /= 15 * (2.5-c_time);
				max_acc_neg /= 15 * (2.5-c_time);
			}
			// Choose the next lane
			double lane_speed = max_vel;
			int next_lane = lane;
			if (!change){ // If we are in a lane (not changing), we try to fing a new better one
				pair<int, double> res = path_plan(lane, car_s, (double)prev_size * DT, cur_speed, sensor_fusion, 4);
				next_lane = res.first;
				if (next_lane != lane){
					// Check for safety and set new lane speed
					double pol_speed = check(car_s, lane, next_lane, cur_speed, sensor_fusion, (double)prev_size * DT);
					if (pol_speed > 1){
						lane_speed = pol_speed;
						cout << next_lane << endl;
						change = true;
					}
					else{
						next_lane = lane;
						lane_speed = cur_speed;
						change = false;
					}
				}
			}
			double speed_back = 0.0;
			// Find the best speed
			if (!change){
				for (int i = 0; i < sensor_fusion.size(); i++){
					double d = sensor_fusion[i][6];
					if (d < 4*(lane+1) && d > 4*lane){ // A car to consider
						double vx = sensor_fusion[i][3];
						double vy = sensor_fusion[i][4];
						double speed = sqrt(vx*vx+vy*vy);
						double s = sensor_fusion[i][5];
						s += (double)prev_size * DT * speed;
						if (s > max_s){ // We have a round track!
							s -= max_s;
						}
						if (s > car_s && (s - car_s) < STEP_SIZE){
							lane_speed = speed;
							if ((s - car_s) < SAFETY_ZONE){
								lane_speed = speed*0.95; // Break if we are too close!
								cout << "Break!" << endl;
							}
							break;
						}
					}
				}
			}
			else{
				int minl = min(lane, next_lane);
				int maxl = max(lane, next_lane);
				double speed_front = MAX_VEL;
				double speed_back = 0.0;
				double s_front = STEP_SIZE;
				double s_back = STEP_SIZE;
				for (int i =0; i < sensor_fusion.size(); i++){
					double d = sensor_fusion[i][6];
					if (d < 4*(maxl+1) && d > 4*minl){
						double vx = sensor_fusion[i][3];
						double vy = sensor_fusion[i][4];
						double speed = sqrt(vx*vx+vy*vy);
						double s = sensor_fusion[i][5];
						s += (double)prev_size * DT * speed;
						if (s > MAX_S){ // We have a round track!
							s -= MAX_S;
						}
						if (s > car_s && (s - car_s) < s_front){
							s_front = s - car_s;
							if (speed < speed_front){
								speed_front = speed;
							}
						}
						if (s < car_s && (car_s - s) < s_back){
							s_back = car_s-s;
							if (speed > speed_back){
								speed_back = speed;
							}
						}	
					}
				}
				lane_speed = max(speed_front, speed_back);
			}
          	double prop_speed = cur_speed;
          	double p_acc = cur_acc;
			if (cur_speed < lane_speed){
				prop_speed += max_acc_pos*t_step;
				cur_acc = max_acc_pos;
				if (prop_speed > lane_speed){ // If it is too much
					prop_speed = lane_speed;
					cur_acc = (lane_speed-cur_speed)/t_step;
				}
				if (prop_speed > max_vel){
					prop_speed = max_vel;
					cur_acc = (max_vel-cur_speed)/t_step;
				}
			}
			else{
				prop_speed += max_acc_neg*t_step;
				cur_acc = max_acc_neg;
			}
          	cur_speed = prop_speed;

          	// Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            vector<double> ptsx;
            vector<double> ptsy;
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);        
            if (prev_size < 2){
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);
				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);
				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			}
			else {
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];
				double prev_ref_x = previous_path_x[prev_size-2];
				double prev_ref_y = previous_path_y[prev_size-2];
				ref_yaw = atan2(ref_y-prev_ref_y,ref_x-prev_ref_x);
				ptsx.push_back(prev_ref_x);
				ptsx.push_back(ref_x);
				ptsy.push_back(prev_ref_y);
				ptsy.push_back(ref_y);
			}
			
			const int s_step = 3;
			const double s_dist = STEP_SIZE;
			for (int i=0; i < s_step; i++){
				vector<double> wp = getXY(car_s+s_dist*(i+1), (2 + 2 * (next_lane + lane)), map_waypoints_s, map_waypoints_x, map_waypoints_y);
				ptsx.push_back(wp[0]);
				ptsy.push_back(wp[1]);
			}
			
			double sin_yaw = sin(0 - ref_yaw);
			double cos_yaw = cos(0 - ref_yaw);
			for (int i = 0; i < ptsx.size(); i++){
				double shift_x = ptsx[i]-ref_x;
				double shift_y = ptsy[i]-ref_y;
				ptsx[i]=shift_x*cos_yaw-shift_y*sin_yaw;
				ptsy[i]=shift_y*cos_yaw+shift_x*sin_yaw;
			}
			tk::spline s;  // Create spline
			s.set_points(ptsx, ptsy); 
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
          	int point_transfer = min(prev_size, POINTS_TO_KEEP);
			for (int i = 0; i < point_transfer; i++){
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}
			if (fabs(4.0 * (double)next_lane - car_d + 2) < 0.25){ // We are ready for a new lane change
				change = false;
			}
			
			lane = next_lane;
			double target_y = s(s_dist);
			double dist = sqrt(target_y*target_y+s_dist*s_dist);
			double x_add_on = 0;
			double N = dist / (DT * cur_speed);
			double cos_y = cos(ref_yaw);
			double sin_y = sin(ref_yaw);
			for (int i = 0; i < POINTS_TO_SEND-point_transfer; i++){
				double x_point = x_add_on+s_dist/N;
				double y_point = s(x_point);
				x_add_on = x_point;
				double x_ref = x_point;
				double y_ref = y_point;
				x_point = x_ref * cos_y - y_ref * sin_y + ref_x;
				y_point = x_ref * sin_y + y_ref * cos_y + ref_y;
				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}
			json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            c_time += t_step;
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
















































































