/*
 * Circle.cpp
 * RVO2 Library
 *
 * Copyright 2008 University of North Carolina at Chapel Hill
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Please send all bug reports to <geom@cs.unc.edu>.
 *
 * The authors may be contacted via:
 *
 * Jur van den Berg, Stephen J. Guy, Jamie Snape, Ming C. Lin, Dinesh Manocha
 * Dept. of Computer Science
 * 201 S. Columbia St.
 * Frederick P. Brooks, Jr. Computer Science Bldg.
 * Chapel Hill, N.C. 27599-3175
 * United States of America
 *
 * <http://gamma.cs.unc.edu/RVO2/>
 */

/*
 * Example file showing a demo with 250 agents initially positioned evenly
 * distributed on a circle attempting to move to the antipodal position on the
 * circle.
 */

#ifndef RVO_OUTPUT_TIME_AND_POSITIONS
#define RVO_OUTPUT_TIME_AND_POSITIONS 1
#endif

#include <cmath>
#include <cstddef>
#include <vector>

#if RVO_OUTPUT_TIME_AND_POSITIONS
#include <iostream>
#endif

#if _OPENMP
#include <omp.h>
#endif

#include <RVO.h>

#ifndef M_PI
const float M_PI = 3.14159265358979323846f;
#endif

#define nLeft 50
#define nRight 100
int posRange = 30;
int posMin = 0;
int getThreLeft(int t, int hurry);
int getThreRight(int t, int hurry);

/* Store illegal crossing parameter. */
int illegalR = 0; // right road segment, including both directions
int illegalL = 0; // left road segment, including both directions
// seqR, seqL, nLeft, and nRight are the only stuff need adjusting when flow rate changes
// seqR, seqL format: {hurryLevel, timeAppear}
int seqR[nRight][2] = {{8, 0}, {7, 30}, {6, 60}, {1, 90}, {8, 120}, {3, 150}, {2, 180}, {8, 210}, {4, 240}, {4, 270},
					{9, 300}, {0, 330}, {5, 360}, {5, 390}, {6, 420}, {8, 450}, {7, 480}, {5, 510}, {0, 540}, {4, 570},
					{2, 600}, {5, 630}, {1, 660}, {5, 690}, {1, 720}, {5, 750}, {4, 780}, {4, 810}, {8, 840}, {10, 870},
					{4, 900}, {7, 907}, {8, 915}, {9, 922}, {9, 930}, {9, 937}, {5, 945}, {10, 952}, {8, 960}, {1, 967},
					{6, 975}, {4, 982}, {0, 990}, {1, 997}, {3, 1005}, {6, 1012}, {2, 1020}, {4, 1027}, {6, 1035}, {7, 1042},
					{1, 1050}, {3, 1080}, {0, 1110}, {8, 1140}, {2, 1170}, {10, 1200}, {4, 1230}, {6, 1260}, {2, 1290}, {3, 1320},
					{1, 1350}, {4, 1380}, {9, 1410}, {5, 1440}, {8, 1470}, {3, 1500}, {1, 1530}, {10, 1560}, {3, 1590}, {4, 1620},
					{7, 1650}, {0, 1680}, {9, 1710}, {8, 1740}, {9, 1770}, {9, 1800}, {8, 1830}, {3, 1860}, {8, 1890}, {4, 1920},
					{4, 1950}, {7, 1980}, {3, 2010}, {3, 2040}, {7, 2070}, {10, 2100}, {10, 2130}, {3, 2160}, {9, 2190}, {0, 2220},
					{0, 2250}, {6, 2280}, {3, 2310}, {3, 2340}, {1, 2370}, {7, 2400}, {10, 2430}, {4, 2460}, {4, 2490}, {1, 2520}};
// int seqL[nLeft][2] = {{8, 0}, {7, 30}, {6, 60}, {1, 90}, {8, 120}, {3, 150}, {2, 180}, {8, 210}, {4, 240}, {4, 270},
// 					{9, 300}, {0, 330}, {5, 360}, {5, 390}, {6, 420}, {8, 450}, {7, 480}, {5, 510}, {0, 540}, {4, 570},
// 					{2, 600}, {5, 630}, {1, 660}, {5, 690}, {1, 720}, {5, 750}, {4, 780}, {4, 810}, {8, 840}, {10, 870},
// 					{4, 900}, {7, 907}, {8, 915}, {9, 922}, {9, 930}, {9, 937}, {5, 945}, {10, 952}, {8, 960}, {1, 967},
// 					{6, 975}, {4, 982}, {0, 990}, {1, 997}, {3, 1005}, {6, 1012}, {2, 1020}, {4, 1027}, {6, 1035}, {7, 1042}};
int seqL[nLeft][2] = {{8, 45}, {4, 105}, {1, 165}, {10, 225}, {8, 285}, {1, 345}, {9, 405}, {4, 465}, {5, 525}, {2, 585},
					{2, 645}, {1, 705}, {4, 765}, {1, 825}, {7, 885}, {8, 915}, {4, 930}, {2, 945}, {9, 960}, {5, 975},
					{0, 990}, {7, 1005}, {9, 1020}, {3, 1035}, {6, 1050}, {3, 1095}, {2, 1155}, {2, 1215}, {7, 1275}, {1, 1335},
					{0, 1395}, {1, 1455}, {10, 1515}, {3, 1575}, {6, 1635}, {4, 1695}, {9, 1755}, {6, 1815}, {6, 1875}, {3, 1935},
					{5, 1995}, {10, 2055}, {5, 2115}, {2, 2175}, {7, 2235}, {7, 2295}, {4, 2355}, {7, 2415}, {7, 2475}, {8, 2535}};

// store the priority of initial waiting positions
float posL[][2] = {{-14.6, 3.1}, {-15.3, -4.0}, {-15.1, 4.6}, {-15.1, 1.9}, {-15.0, -2.0}, {-15.5, -0.7},
					{-14.4, 0.6}, {-14.3, -3.2}, {-15.3, 3.0}, {-14.3, 2.1}, {-14.8, 4.0}, {-14.6, -1.0},
					{-15.4, 0.4}, {-14.2, -0.1}, {-16.0, 1.1}, {-16.3, -2.5}, {-15.6, -3.1}, {-16.5, -3.1},
					{-14.7, 5.5}, {-16.1, -3.7}, {-15.8, 2.3}, {-15.9, 4.0}, {-16.2, -0.1}, {-14.9, -5.5},
					{-14.7, -4.5}, {-14.9, -7.1}, {-15.1, 7.2}, {-14.2, 5.9}, {-14.4, -6.7}, {-15.1, 6.1},
					{-16.6, 0.4}, {-16.8, 4.3}, {-15.1, -8.0}, {-16.0, -6.2}, {-15.6, -2.4}, {-15.5, -5.4},
					{-16.7, 6.0}, {-15.6, 7.9}, {-16.7, -6.4}, {-14.4, -8.4}, {-16.8, -4.5}, {-16.4, -5.3},
					{-16.3, 7.1}, {-17.0, 1.7}, {-14.7, 8.0}, {-16.2, -8.1}, {-17.4, 5.7}, {-17.1, -5.2},
					{-15.9, -7.1}, {-15.8, 5.6}, {-17.7, 2.7}, {-17.3, -2.9}, {-17.1, -7.1}, {-16.8, -1.7}, 
					{-17.1, -7.9}, {-16.5, 8.1}, {-16.7, 3.5}, {-17.0, -0.4}, {-17.1, 7.0}, {-17.4, 1.0}, 
					{-17.4, -3.7}, {-17.4, -6.0}, {-15.7, -1.6}, {-16.9, 2.5}, {-16.3, -1.0}, {-16.7, 5.1}};
float posR[][2] = {{10.4, -3.1}, {11.1, 4.0}, {10.9, -4.6}, {10.9, -1.9}, {10.8, 2.0}, {11.3, 0.7},
					{10.2, -0.6}, {10.1, 3.2}, {11.1, -3.0}, {10.1, -2.1}, {10.6, -4.0}, {10.4, 1.0},
					{11.2, -0.4}, {10.0, 0.1}, {11.8, -1.1}, {12.1, 2.5}, {11.4, 3.1}, {12.3, 3.1},
					{10.5, -5.5}, {11.9, 3.7}, {11.6, -2.3}, {11.7, -4.0}, {12.0, 0.1}, {10.7, 5.5},
					{10.5, 4.5}, {10.7, 7.1}, {10.9, -7.2}, {10.0, -5.9}, {10.2, 6.7}, {10.9, -6.1},
					{12.4, -0.4}, {12.6, -4.3}, {10.9, 8.0}, {11.8, 6.2}, {11.4, 2.4}, {11.3, 5.4},
					{12.5, -6.0}, {11.4, -7.9}, {12.5, 6.4}, {10.2, 8.4}, {12.6, 4.5}, {12.2, 5.3},
					{12.1, -7.1}, {12.8, -1.7}, {10.5, -8.0}, {12.0, 8.1}, {13.2, -5.7}, {12.9, 5.2},
					{11.7, 7.1}, {11.6, -5.6}, {13.5, -2.7}, {13.1, 2.9}, {12.9, 7.1}, {12.6, 1.7},
					{12.9, 7.9}, {12.3, -8.1}, {12.5, -3.5}, {12.8, 0.4}, {12.9, -7.0}, {13.2, -1.0},
					{13.2, 3.7}, {13.2, 6.0}, {11.5, 1.6}, {12.7, -2.5}, {12.1, 1.0}, {12.5, -5.1},
					{13.5, 1.2}, {13.7, 4.8}, {13.7, 0.1}, {13.3, 2.1}, {13.7, -7.3}, {13.4, -3.7},
					{13.5, -4.7}, {13.6, -6.5}, {13.7, 7.8}, {13.7, 6.7}, {11.8, 4.6}, {13.7, -1.8}};
// storing the goals on buffer and the other side with the same index
float dirR_buf[][2] = {{0.5, -3.4}, {-0.8, 3.5}, {0.4, -4.1}, {0.0, 0.7}, {0.1, 2.1}, {1.1, 1.1}, 
						{-0.8, 0.1}, {-0.2, 3.1}, {0.0, -2.1}, {-1.0, 1.8}, {-0.1, -2.9}, {0.7, -0.2}, 
						{-0.1, -0.1}, {-0.8, 0.1}, {1.1, -1.0}, {-0.7, 1.1}, {-0.8, 2.5}, {-0.2, 3.1}, 
						{-0.5, -3.9}, {-0.4, 4.1}, {0.4, -1.0}, {-0.6, -2.4}, {0.5, 1.5}, {0.4, 4.4}, 
						{-0.8, 2.5}, {0.8, 3.6}, {-0.5, -3.9}, {-0.5, -3.9}, {0.6, 5.1}, {0.4, -4.1}, 
						{-0.7, -1.0}, {-0.1, -2.9}, {-0.4, 5.0}, {0.4, 4.4}, {-0.8, 2.5}, {0.7, 6.0}, 
						{0.4, -4.1}, {-0.1, -5.2}, {1.1, 4.3}, {0.0, 6.0}, {-0.2, 3.1}, {-1.1, 4.7}, 
						{-0.8, -5.0}, {-0.7, -1.0}, {-0.5, -5.8}, {1.1, 4.6}, {0.5, -3.4}, {-0.4, 4.1},
						{0.8, 3.6}, {0.0, -2.1}, {0.5, -3.4}, {-0.7, 1.1}, {0.8, 3.6}, {0.7, -0.2},
						{-0.4, 5.0}, {-0.5, -3.9}, {0.4, -1.0}, {-0.7, -1.0}, {0.4, -4.1}, {1.1, -1.0},
						{-0.4, 4.1}, {-0.5, -5.8}, {0.7, -0.2}, {-1.0, 1.8}, {-0.8, 2.5}, {-0.5, -3.9},
						{0.1, 2.1}, {-0.8, 2.5}, {-0.8, 0.1}, {-0.7, 1.1}, {-0.5, -3.9}, {0.5, -3.4},
						{-0.6, -2.4}, {0.4, -4.1}, {0.6, 5.1}, {0.6, 5.1}, {-0.8, 3.5}, {-1.0, 1.8}};
float dirR_goal[][2] = {{-14.2, -0.1}, {-14.6, 3.1}, {-15.5, -0.7}, {-15.1, 1.9}, {-14.4, 0.6}, {-14.6, -1.0}, 
						{-14.4, 0.6}, {-15.3, 3.0}, {-16.3, -2.5}, {-15.4, 0.4}, {-15.0, -2.0}, {-16.0, 1.1}, 
						{-14.2, -0.1}, {-16.6, 0.4}, {-15.5, -0.7}, {-17.0, 1.7}, {-15.8, 2.3}, {-14.7, 5.5}, 
						{-16.1, -3.7}, {-14.7, 5.5}, {-15.5, -0.7}, {-14.6, -1.0}, {-17.0, 1.7}, {-16.8, 4.3}, 
						{-16.0, 1.1}, {-16.0, 1.1}, {-15.6, -3.1}, {-14.9, -5.5}, {-15.1, 4.6}, {-15.6, -2.4}, 
						{-16.3, -2.5}, {-15.0, -2.0}, {-14.8, 4.0}, {-17.0, 1.7}, {-14.3, 2.1}, {-15.1, 6.1}, 
						{-16.1, -3.7}, {-14.7, -4.5}, {-15.9, 4.0}, {-14.8, 4.0}, {-15.8, 2.3}, {-15.3, 3.0}, 
						{-14.9, -5.5}, {-15.5, -0.7}, {-16.1, -3.7}, {-14.3, 2.1}, {-15.0, -2.0}, {-14.8, 4.0},
						{-16.0, 1.1}, {-16.3, -2.5}, {-15.0, -2.0}, {-17.0, 1.7}, {-15.0, -2.0}, {-16.0, 1.1},
						{-14.8, 4.0}, {-16.1, -3.7}, {-15.5, -0.7}, {-16.3, -2.5}, {-15.6, -2.4}, {-15.5, -0.7},
						{-14.7, 5.5}, {-15.9, 4.0}, {-16.0, 1.1}, {-15.4, 0.4}, {-14.3, 2.1}, {-16.1, -3.7},
						{-14.4, 0.6}, {-16.0, 1.1}, {-16.6, 0.4}, {-17.0, 1.7}, {-16.1, -3.7}, {-14.2, -0.1},
						{-14.6, -1.0}, {-15.6, -2.4}, {-15.1, 4.6}, {-15.1, 4.6}, {-14.6, 3.1}, {-15.4, 0.4}};
float dirL_buf[][2] = {{-0.4, 4.1}, {-0.6, -2.4}, {0.4, 4.4}, {-0.8, 2.5}, {-0.8, 0.1}, {-0.7, -1.0}, 
						{0.1, 2.1}, {0.4, -1.0}, {0.8, 3.6}, {-0.2, 3.1}, {-0.4, 5.0}, {-0.1, -0.1}, 
						{0.7, 0.5}, {-0.8, -1.7}, {0.6, 2.7}, {0.7, -0.2}, {0.6, -1.7}, {0.0, 0.7}, 
						{0.6, 5.1}, {-0.1, -2.9}, {-0.2, 3.1}, {-0.2, 3.1}, {-0.7, 1.1}, {0.5, -3.4},
						{-0.1, -2.9}, {-0.8, -5.0}, {-0.4, 5.0}, {-0.6, 6.4}, {0.4, -4.1}, {-0.4, 5.0}, 
						{0.7, -0.2}, {-0.2, 3.1}, {-0.5, -3.9}, {-0.5, -3.9}, {0.7, -0.2}, {-0.5, -3.9}, 
						{-0.4, 4.1}, {1.0, 2.0}, {0.4, -4.1}, {-0.5, -5.8}, {-0.6, -2.4}, {-0.5, -3.9}, 
						{0.8, 3.6}, {-0.7, 1.1}, {0.0, 6.0}, {-0.1, -5.2}, {0.4, 4.4}, {-0.1, -2.9},
						{-0.8, -5.0}, {-0.4, 5.0}, {-0.2, 3.1}, {0.7, -0.2}, {-0.8, -5.0}, {-0.1, -0.1},
						{-0.5, -3.9}, {0.0, 6.0}, {0.8, 3.6}, {-0.7, 1.1}, {-0.2, 3.1}, {0.6, 2.7},
						{0.4, -1.0}, {0.4, -4.1}, {-0.1, -0.1}, {-0.2, 3.1}, {-0.1, -0.1}, {-0.6, 6.4}};
float dirL_goal[][2] = {{11.1, 4.0}, {10.0, 0.1}, {11.9, 3.7}, {10.1, 3.2}, {10.8, 2.0}, {12.8, -1.7}, 
						{10.8, 2.0}, {12.4, -0.4}, {11.3, 5.4}, {11.3, 5.4}, {11.4, 2.4}, {10.0, 0.1}, 
						{11.2, -0.4}, {10.4, 1.0}, {12.3, 3.1}, {12.4, -0.4}, {11.3, 0.7}, {11.4, 3.1}, 
						{10.2, 6.7}, {12.8, -1.7}, {11.9, 3.7}, {10.1, 3.2}, {12.1, 2.5}, {10.2, -0.6}, 
						{10.2, -0.6}, {13.2, -5.7}, {10.5, 4.5}, {12.1, 2.5}, {11.8, -1.1}, {12.2, 5.3}, 
						{11.4, 2.4}, {10.5, 4.5}, {10.1, -2.1}, {13.2, -5.7}, {12.4, -0.4}, {11.1, -3.0}, 
						{12.6, 4.5}, {12.3, 3.1}, {11.1, -3.0}, {11.1, -3.0}, {11.2, -0.4}, {12.4, -0.4}, 
						{11.4, 2.4}, {12.1, 2.5}, {10.5, 4.5}, {10.6, -4.0}, {12.9, 5.2}, {10.2, -0.6},
						{13.2, -5.7}, {12.2, 5.3}, {11.9, 3.7}, {12.4, -0.4}, {13.2, -5.7}, {10.0, 0.1},
						{10.1, -2.1}, {10.5, 4.5}, {11.3, 5.4}, {11.8, -1.1}, {11.9, 3.7}, {12.3, 3.1},
						{12.4, -0.4}, {11.8, -1.1}, {10.0, 0.1}, {11.3, 5.4}, {10.0, 0.1}, {12.1, 2.5}};
// used for altering goal if occupied
float buff_sort[42][2] = {{0.0, -6.3}, {0.5, -5.8}, {-0.5, -5.8}, {-0.1, -5.2}, {-0.8, -5.0}, {0.4, -4.1}, 
						{-0.5, -3.9}, {0.5, -3.4}, {-0.1, -2.9}, {-0.6, -2.4}, {0.0, -2.1}, {0.6, -1.7}, 
						{-0.8, -1.7}, {0.4, -1.0}, {-0.7, -1.0}, {1.1, -1.0}, {0.7, -0.2}, {-0.1, -0.1}, 
						{-0.8, 0.1}, {0.7, 0.5}, {0.0, 0.7}, {-0.7, 1.1}, {1.1, 1.1}, {0.5, 1.5}, 
						{-1.0, 1.8}, {1.0, 2.0}, {0.1, 2.1}, {-0.8, 2.5}, {0.6, 2.7}, {-0.2, 3.1}, 
						{-0.8, 3.5}, {0.8, 3.6}, {-0.4, 4.1}, {1.1, 4.3}, {0.4, 4.4}, {1.1, 4.6}, 
						{-1.1, 4.7}, {-0.4, 5.0}, {0.6, 5.1}, {0.0, 6.0}, {0.7, 6.0}, {-0.6, 6.4}};
struct agt_fromL
{
	int hurry;
	int app_time;
	float app_x;
	float app_y;
	bool waitingL;
	bool running_fromL;
	bool waiting_buf;
	float buf_time;
	float buf_x;
	float buf_y;
	bool running_from_buf;
	float goal_x;
	float goal_y;
	int sim_index;
};
struct agt_fromR
{
	int hurry;
	int app_time;
	float app_x;
	float app_y;
	bool waitingR;
	bool running_fromR;
	bool waiting_buf;
	float buf_time;
	float buf_x;
	float buf_y;
	bool running_from_buf;
	float goal_x;
	float goal_y;
	int sim_index;
};
/* init waiting agent structure. */
struct agt_fromL agt_fromL[nLeft];
struct agt_fromR agt_fromR[nRight];

/* Store the goals of the agents. */
std::vector<RVO::Vector2> goals;

void setupScenario(RVO::RVOSimulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(1.0f);

	/* neighborDist, maxNeighbors, timeHorizon, timeHorizonObst, radius, maxSpeed */
	sim->setAgentDefaults(10.0f, 10, 0.1f, 0.1f, 0.3f, 0.09f);

	/* Init with all agents as not waiting. */
	for(int i = 0; i < nLeft; i++) {
		agt_fromL[i].waitingL = false;
	}
	for(int i = 0; i < nRight; i++) {
		agt_fromR[i].waitingR = false;
	}
}

#if RVO_OUTPUT_TIME_AND_POSITIONS
void updateVisualization(RVO::RVOSimulator *sim)
{
	/* Output the current global time. */
	std::cout << sim->getGlobalTime();

	/* Output the current position of all the agents. */
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		std::cout << " " << sim->getAgentPosition(i);
	}

	std::cout << std::endl;
}
#endif

void setPreferredVelocities(RVO::RVOSimulator *sim)
{
	/*
	 * Set the preferred velocity to be a vector of unit magnitude (speed) in the
	 * direction of the goal.
	 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < static_cast<int>(sim->getNumAgents()); ++i) {
		RVO::Vector2 goalVector = goals[i] - sim->getAgentPosition(i);

		if (RVO::absSq(goalVector) > 1.0f) {
			goalVector = RVO::normalize(goalVector);
		}

		sim->setAgentPrefVelocity(i, goalVector);
	}
}

void thesisManipulation(RVO::RVOSimulator *sim)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	/* Agents appear and wait - left */
	for(int i = 0; i < sizeof(seqL)/sizeof(*seqL); i++) {
		if (sim->getGlobalTime() == seqL[i][1]) {
			// check whether the first 30 are all occupied
			bool first_30_occupy = true;
			for (int j = 0; j < 30; j++) {
				bool this_occupy = false;
				for (int k = 0; k < nLeft; k++) {
					if (agt_fromL[k].waitingL && goals[agt_fromL[k].sim_index].x() == posL[j][0] &&
						goals[agt_fromL[k].sim_index].y() == posL[j][1]) {
						this_occupy = true;
					}
				}
				if (!this_occupy) {
					first_30_occupy = false;
					break;
				}
			}
			// if all occupied, the range of ramdom generated position index changes
			if (first_30_occupy) {
				posRange = 18;
				posMin = 30;
			}
			// unoccupied position assigned to a new agent
			int p = rand() % posRange + posMin;
			bool occupy = true;
			while (occupy) {
				occupy = false;
				for (int k = 0; k < nLeft; k++) {
					if (agt_fromL[k].waitingL && goals[agt_fromL[k].sim_index].x() == posL[p][0] &&
						goals[agt_fromL[k].sim_index].y() == posL[p][1]) {
						occupy = true;
					}
				}
				if (occupy == true) { p = rand() % posRange + posMin; }
			}
			// init new agent structure and set state as not-moving
			agt_fromL[i].hurry = seqL[i][0];
			agt_fromL[i].app_time = seqL[i][1];
			agt_fromL[i].app_x = posL[p][0];
			agt_fromL[i].app_y = posL[p][1];
			agt_fromL[i].buf_x = dirL_buf[p][0];
			agt_fromL[i].buf_y = dirL_buf[p][1];
			agt_fromL[i].goal_x = dirL_goal[p][0];
			agt_fromL[i].goal_y = dirL_goal[p][1];
			agt_fromL[i].waitingL = true;
			agt_fromL[i].sim_index = goals.size();
			sim->addAgent(RVO::Vector2(agt_fromL[i].app_x, agt_fromL[i].app_y));
			goals.push_back(RVO::Vector2(agt_fromL[i].app_x, agt_fromL[i].app_y));
		}
	}
	/* Agents appear and wait - right */
	for (int i = 0; i < sizeof(seqR)/sizeof(*seqR); i++) {
		if (sim->getGlobalTime() == seqR[i][1]) {
			// check whether the first 30 are all occupied
			bool first_30_occupy = true;
			for (int j = 0; j < 30; j++) {
				bool this_occupy = false;
				for (int k = 0; k < nRight; k++) {
					if (agt_fromR[k].waitingR && goals[agt_fromR[k].sim_index].x() == posR[j][0] &&
						goals[agt_fromR[k].sim_index].y() == posR[j][1]) {
						this_occupy = true;
					}
				}
				if (!this_occupy) {
					first_30_occupy = false;
					break;
				}
			}
			// if all occupied, the range of ramdom generated position index changes
			if (first_30_occupy) {
				posRange = 48;
				posMin = 30;
			}
			// unoccupied position assigned to a new agent
			int p = rand() % posRange + posMin;
			bool occupy = true;
			while (occupy) {
				occupy = false;
				for (int k = 0; k < nRight; k++) {
					if (agt_fromR[k].waitingR && goals[agt_fromR[k].sim_index].x() == posR[p][0] &&
						goals[agt_fromR[k].sim_index].y() == posR[p][1]) {
						occupy = true;
					}
				}
				if (occupy == true) { p = rand() % posRange + posMin; }
			}
			// init new agent structure and set state as not-moving
			agt_fromR[i].hurry = seqR[i][0];
			agt_fromR[i].app_time = seqR[i][1];
			agt_fromR[i].app_x = posR[p][0];
			agt_fromR[i].app_y = posR[p][1];
			agt_fromR[i].buf_x = dirR_buf[p][0];
			agt_fromR[i].buf_y = dirR_buf[p][1];
			agt_fromR[i].goal_x = dirR_goal[p][0];
			agt_fromR[i].goal_y = dirR_goal[p][1];
			agt_fromR[i].waitingR = true;
			agt_fromR[i].sim_index = goals.size();
			sim->addAgent(RVO::Vector2(agt_fromR[i].app_x, agt_fromR[i].app_y));
			goals.push_back(RVO::Vector2(agt_fromR[i].app_x, agt_fromR[i].app_y));
		}
	}

	/* Check running light - from left */
	for (int i = 0; i < nLeft; i++) {
		if (agt_fromL[i].waitingL == true && (
		(sim->getGlobalTime() < 460) ||
		(sim->getGlobalTime() > 1500))) {
			int waiting_time = sim->getGlobalTime() - agt_fromL[i].app_time;
			// reset initial waiting time after 30-50 vehicle stream
			if (sim->getGlobalTime() > 1500 && agt_fromR[i].app_time > 390 && agt_fromR[i].app_time < 1500) {
				waiting_time = sim->getGlobalTime() - 1470;
			}
			if (illegalR + illegalL >= getThreLeft(waiting_time, agt_fromL[i].hurry)) {
				goals[agt_fromL[i].sim_index] = RVO::Vector2(agt_fromL[i].buf_x, agt_fromL[i].buf_y);
				illegalL++;
				agt_fromL[i].waitingL = false;
				agt_fromL[i].running_fromL = true;
			}
		}
	}
	/* Check running light - from right */
	for (int i = 0; i < nRight; i++) {
		if (agt_fromR[i].waitingR && (
		(sim->getGlobalTime() < 210) ||
		(sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810) ||
		(sim->getGlobalTime() > 1035))) {
			int waiting_time = sim->getGlobalTime() - agt_fromR[i].app_time;
			// adjust waiting time if have met the vehicle flows
			if (sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810 && agt_fromR[i].app_time < 300) {
				waiting_time -= 150;
			}
			if (sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810 && agt_fromR[i].app_time > 300 && agt_fromR[i].app_time < 600) {
				waiting_time -= 90;
			}
			if (sim->getGlobalTime() > 1035 && agt_fromR[i].app_time > 690 && agt_fromR[i].app_time < 1035) {
				waiting_time = sim->getGlobalTime() - 1020; // reset initial waiting time after 30-50 vehicle stream
			}
			// look up the state table about tolerable illegals in view
			if (illegalR + illegalL >= getThreRight(waiting_time, agt_fromR[i].hurry)) {
				goals[agt_fromR[i].sim_index] = RVO::Vector2(agt_fromR[i].buf_x, agt_fromR[i].buf_y);
				illegalR++;
				agt_fromR[i].waitingR = false;
				agt_fromR[i].running_fromR = true;
			}
		}
	}
	/* Goals at buffer avoid overlapping */
	for (size_t i = 0; i < sim->getNumAgents(); i++) {
		if (-1.8f <sim->getAgentPosition(i).x() < 1.8f) {
			for (size_t j = 0; j < sim->getNumAgents(); j++) {
				if (i != j && goals[i].x() == sim->getAgentPosition(j).x() && goals[i].y() == sim->getAgentPosition(j).y()) {
					for (int m = 0; m < 42; m++) {
						if (goals[i].x() == buff_sort[m][0] && goals[i].y() == buff_sort[m][1]) {
							// match goals[index] and sim_index to reset goal
							// sim_index either find in agt_fromL or agt_fromR
							for (int n = 0; n < nLeft; n++) {
								if (agt_fromL[n].sim_index == i) {
									if (rand() % 2 == 0) {
										if (m<41) {
											agt_fromL[n].buf_x = buff_sort[m+1][0];
											agt_fromL[n].buf_y = buff_sort[m+1][1];
										} else {
											agt_fromL[n].buf_x = buff_sort[m-1][0];
											agt_fromL[n].buf_y = buff_sort[m-1][1];
										}
									} else {
										if (m>0) {
											agt_fromL[n].buf_x = buff_sort[m-1][0];
											agt_fromL[n].buf_y = buff_sort[m-1][1];
										} else {
											agt_fromL[n].buf_x = buff_sort[m+1][0];
											agt_fromL[n].buf_y = buff_sort[m+1][1];
										}
									}
									goals[i] = RVO::Vector2(agt_fromL[n].buf_x, agt_fromL[n].buf_y);
								}
							}
							for (int n = 0; n < nRight; n++) {
								if (agt_fromR[n].sim_index == i) {
									if (rand() % 2 == 0) {
										if (m<41) {
											agt_fromR[n].buf_x = buff_sort[m+1][0];
											agt_fromR[n].buf_y = buff_sort[m+1][1];
										} else {
											agt_fromR[n].buf_x = buff_sort[m-1][0];
											agt_fromR[n].buf_y = buff_sort[m-1][1];
										}
									} else {
										if (m>0) {
											agt_fromR[n].buf_x = buff_sort[m-1][0];
											agt_fromR[n].buf_y = buff_sort[m-1][1];
										} else {
											agt_fromR[n].buf_x = buff_sort[m+1][0];
											agt_fromR[n].buf_y = buff_sort[m+1][1];
										}
									}
									goals[i] = RVO::Vector2(agt_fromR[n].buf_x, agt_fromR[n].buf_y);
								}
							}
						}
					}
				}
			}
		}
	}

	/* Reach buffer reset state- from left */
	for (int i = 0; i < nLeft; i++) {
		if (agt_fromL[i].running_fromL) {
			if (RVO::absSq(sim->getAgentPosition(i) - goals[i]) <= sim->getAgentRadius(i) * sim->getAgentRadius(i)) {
				agt_fromL[i].running_fromL = false;
				agt_fromL[i].waiting_buf = true;
				agt_fromL[i].buf_time = sim->getGlobalTime();
				illegalL--;
			}
		}
	}
	//* Reach buffer reset state- from right */
	for (int i = 0; i < nRight; i++) {
		if (agt_fromR[i].running_fromR) {
			if (RVO::absSq(sim->getAgentPosition(i) - goals[i]) <= sim->getAgentRadius(i) * sim->getAgentRadius(i)) {
				agt_fromR[i].running_fromR = false;
				agt_fromR[i].waiting_buf = true;
				agt_fromR[i].buf_time = sim->getGlobalTime();
				illegalR--;
			}
		}
	}

	/* Check running light - buffer to left */
	for (int i = 0; i < nRight; i++) {
		if ((agt_fromR[i].waiting_buf == true) && (
		(sim->getGlobalTime() < 540) ||
		(sim->getGlobalTime() > 1500))) {
			int waiting_time = sim->getGlobalTime() - agt_fromR[i].buf_time;
			// reset initial waiting time after 30-50 vehicle stream
			if (sim->getGlobalTime() > 1500 && agt_fromR[i].buf_time > 390 && agt_fromR[i].buf_time < 1500) {
				waiting_time = sim->getGlobalTime() - 1485;
			}
			if (illegalL >= getThreLeft(waiting_time, agt_fromR[i].hurry)) {
				goals[agt_fromR[i].sim_index] = RVO::Vector2(agt_fromR[i].goal_x, agt_fromR[i].goal_y);
				illegalL++;
				agt_fromR[i].waiting_buf = false;
				agt_fromR[i].running_from_buf = true;
			}
		}
	}
	/* Check running light - buffer to right */
	for (int i = 0; i < nLeft; i++) {
		if ((agt_fromL[i].waiting_buf == true) && (
		(sim->getGlobalTime() < 210) ||
		(sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810) ||
		(sim->getGlobalTime() > 1035))) {
			int waiting_time = sim->getGlobalTime() - agt_fromL[i].buf_time;
			// adjust waiting time if have met the vehicle flows
			if (sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810 && agt_fromL[i].buf_time < 300) {
				waiting_time -= 150;
			}
			if (sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810 && agt_fromL[i].buf_time > 300 && agt_fromL[i].buf_time < 600) {
				waiting_time -= 90;
			}
			if (sim->getGlobalTime() > 1035 && agt_fromL[i].buf_time > 690 && agt_fromL[i].buf_time < 1035) {
				waiting_time = sim->getGlobalTime() - 1020; // reset initial waiting time after 30-50 vehicle stream
			}
			if (illegalR >= getThreRight(waiting_time, agt_fromL[i].hurry)) {
				goals[agt_fromL[i].sim_index] = RVO::Vector2(agt_fromL[i].goal_x, agt_fromL[i].goal_y);
				illegalR++;
				agt_fromL[i].waiting_buf = false;
				agt_fromL[i].running_from_buf = true;
			}
		}
	}

	/* Light green - from left */
	for (int i = 0; i < nLeft; i++) {
		// the agents waiting at the left side
		if ((agt_fromL[i].waitingL || agt_fromL[i].running_fromL) && (sim->getGlobalTime() > 1800)) {
			goals[agt_fromL[i].sim_index] = RVO::Vector2(agt_fromL[i].goal_x, agt_fromL[i].goal_y);
			agt_fromL[i].waitingL = false;
		}
		// the agents waiting at buffer to go right
		if ((agt_fromL[i].waiting_buf) && (sim->getGlobalTime() > 1800)) {
			goals[agt_fromL[i].sim_index] = RVO::Vector2(agt_fromL[i].goal_x, agt_fromL[i].goal_y);
			agt_fromL[i].waiting_buf = false;
		}
	}
	/* Light green - from right */
	for (int i = 0; i < nRight; i++) {
		// the agents waiting at the right side
		if ((agt_fromR[i].waitingR || agt_fromR[i].running_fromR) && (sim->getGlobalTime() > 1800)) {
			goals[agt_fromR[i].sim_index] = RVO::Vector2(agt_fromR[i].goal_x, agt_fromR[i].goal_y);
			agt_fromR[i].waitingR = false;
		}
		// the agents waiting at buffer to go left
		if ((agt_fromR[i].waiting_buf) && (sim->getGlobalTime() > 1800)) {
			goals[agt_fromR[i].sim_index] = RVO::Vector2(agt_fromR[i].goal_x, agt_fromR[i].goal_y);
			agt_fromR[i].waiting_buf = false;
		}
	}

	/* Remove crossed agents - from left */
	for (int i = 0; i < nLeft; i++) {
		if (RVO::absSq(sim->getAgentPosition(agt_fromL[i].sim_index) - goals[agt_fromL[i].sim_index]) <= 0.09f &&
			goals[agt_fromL[i].sim_index].x() == agt_fromL[i].goal_x && 
			goals[agt_fromL[i].sim_index].y() == agt_fromL[i].goal_y) {
				sim->setAgentPosition(agt_fromL[i].sim_index, RVO::Vector2(100.0f, 100.0f));
				sim->setAgentMaxSpeed(agt_fromL[i].sim_index, 0.0f);
		}
	}
	/* Remove crossed agents - from right */
	for (int i = 0; i < nRight; i++) {
		if (RVO::absSq(sim->getAgentPosition(agt_fromR[i].sim_index) - goals[agt_fromR[i].sim_index]) <= 0.09f &&
			goals[agt_fromR[i].sim_index].x() == agt_fromR[i].goal_x && 
			goals[agt_fromR[i].sim_index].y() == agt_fromR[i].goal_y) {
				sim->setAgentPosition(agt_fromR[i].sim_index, RVO::Vector2(100.0f, 100.0f));
				sim->setAgentMaxSpeed(agt_fromR[i].sim_index, 0.0f);
		}
	}
	
}

/* Max amount of illegal crossings an agent is able to resist. */
int getThreLeft(int t, int hurry) {
	if (t == 30 * 0) {
		switch (hurry)
		{
			case 0: return 50;
			case 1: return 35;
			case 2: return 25;
			case 3: return 16;
			case 4: return 10;
			case 5: return 7;
			case 6: return 4;
			case 7: return 2;
			case 8: return 1;
			case 9: return 1;
			default: return 0;
		}
	} else if (t == 30 * 1) {
		switch (hurry)
		{
			case 0: return 35;
			case 1: return 22;
			case 2: return 16;
			case 3: return 9;
			case 4: return 6;
			case 5: return 4;
			case 6: return 2;
			case 7: return 1;
			case 8: return 1;
			default: return 0;
		}
	} else if (t == 30 * 2) {
		switch (hurry)
		{
			case 0: return 25;
			case 1: return 17;
			case 2: return 11;
			case 3: return 7;
			case 4: return 4;
			case 5: return 2;
			case 6: return 1;
			default: return 0;
		} 
	} else if (t == 30 * 3) {
		switch (hurry)
		{
			case 0: return 22;
			case 1: return 14;
			case 2: return 9;
			case 3: return 6;
			case 4: return 2;
			case 5: return 1;
			default: return 0;
		}
	} else if (t == 30 * 4) {
		switch (hurry)
		{
			case 0: return 19;
			case 1: return 12;
			case 2: return 7;
			case 3: return 5;
			case 4: return 2;
			default: return 0;
		}
	} else if (t == 30 * 5) {
		switch (hurry)
		{
			case 0: return 15;
			case 1: return 9;
			case 2: return 5;
			case 3: return 4;
			case 4: return 1;
			default: return 0;
		}
	} else if (t == 30 * 6) {
		switch (hurry)
		{
			case 0: return 12;
			case 1: return 8;
			case 2: return 3;
			case 3: return 2;
			default: return 0;
		}
	} else if (t == 30 * 7) {
		switch (hurry)
		{
			case 0: return 10;
			case 1: return 7;
			case 2: return 2;
			case 3: return 1;
			default: return 0;
		}
	}
	return 1000;
}

int getThreRight(int t, int hurry) {
	if (t == 30 * 0) {
		switch (hurry)
		{
			case 0: return 40;
			case 1: return 25;
			case 2: return 20;
			case 3: return 13;
			case 4: return 8;
			case 5: return 5;
			case 6: return 3;
			case 7: return 2;
			case 8: return 1;
			default: return 0;
		}
	} else if (t == 30 * 1) {
		switch (hurry)
		{
			case 0: return 30;
			case 1: return 18;
			case 2: return 13;
			case 3: return 7;
			case 4: return 5;
			case 5: return 3;
			case 6: return 1;
			default: return 0;
		}
	} else if (t == 30 * 2) {
		switch (hurry)
		{
			case 0: return 20;
			case 1: return 13;
			case 2: return 9;
			case 3: return 6;
			case 4: return 3;
			case 5: return 2;
			default: return 0;
		} 
	} else if (t == 30 * 3) {
		switch (hurry)
		{
			case 0: return 18;
			case 1: return 11;
			case 2: return 7;
			case 3: return 4;
			case 4: return 1;
			default: return 0;
		}
	} else if (t == 30 * 4) {
		switch (hurry)
		{
			case 0: return 15;
			case 1: return 9;
			case 2: return 5;
			case 3: return 2;
			default: return 0;
		}
	} else if (t == 30 * 5) {
		switch (hurry)
		{
			case 0: return 13;
			case 1: return 7;
			case 2: return 4;
			case 3: return 1;
			default: return 0;
		}
	} else if (t == 30 * 6) {
		switch (hurry)
		{
			case 0: return 10;
			case 1: return 6;
			case 2: return 2;
			default: return 0;
		}
	} else if (t == 30 * 7) {
		switch (hurry)
		{
			case 0: return 9;
			case 1: return 4;
			case 2: return 1;
			default: return 0;
		}
	}
	return 1000;
}

bool timeUp(RVO::RVOSimulator *sim)
{
	return sim->getGlobalTime() > 2550.0f;
}

int main()
{
	/* Create a new simulator instance. */
	RVO::RVOSimulator *sim = new RVO::RVOSimulator();

	/* Set up the scenario. */
	setupScenario(sim);

	/* Perform (and manipulate) the simulation. */
	do {
#if RVO_OUTPUT_TIME_AND_POSITIONS
		updateVisualization(sim);
#endif
		setPreferredVelocities(sim);
		thesisManipulation(sim);
		sim->doStep();
	}
	while (!timeUp(sim));

	delete sim;

	return 0;
}
