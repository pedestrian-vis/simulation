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

#define nLeft 15
#define nRight 30
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
					{2, 600}, {5, 630}, {1, 660}, {5, 690}, {1, 720}, {5, 750}, {4, 780}, {4, 810}, {8, 840}, {10, 870}};
int seqL[nLeft][2] = {{8, 45}, {4, 105}, {1, 165}, {10, 225}, {8, 285}, {1, 345}, {9, 405}, {4, 465}, {5, 525}, {2, 585},
					{2, 645}, {1, 705}, {4, 765}, {1, 825}, {7, 885}};
// store the priority of initial waiting positions
float posL[][2] = {{-14.6, 3.1}, {-15.3, -4.0}, {-15.1, 4.6}, {-15.1, 1.9}, {-15.0, -2.0}, {-15.5, -0.7},
					{-14.4, 0.6}, {-14.3, -3.2}, {-15.3, 3.0}, {-14.3, 2.1}, {-14.8, 4.0}, {-14.6, -1.0},
					{-15.4, 0.4}, {-14.2, -0.1}, {-16.0, 1.1}, {-16.3, -2.5}, {-15.6, -3.1}, {-16.5, -3.1},
					{-14.7, 5.5}, {-16.1, -3.7}, {-15.8, 2.3}, {-15.9, 4.0}, {-16.2, -0.1}, {-14.9, -5.5},
					{-14.7, -4.5}, {-14.9, -7.1}, {-15.1, 7.2}, {-14.2, 5.9}, {-14.4, -6.7}, {-15.1, 6.1},
					{-16.6, 0.4}, {-16.8, 4.3}, {-15.1, -8.0}, {-16.0, -6.2}, {-15.6, -2.4}, {-15.5, -5.4},
					{-16.7, 6.0}, {-15.6, 7.9}, {-16.7, -6.4}, {-14.4, -8.4}, {-16.8, -4.5}, {-16.4, -5.3},
					{-16.3, 7.1}, {-17.0, 1.7}, {-14.7, 8.0}, {-16.2, -8.1}, {-17.4, 5.7}, {-17.1, -5.2}};
float posR[][2] = {{10.4, -3.1}, {11.1, 4.0}, {10.9, -4.6}, {10.9, -1.9}, {10.8, 2.0}, {11.3, 0.7},
					{10.2, -0.6}, {10.1, 3.2}, {11.1, -3.0}, {10.1, -2.1}, {10.6, -4.0}, {10.4, 1.0},
					{11.2, -0.4}, {10.0, 0.1}, {11.8, -1.1}, {12.1, 2.5}, {11.4, 3.1}, {12.3, 3.1},
					{10.5, -5.5}, {11.9, 3.7}, {11.6, -2.3}, {11.7, -4.0}, {12.0, 0.1}, {10.7, 5.5},
					{10.5, 4.5}, {10.7, 7.1}, {10.9, -7.2}, {10.0, -5.9}, {10.2, 6.7}, {10.9, -6.1},
					{12.4, -0.4}, {12.6, -4.3}, {10.9, 8.0}, {11.8, 6.2}, {11.4, 2.4}, {11.3, 5.4},
					{12.5, -6.0}, {11.4, -7.9}, {12.5, 6.4}, {10.2, 8.4}, {12.6, 4.5}, {12.2, 5.3},
					{12.1, -7.1}, {12.8, -1.7}, {10.5, -8.0}, {12.0, 8.1}, {13.2, -5.7}, {12.9, 5.2}};
// storing the goals on buffer and the other side with the same index
float dirR_buf[][2] = {{0.5, -3.4}, {-0.8, 3.5}, {0.4, -4.1}, {0.0, 0.7}, {0.1, 2.1}, {1.1, 1.1}, 
						{-0.8, 0.1}, {-0.2, 3.1}, {0.0, -2.1}, {-1.0, 1.8}, {-0.1, -2.9}, {0.7, -0.2}, 
						{-0.1, -0.1}, {-0.8, 0.1}, {1.1, -1.0}, {-0.7, 1.1}, {-0.8, 2.5}, {-0.2, 3.1}, 
						{-0.5, -3.9}, {-0.4, 4.1}, {0.4, -1.0}, {-0.6, -2.4}, {0.5, 1.5}, {0.4, 4.4}, 
						{-0.8, 2.5}, {0.8, 3.6}, {-0.5, -3.9}, {-0.5, -3.9}, {0.6, 5.1}, {0.4, -4.1}, 
						{-0.7, -1.0}, {-0.1, -2.9}, {-0.4, 5.0}, {0.4, 4.4}, {-0.8, 2.5}, {0.7, 6.0}, 
						{0.4, -4.1}, {-0.1, -5.2}, {1.1, 4.3}, {0.0, 6.0}, {-0.2, 3.1}, {-1.1, 4.7}, 
						{-0.8, -5.0}, {-0.7, -1.0}, {-0.5, -5.8}, {1.1, 4.6}, {0.5, -3.4}, {-0.4, 4.1}};
float dirR_goal[][2] = {{-14.2, -0.1}, {-14.6, 3.1}, {-15.5, -0.7}, {-15.1, 1.9}, {-14.4, 0.6}, {-14.6, -1.0}, 
						{-14.4, 0.6}, {-15.3, 3.0}, {-16.3, -2.5}, {-15.4, 0.4}, {-15.0, -2.0}, {-16.0, 1.1}, 
						{-14.2, -0.1}, {-16.6, 0.4}, {-15.5, -0.7}, {-17.0, 1.7}, {-15.8, 2.3}, {-14.7, 5.5}, 
						{-16.1, -3.7}, {-14.7, 5.5}, {-15.5, -0.7}, {-14.6, -1.0}, {-17.0, 1.7}, {-16.8, 4.3}, 
						{-16.0, 1.1}, {-16.0, 1.1}, {-15.6, -3.1}, {-14.9, -5.5}, {-15.1, 4.6}, {-15.6, -2.4}, 
						{-16.3, -2.5}, {-15.0, -2.0}, {-14.8, 4.0}, {-17.0, 1.7}, {-14.3, 2.1}, {-15.1, 6.1}, 
						{-16.1, -3.7}, {-14.7, -4.5}, {-15.9, 4.0}, {-14.8, 4.0}, {-15.8, 2.3}, {-15.3, 3.0}, 
						{-14.9, -5.5}, {-15.5, -0.7}, {-16.1, -3.7}, {-14.3, 2.1}, {-15.0, -2.0}, {-14.8, 4.0}};
float dirL_buf[][2] = {{-0.4, 4.1}, {-0.6, -2.4}, {0.4, 4.4}, {-0.8, 2.5}, {-0.8, 0.1}, {-0.7, -1.0}, 
						{0.1, 2.1}, {0.4, -1.0}, {0.8, 3.6}, {-0.2, 3.1}, {-0.4, 5.0}, {-0.1, -0.1}, 
						{0.7, 0.5}, {-0.8, -1.7}, {0.6, 2.7}, {0.7, -0.2}, {0.6, -1.7}, {0.0, 0.7}, 
						{0.6, 5.1}, {-0.1, -2.9}, {-0.2, 3.1}, {-0.2, 3.1}, {-0.7, 1.1}, {0.5, -3.4},
						{-0.1, -2.9}, {-0.8, -5.0}, {-0.4, 5.0}, {-0.6, 6.4}, {0.4, -4.1}, {-0.4, 5.0}, 
						{0.7, -0.2}, {-0.2, 3.1}, {-0.5, -3.9}, {-0.5, -3.9}, {0.7, -0.2}, {-0.5, -3.9}, 
						{-0.4, 4.1}, {1.0, 2.0}, {0.4, -4.1}, {-0.5, -5.8}, {-0.6, -2.4}, {-0.5, -3.9}, 
						{0.8, 3.6}, {-0.7, 1.1}, {0.0, 6.0}, {-0.1, -5.2}, {0.4, 4.4}, {-0.1, -2.9}};
float dirL_goal[][2] = {{11.1, 4.0}, {10.0, 0.1}, {11.9, 3.7}, {10.1, 3.2}, {10.8, 2.0}, {12.8, -1.7}, 
						{10.8, 2.0}, {12.4, -0.4}, {11.3, 5.4}, {11.3, 5.4}, {11.4, 2.4}, {10.0, 0.1}, 
						{11.2, -0.4}, {10.4, 1.0}, {12.3, 3.1}, {12.4, -0.4}, {11.3, 0.7}, {11.4, 3.1}, 
						{10.2, 6.7}, {12.8, -1.7}, {11.9, 3.7}, {10.1, 3.2}, {12.1, 2.5}, {10.2, -0.6}, 
						{10.2, -0.6}, {13.2, -5.7}, {10.5, 4.5}, {12.1, 2.5}, {11.8, -1.1}, {12.2, 5.3}, 
						{11.4, 2.4}, {10.5, 4.5}, {10.1, -2.1}, {13.2, -5.7}, {12.4, -0.4}, {11.1, -3.0}, 
						{12.6, 4.5}, {12.3, 3.1}, {11.1, -3.0}, {11.1, -3.0}, {11.2, -0.4}, {12.4, -0.4}, 
						{11.4, 2.4}, {12.1, 2.5}, {10.5, 4.5}, {10.6, -4.0}, {12.9, 5.2}, {10.2, -0.6}};
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
	// agents appears and waiting - left
	for(int i = 0; i < sizeof(seqL)/sizeof(*seqL); i++) {
		if (sim->getGlobalTime() == seqL[i][1]) {
			// check whether the first 30 are all occupied
			bool first_30_occupy = true;
			for (int j = 0; j < 30; j++) {
				bool this_occupy = false;
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == posL[j][0] && sim->getAgentPosition(k).y() == posL[j][1]) {
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
			while (occupy == true && p < sizeof(posL)/sizeof(*posL)) {
				occupy = false;
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == posL[p][0] && sim->getAgentPosition(k).y() == posL[p][1]) {
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
			agt_fromL[i].waitingL = true;
			agt_fromL[i].sim_index = goals.size() - 1;
			sim->addAgent(RVO::Vector2(agt_fromL[i].app_x, agt_fromL[i].app_y));
			goals.push_back(RVO::Vector2(agt_fromL[i].app_x, agt_fromL[i].app_y));
		}
	}
	// agents appears and waiting - right
	for (int i = 0; i < sizeof(seqR)/sizeof(*seqR); i++) {
		if (sim->getGlobalTime() == seqR[i][1]) {
			// check whether the first 30 are all occupied
			bool first_30_occupy = true;
			for (int j = 0; j < 30; j++) {
				bool this_occupy = false;
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == posL[j][0] && sim->getAgentPosition(k).y() == posL[j][1]) {
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
			while (occupy == true && p < sizeof(posR)/sizeof(*posR)) {
				occupy = false;
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == posR[p][0] && sim->getAgentPosition(k).y() == posR[p][1]) {
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
			agt_fromR[i].waitingR = true;
			agt_fromR[i].sim_index = goals.size() - 1;
			sim->addAgent(RVO::Vector2(agt_fromR[i].app_x, agt_fromR[i].app_y));
			goals.push_back(RVO::Vector2(agt_fromR[i].app_x, agt_fromR[i].app_y));
		}
	}

	// check waiting agents whether run the light - from left
	for (int i = 0; i < nLeft; i++) {
		if ((agt_fromL[i].waitingL == true && (sim->getGlobalTime() < 460))) {
			if (illegalR + illegalL >= getThreLeft(sim->getGlobalTime() - agt_fromL[i].app_time, agt_fromL[i].hurry)) {
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == agt_fromL[i].app_x && sim->getAgentPosition(k).y() == agt_fromL[i].app_y) {
						// goal at buffer and final set via the same index
						for(int j = 0; j < sizeof(posL)/sizeof(*posL); j++) {
							if (agt_fromL[i].app_x == posL[j][0] && agt_fromL[i].app_y == posL[j][1]) {
								agt_fromL[i].buf_x = dirL_buf[j][0];
								agt_fromL[i].buf_y = dirL_buf[j][1];
								agt_fromL[i].goal_x = dirL_goal[j][0];
								agt_fromL[i].goal_y = dirL_goal[j][1];
								goals[k] = RVO::Vector2(agt_fromL[i].buf_x, agt_fromL[i].buf_y);
							}
						}
					}
				}
				illegalL++;
				agt_fromL[i].waitingL = false;
				agt_fromL[i].running_fromL = true;
			}
		}
	}
	// check waiting agents whether run the light - from right
	for (int i = 0; i < nRight; i++) {
		if ((agt_fromR[i].waitingR == true) && (
			(sim->getGlobalTime() < 210) ||
			(sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810))) {
			if (illegalR + illegalL >= getThreRight(sim->getGlobalTime() - agt_fromR[i].app_time, agt_fromR[i].hurry)) {
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == agt_fromR[i].app_x && sim->getAgentPosition(k).y() == agt_fromR[i].app_y) {
						// goal at buffer and final set via the same index
						for(int j = 0; j < sizeof(posR)/sizeof(*posR); j++) {
							if (agt_fromR[i].app_x == posR[j][0] && agt_fromR[i].app_y == posR[j][1]) {
								agt_fromR[i].buf_x = dirR_buf[j][0];
								agt_fromR[i].buf_y = dirR_buf[j][1];
								agt_fromR[i].goal_x = dirR_goal[j][0];
								agt_fromR[i].goal_y = dirR_goal[j][1];
								goals[k] = RVO::Vector2(agt_fromR[i].buf_x, agt_fromR[i].buf_y);
							}
						}
					}
				}
				illegalR++;
				agt_fromR[i].waitingR = false;
				agt_fromR[i].running_fromR = true;
			}
		}
	}
	// adjust goal at buffer if the original goal is occupied
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

	// if reaches buffer, reset state as waiting - from left
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
	// if reaches buffer, reset state as waiting - from right
	for (int i = 0; i < nRight; i++) {
		if (agt_fromR[i].running_fromR) {
			if (RVO::absSq(sim->getAgentPosition(i) - goals[i]) <= sim->getAgentRadius(i) * sim->getAgentRadius(i)) {
				agt_fromR[i].running_fromR = false;
				agt_fromR[i].waiting_buf = true;
				agt_fromL[i].buf_time = sim->getGlobalTime();
				illegalR--;
			}
		}
	}

	// check buffer waiting agents whether run the light - to left
	for (int i = 0; i < nRight; i++) {
		if ((agt_fromR[i].waiting_buf == true) && (
			(sim->getGlobalTime() < 540))) {
			if (illegalL >= getThreLeft(sim->getGlobalTime() - agt_fromR[i].buf_time, agt_fromR[i].hurry)) {
				// find corresponding k and reset goal
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (goals[k].x() == agt_fromR[i].buf_x && goals[k].y() == agt_fromR[i].buf_y) {
						goals[k] = RVO::Vector2(agt_fromR[i].goal_x, agt_fromR[i].goal_y);
					}
				}
				illegalL++;
				agt_fromR[i].waiting_buf = false;
				agt_fromR[i].running_from_buf = true;
			}
		}
	}
	// check buffer waiting agents whether run the light - to right
	for (int i = 0; i < nLeft; i++) {
		if ((agt_fromL[i].waiting_buf == true) && (
			(sim->getGlobalTime() < 210) ||
			(sim->getGlobalTime() > 540 && sim->getGlobalTime() < 810))) {
			if (illegalR >= getThreRight(sim->getGlobalTime() - agt_fromL[i].buf_time, agt_fromL[i].hurry)) {
				// find corresponding k and reset goal
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (goals[k].x() == agt_fromL[i].buf_x && goals[k].y() == agt_fromL[i].buf_y) {
						goals[k] = RVO::Vector2(agt_fromL[i].goal_x, agt_fromL[i].goal_y);
					}
				}
				illegalR++;
				agt_fromL[i].waiting_buf = false;
				agt_fromL[i].running_from_buf = true;
			}
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


bool reachedGoal(RVO::RVOSimulator *sim)
{
	/* Check if all agents have reached their goals. */
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		if (RVO::absSq(sim->getAgentPosition(i) - goals[i]) > sim->getAgentRadius(i) * sim->getAgentRadius(i)) {
			return false;
		}
	}

	return true;
}

bool timeUp(RVO::RVOSimulator *sim)
{
	return sim->getGlobalTime() > 1000.0f;
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
	// while (!reachedGoal(sim));
	/* Force stop if agents never reach their goals. */
	while (!timeUp(sim));

	delete sim;

	return 0;
}
