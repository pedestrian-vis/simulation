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

#define nLeft 10
int getThreLeft(int t, int hurry);

/* Store illegal crossing parameter. */
int seqL[nLeft][2] = {{5,30}, {3,60}, {1,90}, {10,120}, {3,150}, {10,180}, {6,210}, {1,240}, {10,270}, {4,300}};
int illegal = 0;
struct waitL
{
	int hurry;
	int app_time;
	float app_x;
	float app_y;
	bool waiting;
};
float posL[][2] = {{-14.6, 3.1}, {-15.3, -4.0}, {-15.1, 4.6}, {-15.1, 1.9}, {-15.0, -2.0}, {-15.5, -0.7},
									{-14.4, 0.6}, {-14.3, -3.2}, {-15.3, 3.0}, {-14.3, 2.1}, {-14.8, 4.0}, {-14.6, -1.0},
									{-15.4, 0.4}, {-14.2, -0.1}, {-16.0, 1.1}, {-16.3, -2.5}, {-15.6, -3.1}, {-16.5, -3.1},
									{-14.7, 5.5}, {-16.1, -3.7}, {-15.8, 2.3}, {-15.9, 4.0}, {-16.2, -0.1}, {-14.9, -5.5},
									{-14.7, -4.5}, {-14.9, -7.1}, {-15.1, 7.2}, {-14.2, 5.9}, {-14.4, -6.7}, {-15.1, 6.1},
									{-16.6, 0.4}, {-16.8, 4.3}, {-15.1, -8.0}, {-16.0, -6.2}, {-15.6, -2.4}, {-15.5, -5.4},
									{-16.7, 6.0}, {-15.6, 7.9}, {-16.7, -6.4}, {-14.4, -8.4}, {-16.8, -4.5}, {-16.4, -5.3},
									{-16.3, 7.1}, {-17.0, 1.7}, {-14.7, 8.0}, {-16.2, -8.1}, {-17.4, 5.7}, {-17.1, -5.2}};


/* Store the goals of the agents. */
std::vector<RVO::Vector2> goals;

void setupScenario(RVO::RVOSimulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(1.0f);

	/* neighborDist, maxNeighbors, timeHorizon, timeHorizonObst, radius, maxSpeed */
	sim->setAgentDefaults(10.0f, 5, 0.1f, 0.1f, 0.3f, 0.09f);
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
	// if (sim->getGlobalTime() == 30.0f) {
	// 	goals[1] = RVO::Vector2(9.9f, -5.0f);
	// }
	struct waitL wait_left[nLeft];
	for(int i = 0; i < nLeft; i++) {
		wait_left[i].waiting = false;
	}
	for(int i = 0; i < sizeof(seqL)/sizeof(*seqL); i++) {
		// agents appears and waiting
		if (sim->getGlobalTime() == seqL[i][1]) {
			int p = 0;
			bool occupy = true;
			// unoccupied position with highest priority is assigned to new agent
			while (occupy == true && p < sizeof(posL)/sizeof(*posL)) {
				occupy = false;
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == posL[p][0] && sim->getAgentPosition(k).y() == posL[p][1]) {
						occupy = true;
					}
				}
				if (occupy == true) { p++; }
			}
			sim->addAgent(RVO::Vector2(posL[p][0], posL[p][1]));
			goals.push_back(RVO::Vector2(posL[p][0], posL[p][1]));
			wait_left[i].hurry = seqL[i][0];
			wait_left[i].app_time = seqL[i][1];
			wait_left[i].app_x = posL[p][0];
			wait_left[i].app_y = posL[p][1];
			wait_left[i].waiting = true;
		}
	}
	// chech waiting agents whether run the light
	for (int i = 0; i < nLeft; i++) {
		if (wait_left[i].waiting == true) {
			if (illegal >= getThreLeft(sim->getGlobalTime() - wait_left[i].app_time, wait_left[i].hurry)) {
				for (size_t k = 0; k < sim->getNumAgents(); k++) {
					if (sim->getAgentPosition(k).x() == wait_left[i].app_x && sim->getAgentPosition(k).y() == wait_left[i].app_y) {
						goals[k] = RVO::Vector2(9.9f, 0.0f);
					}
				}
				illegal++;
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
