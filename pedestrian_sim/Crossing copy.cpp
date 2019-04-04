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

/* Store the parameters for illegal crossing judgement. */
int illeal = 0;
std::vector<list> seqL; // [[hurrylevel, appearingtime]...] [[1,2]]
// std::vector<list> seqR;
// std::vector<list> posL; // [[[x,y],[x_goalBuffer,y_goalBuffer],[x_goal,y_goal]...]
// std::vector<list> posR;
// std::vector<list> posBuffer;
// std::vector<int> waitL; // [[hurrylevel, appearingtime]...]
// std::vector<int> waitR;
// std::vector<int> waitBuffer; // [[hurrylevel, reachbuffertime]...]

/* Store the goals of the agents. */
std::vector<RVO::Vector2> goals;

void setupScenario(RVO::RVOSimulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(0.25f);

	/* Specify the default parameters for agents that are subsequently added. */
	sim->setAgentDefaults(10.0f, 5, 0.1f, 0.1f, 0.3f, 0.09f);

	/*
	 * Add agents, specifying their start position, and store their goals on the
	 * opposite side of the environment.
	 */
	sim->addAgent(RVO::Vector2(9.9f, -6.0f));
	goals.push_back(RVO::Vector2(-14.1f, -6.0f));
	// sim->addAgent(RVO::Vector2(9.9f, -5.0f));
	// goals.push_back(RVO::Vector2(-14.1f, -5.0f));
	// sim->addAgent(RVO::Vector2(9.9f, -4.0f));
	// goals.push_back(RVO::Vector2(-14.1f, -4.0f));
	// sim->addAgent(RVO::Vector2(9.9f, 2.0f));
	// goals.push_back(RVO::Vector2(-14.1f, 2.0f));
	// sim->addAgent(RVO::Vector2(-14.1f, 0.0f));
	// goals.push_back(RVO::Vector2(9.9f, 6.0f));
	// sim->addAgent(RVO::Vector2(-14.1f, 3.0f));
	// goals.push_back(RVO::Vector2(9.9f, 8.0f));
	// sim->addAgent(RVO::Vector2(-14.1f, -2.0f));
	// goals.push_back(RVO::Vector2(9.9f, -5.0f));
	sim->addAgent(RVO::Vector2(-14.1f, -2.0f));
	goals.push_back(RVO::Vector2(-14.1f, -2.0f));
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
	// test changing goal
	if (sim->getGlobalTime() == 30.0f) {
		goals[1] = RVO::Vector2(9.9f, -5.0f);
	}
	// ACTUAL DRAFT STARTS HERE
	// agents start waiting - goal the same as original position
	for(int i = 0; i < seqL.size(); i++) {
		if (sim->getGlobalTime() == seqL[i][1]) {
			for (el in posL) {
				for (size_t j = 0; j < sim->getNumAgents(); ++j) {
					if (sim->getAgentPosition(j) != el) {
						position = el
					}
					break // only get the fist fit
				}
			}
			sim->addAgent(RVO::Vector2(position));
			goals.push_back(RVO::Vector2(position));
			waitL.append(seqL[i] + position); // track time and hurry level
		}
	}
	// check each of the waiting agent, shall they run the light or not
	// do the same for left, right, buffer slightly different
	for agent in waitL: // track index by waitL[i]
		if pythonfunction(HurryLevel, currentTime, AppearedTime): // return True if runs the light
			sim->addAgent(RVO::Vector2(posL[i][0]));
			goals.push_back(RVO::Vector2(posL[i][1]));
			waitL.pop(this_agent)
			illegal += 1
	for agent in waitBuffer: // track index by waitL[i]
		if pythonfunction(HurryLevel, currentTime, reachbuffertime): // return True if runs the light
			// change goal, do not know how
			goals.push_back(RVO::Vector2(????);
			waitBuffer.pop(this_agent)
			illegal += 1
	// check all running agent whether the have finished running
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		if not	1.2f < sim->getAgentPosition(i)[x] < 9.9f and not	-14.1f < sim->getAgentPosition(i)[x] < -1.2f:
			illegal -= 1
	}
	// reach buffer area, start another waiting
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		if	-1.2f < sim->getAgentPosition(i)[x] < 1.2f:
			waitBuffer.append(this_agent) // this_agent = [hurrylevel, sim->getGlobalTime()]
	}
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
	return sim->getGlobalTime() > 500.0f;
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
	while (!timeUp(sim));

	delete sim;

	return 0;
}
