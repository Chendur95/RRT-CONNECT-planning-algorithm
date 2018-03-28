#include "planner.h"
#include <time.h>

static void plannerRRTCONNECT(
	double*  map,
	int x_size,
	int y_size,
	double* armstart_anglesV_rad,
	double* armgoal_anglesV_rad,
	int numofDOFs,
	double*** plan,
	int* planlength)
{
	*plan = NULL;
	*planlength = 0;
	int n_samples = 0;
	srand(time(NULL)); 

	RRTTree startTree(numofDOFs, armstart_anglesV_rad);
	RRTTree goalTree(numofDOFs, armgoal_anglesV_rad);

	bool pathFound = false; //path not found by default
	double* rand_angles = (double*)malloc(numofDOFs * sizeof(double));
	bool extend_stree = true; //To change trees

	while (!pathFound && n_samples < 50000)
	{
		n_samples++;
		if (extend_stree)
		{
			//Code to extend startTree to random_angles
			GenRandomConfig(rand_angles, map, numofDOFs, x_size, y_size);
			int start_nbid = startTree.GetNearestNeighbour(rand_angles);
			double* start_extend_angles_rad = Extend(startTree.get_vertice(start_nbid), rand_angles, map, numofDOFs, x_size, y_size, Ephsilon);
			if (start_extend_angles_rad == NULL) 
				continue;
			int start_extendid = startTree.add_vertex(start_extend_angles_rad);
			startTree.add_edge(start_nbid, start_extendid);

			//Code to extend goalTree to extend_angless
			int goal_nb_id = goalTree.GetNearestNeighbour(start_extend_angles_rad);
			double* goal_extend_angles_rad = Extend(goalTree.get_vertice(goal_nb_id), start_extend_angles_rad, map, numofDOFs, x_size, y_size, Ephsilon);
			if (goal_extend_angles_rad == NULL) 
				continue;
			int goal_extend_id = goalTree.add_vertex(goal_extend_angles_rad);
			goalTree.add_edge(goal_nb_id, goal_extend_id);
			if (angles_equal(start_extend_angles_rad, goal_extend_angles_rad, numofDOFs))
			{
				mexPrintf("Reached\n");
				pathFound = true;
				break;
			}
		}

		else
		{
			//Extend goal tree to random_angles
			GenRandomConfig(rand_angles, map, numofDOFs, x_size, y_size);
			int goal_nb_id = goalTree.GetNearestNeighbour(rand_angles);
			double* goal_extend_angles_rad = Extend(goalTree.get_vertice(goal_nb_id), rand_angles, map, numofDOFs, x_size, y_size, Ephsilon);
			if (goal_extend_angles_rad == NULL) 
				continue;
			int goal_extend_id = goalTree.add_vertex(goal_extend_angles_rad);
			goalTree.add_edge(goal_nb_id, goal_extend_id);
			
			//Extend start tree to extend_angles
			int start_nb_id = startTree.GetNearestNeighbour(goal_extend_angles_rad);
			double* start_extend_angles_rad = Extend(startTree.get_vertice(start_nb_id), goal_extend_angles_rad, map, numofDOFs, x_size, y_size, Ephsilon);
			if (start_extend_angles_rad == NULL) 
				continue;
			int start_extend_id = startTree.add_vertex(start_extend_angles_rad);
			startTree.add_edge(start_nb_id, start_extend_id);
			if (angles_equal(goal_extend_angles_rad, start_extend_angles_rad, numofDOFs))
			{
				mexPrintf("Reached\n");
				pathFound = true;
				break;
			}
		}
		extend_stree = !extend_stree;
	}

	if (pathFound)
	{
		vector<int> path;
		int next_Index = startTree.get_last_index();
		while (next_Index != 0)
		{
			path.insert(path.begin(), next_Index);
			next_Index = startTree.get_prnt(next_Index);
		}
		path.insert(path.begin(), startTree.get_root_index());

		int startTree_size = path.size();

		next_Index = goalTree.get_last_index();

		while (next_Index != 0)
		{
			next_Index = goalTree.get_prnt(next_Index);
			path.push_back(next_Index);
		}

		*planlength = path.size();

		*plan = (double**)malloc(path.size() * sizeof(double*));

		for (int i = 0; i<path.size(); i++)
		{
			(*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
			if (i<startTree_size)
				memcpy((*plan)[i], startTree.get_vertice(path[i]), numofDOFs * sizeof(double));
			else
				memcpy((*plan)[i], goalTree.get_vertice(path[i]), numofDOFs * sizeof(double));
		}

	}

	else
	{
		mexPrintf("Path not found!\n");
	}

	mexPrintf("Num of Samples: %d\n", n_samples);
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	clock_t t, tend;
	/* Check for proper number of arguments */
	if (nrhs != 3) {
		mexErrMsgIdAndTxt("MATLAB:planner:invalidNumInputs",
			"Three input arguments required.");
	}
	else if (nlhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:planner:maxlhs",
			"One output argument required.");
	}

	/* get the dimensions of the map and the map matrix itself*/
	int x_size = (int)mxGetM(MAP_IN);
	int y_size = (int)mxGetN(MAP_IN);
	double* map = mxGetPr(MAP_IN);

	/* get the start and goal angles*/
	int numofDOFs = (int)(MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
	if (numofDOFs <= 1) {
		mexErrMsgIdAndTxt("MATLAB:planner:invalidnumofdofs",
			"it should be at least 2");
	}
	double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
	if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))) {
		mexErrMsgIdAndTxt("MATLAB:planner:invalidnumofdofs",
			"numofDOFs in startangles is different from goalangles");
	}
	double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);

	double** plan = NULL;
	int planlength = 0;

	bool valid_startgoal = true;

	//Check the start and goal configuration.
	if (!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size) ||
		!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size))
	{
		mexPrintf("Not a valid start or goal config. \n");
		valid_startgoal = false;
	}

	if (valid_startgoal)
	{
		double t = clock();
		plannerRRTCONNECT(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
		double tend = clock(); 
		double time_spent = (double)(tend - t);
		mexPrintf("The planning takes %f \n", time_spent );
		printf("Plan is of length=%d \n", planlength);
	}

	/* Create return values */
	if (planlength > 0)
	{
		double plan_distance = 0;

		for (int i = 0; i<planlength - 1; i++)
		{
			plan_distance = plan_distance + dist_angles(plan[i], plan[i + 1], numofDOFs);
		}

		mexPrintf("Total distance of plan = %f\n", plan_distance);

		PLAN_OUT = mxCreateNumericMatrix((mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL);
		double* plan_out = mxGetPr(PLAN_OUT);
		//copy the values
		int i, j;
		for (i = 0; i < planlength; i++)
		{
			for (j = 0; j < numofDOFs; j++)
			{
				plan_out[j*planlength + i] = plan[i][j];
			}
		}
	}
	else
	{
		PLAN_OUT = mxCreateNumericMatrix((mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL);
		double* plan_out = mxGetPr(PLAN_OUT);
		//copy the values
		int j;
		for (j = 0; j < numofDOFs; j++)
		{
			plan_out[j] = armstart_anglesV_rad[j];
		}
	}
	PLANLENGTH_OUT = mxCreateNumericMatrix((mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL);
	int* planlength_out = (int*)mxGetPr(PLANLENGTH_OUT);
	*planlength_out = planlength;


	return;

}
