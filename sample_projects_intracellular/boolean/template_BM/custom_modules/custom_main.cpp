#include "./custom_main.h"

/**
 *	\main Template-BM custom main file
 *	\brief Custom module file for template-BM sample_project
 * 
 *	\details Modules needed for the template-BM sample_project. This custom module can be used as template to generate other PhysiBoSS examples.
 *
 *	\date 19/10/2020
 *	\author Gerard Pradas and Miguel Ponce de Leon
 */

void inject_density_sphere(int density_index, double concentration, double membrane_lenght) 
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration; 	
	}
}

void remove_density( int density_index )
{	
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
		microenvironment.density_vector(n)[density_index] = 0; 	
	std::cout << "Removal done" << std::endl;
}