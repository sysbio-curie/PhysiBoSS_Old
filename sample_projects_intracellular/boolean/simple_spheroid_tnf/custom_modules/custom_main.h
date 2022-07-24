#include "../core/PhysiCell.h"

/**
 *	\main Template-BM custom main file
 *	\brief Custom module file for template-BM sample_project
 * 
 *	\details Modules needed for the template-BM sample_project. This custom module can be used as template to generate other PhysiBoSS examples.
 *
 *	\date 19/10/2020
 *	\author Gerard Pradas and Miguel Ponce de Leon
 */

using namespace BioFVM; 

void inject_density_sphere(int density_index, double concentration, double membrane_lenght);
void remove_density( int density_index );