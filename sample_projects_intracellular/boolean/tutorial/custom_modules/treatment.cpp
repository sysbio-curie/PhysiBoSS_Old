#include "treatment.h"

void treatment_function () 
{
	if (PhysiCell::parameters.bools.find_index("treatment") != -1) 
	{
		int treatment_substrate_index = BioFVM::microenvironment.find_density_index(PhysiCell::parameters.strings("treatment_substrate"));

		if (PhysiCell::parameters.bools("treatment")){
		
			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == 0 
				&& !BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("treatment_substrate") << " activation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, true);	
			}

			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == PhysiCell::parameters.ints("treatment_duration") 
				&& BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
			}
			
		} else if ( BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index) ){
			std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation (NO TREATMENT) at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
			BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
		}
	}
}