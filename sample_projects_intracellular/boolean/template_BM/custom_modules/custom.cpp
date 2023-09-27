#include "./custom.h"
#include <sstream>

/**
 *	\main template-BM custom
 *	\brief Custom module file for template-BM sample_project
 * 
 *	\details Modules needed for the template-BM sample_project. This custom module can be used as template to generate other PhysiBoSS examples.
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	
	// set molecular properties 
	int ainhib_substrate_index = microenvironment.find_density_index( "Ainhib" ); 
	cell_defaults.phenotype.molecular.fraction_released_at_death[ainhib_substrate_index] = 0.0;
	int binhib_substrate_index = microenvironment.find_density_index( "Binhib" ); 
	cell_defaults.phenotype.molecular.fraction_released_at_death[binhib_substrate_index] = 0.0;

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 

	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	

	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}
void update_custom_variables( Cell* pCell )
{
	static int ainhib_index = microenvironment.find_density_index( "Ainhib" ); 
	static int index_ainhib_concentration = pCell->custom_data.find_variable_index("ainhib_concentration");
	static int index_ainhib_node = pCell->custom_data.find_variable_index("ainhib_node");
	pCell->custom_data.variables.at(index_ainhib_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[ainhib_index];
	pCell->custom_data.variables.at(index_ainhib_node).value = pCell->phenotype.intracellular->get_boolean_variable_value("anti_A");

	static int binhib_index = microenvironment.find_density_index( "Binhib" ); 
	static int index_binhib_concentration = pCell->custom_data.find_variable_index("binhib_concentration");
	static int index_binhib_node = pCell->custom_data.find_variable_index("binhib_node");
	pCell->custom_data.variables.at(index_binhib_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[binhib_index];
	pCell->custom_data.variables.at(index_binhib_node).value = pCell->phenotype.intracellular->get_boolean_variable_value("anti_B");
}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	
	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = create_cell(get_cell_definition("default")); 
		pC->assign_position( x, y, z );

		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		update_custom_variables(pC);
	}
	return; 
}

// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	static int index_next_physiboss_run = pCell->custom_data.find_variable_index("next_physiboss_run");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (pCell->phenotype.intracellular->need_update())
	{
		set_input_nodes(pCell);

		pCell->phenotype.intracellular->update();
		
		update_custom_variables(pCell);
		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with live coloring 
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 
	return output; 
}

void set_input_nodes(Cell* pCell) {
	static int ainhib_index = microenvironment.find_density_index( "Ainhib" );
	static double ainhib_threshold = parameters.doubles("ainhib_threshold");
	static int binhib_index = microenvironment.find_density_index( "Binhib" );
	static double binhib_threshold = parameters.doubles("binhib_threshold");

	if (ainhib_index != -1)
	{
		double ainhib_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[ainhib_index];
		pCell->phenotype.intracellular->set_boolean_variable_value("anti_A", ainhib_cell_concentration >= ainhib_threshold);
	}
	
	if (binhib_index != -1)
	{
		double binhib_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[binhib_index];
		pCell->phenotype.intracellular->set_boolean_variable_value("anti_B", binhib_cell_concentration >= binhib_threshold);
	}
}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	double prosurvival_value = pCell->phenotype.intracellular->get_boolean_variable_value("C") ? 1.0 : 0.0;

	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	double multiplier = 1.0;

	// live model 
			
	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );

		multiplier = ( ( prosurvival_value * 20 ) + 1 ); //[1, 21]
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
	}

	pCell->set_internal_uptake_constants(dt);
}

// ***********************************************************
// * NOTE: Funtion replicated from PhysiBoSS, but not used   *
// *       as we use a live cycle model instead a Ki67 model *
// ***********************************************************
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{

}

// ***********************************************************
// * NOTE: Funtion to read init files created with PhysiBoSS *
// ***********************************************************
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header) 
{ 
	// File pointer 
	std::fstream fin; 
	std::vector<init_record> result;

	// Open an existing file 
	fin.open(filename, std::ios::in); 

	// Read the Data from the file 
	// as String Vector 
	std::vector<std::string> row; 
	std::string line, word;

	if(header)
		getline(fin, line);

	do 
	{
		row.clear(); 

		// read an entire row and 
		// store it in a string variable 'line' 
		getline(fin, line);

		// used for breaking words 
		std::stringstream s(line); 

		// read every column data of a row and 
		// store it in a string variable, 'word' 
		while (getline(s, word, delimiter)) { 

			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}
