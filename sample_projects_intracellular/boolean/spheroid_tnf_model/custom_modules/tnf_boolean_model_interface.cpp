/*
 * tnf_boolean_model_interface.cpp
 *
 *  Created on: 15 jun. 2020
 *  Author: Miguel Ponce-de-Leon (miguel.ponce@bsc.es)
 *  Contributor: Gerard Pradas
 *  Contributor: Arnau Montagud
 *  Contributor: Thalia Diniaco
 *  Cite as: arXiv:2103.14132 [q-bio.QM]
 *  Description: 
 *      Submodel that work as an interface 
 *      between the Boolean Network (BN) and PhysiCell (PC). 
 *      The submodel run the following steps:
 *      1- updates BN input nodes based on custom cell variables (see receptor model)
 *      2- updates the BN intracellular model by running MaBoSS
 *      3- updates cell state/behaviour based on the state of the BN readout nodes
 *  
 *      The update_monitor_variables funtion is just used to keep track of some relevand
 *      BN nodes' values that are stored as custom variables
 */


#include "./tnf_boolean_model_interface.h"

using namespace PhysiCell; 

Submodel_Information tnf_bm_interface_info;

void tnf_boolean_model_interface_setup()
{
    tnf_bm_interface_info.name = "TNF Boolean model interface"; 
	tnf_bm_interface_info.version = "0.2.0";
	
    tnf_bm_interface_info.main_function= tnf_bm_interface_main; 

	// These are just auxiliary variables to keep track of some BN nodes
	tnf_bm_interface_info.cell_variables.push_back( "tnf_node" );
    tnf_bm_interface_info.cell_variables.push_back( "fadd_node" );
    tnf_bm_interface_info.cell_variables.push_back( "nfkb_node" );

	tnf_bm_interface_info.register_model();
}


void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{
    if( pCell->phenotype.death.dead == true )
	{ return; } 

    static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_TNFR" ); 
    static int nTNF_threshold = pCell->custom_data.find_variable_index( "TNFR_activation_threshold" );

    // This if the step transfer function used to update the state of boolean model inputs
    // using the state of the receptor dynamics model. The continuos value thresholded is
    // the total TNF-recptor complex (doi:10.1016/j.cellsig.2010.08.016)
    if ( pCell->custom_data[nR_EB] > pCell->custom_data[nTNF_threshold] )
	{ pCell->phenotype.intracellular->set_boolean_variable_value("TNF", 1); }
	else
    { pCell->phenotype.intracellular->set_boolean_variable_value("TNF", 0); }

    return;

}


void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{	
    static int nTNF_external = microenvironment.find_density_index( "tnf" );
    static int nTNF_export_rate = pCell->custom_data.find_variable_index( "TFN_net_production_rate" );

    static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
    
    // Getting the state of the boolean model readouts (Readout can be in the XML)
    bool apoptosis = pCell->phenotype.intracellular->get_boolean_variable_value( "Apoptosis" );
    bool nonACD = pCell->phenotype.intracellular->get_boolean_variable_value( "NonACD" );
    bool survival =pCell->phenotype.intracellular->get_boolean_variable_value( "Survival" );
    bool NFkB = pCell->phenotype.intracellular->get_boolean_variable_value( "NFkB" );
	
	if ( apoptosis ) 
    {
		pCell->start_death(apoptosis_model_index);
		return;
	}

	if ( nonACD ) 
    {
		pCell->start_death(necrosis_model_index);
		return;
	}

	if ( survival && pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative ) 
    { 
        pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt); 
    }

    // // If NFkB node is active produce some TNF
	if ( NFkB )	
    { 
        phenotype.secretion.net_export_rates[nTNF_external] = pCell->custom_data[nTNF_export_rate]; 
    } else 
    {
        phenotype.secretion.net_export_rates[nTNF_external] = 0;
    }
    
    return;
}


void update_monitor_variables(Cell* pCell )
{
	static int index_tnf_node = pCell->custom_data.find_variable_index("tnf_node");
	static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
	static int index_nfkb_node = pCell->custom_data.find_variable_index("nfkb_node");

	pCell->custom_data[index_nfkb_node] = pCell->phenotype.intracellular->get_boolean_variable_value( "NFkB" ) ;
	pCell->custom_data[index_tnf_node] =pCell->phenotype.intracellular->get_boolean_variable_value("TNF");
	pCell->custom_data[index_fadd_node] = pCell->phenotype.intracellular->get_boolean_variable_value("FADD");

    return;
}


void tnf_bm_interface_main(Cell* pCell, Phenotype& phenotype, double dt)
{
    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

    if ( pCell->phenotype.intracellular->need_update() )
    {
        // First we update the Boolean Model inputs
        update_boolean_model_inputs(pCell, phenotype, dt );
    
        // Run maboss to update the boolean state of the cell
        pCell->phenotype.intracellular->update();

        // update the cell fate based on the boolean outputs
        update_cell_from_boolean_model(pCell, phenotype, dt);

        // Get track of some boolean node values for debugging
        update_monitor_variables(pCell);
    }

    return;
}
