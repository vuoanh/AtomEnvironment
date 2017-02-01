

// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users:
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users:
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University:
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University.
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

// include example header
#include "example.h"

// include the header of the class which this example is for
#include "chemistry/bcl_chemistry_molecule_environment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment_two.h"
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecule_environment.cpp
  //!
  //! @author vuot2
  //! @date   09/01/2016
  //! @remarks
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMoleculeEnvironment :
public ExampleInterface
{
public:

  ExampleChemistryMoleculeEnvironment *Clone() const
  {
    return new ExampleChemistryMoleculeEnvironment( *this);
  }

  /////////////////
  // data access //
  /////////////////

  //! @brief returns class name
  //! @return the class name as const ref std::string
  const std::string &GetClassIdentifier() const
  {
    return GetStaticClassName( *this);
  }

  int Run() const
  {
    // read sdf file
    io::IFStream input_sdf;
    //const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "1798_inactives_clean.sdf.gz"));
    const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "463087_actives_clean.sdf" ) );
    //const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "azetin3ol.sdf"));
    //const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
    //const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
    BCL_ExampleMustOpenInputFile( input_sdf, filename);

    // load information into small_mol_conformation
    chemistry::FragmentEnsemble molecules( input_sdf);
    // close the input stream
    io::File::CloseClearFStream( input_sdf);

    // read sdf file

    io::IFStream input_sdf1;
    const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
    BCL_ExampleMustOpenInputFile( input_sdf1, hexane_filename);

    // load information into small_mol_conformation
    chemistry::FragmentComplete cyclohexane (sdf::FragmentFactory::MakeFragment( input_sdf1) );
    // close the input stream
    io::File::CloseClearFStream( input_sdf1);
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    //chemistry::MoleculeEnvironment cyclohexane_elem( chemistry::AtomEnvironment2::e_Element, cyclohexane);
    //chemistry::MoleculeEnvironment cyclohexane_elem( chemistry::AtomEnvironment2::e_ElemRC);
    // CheckMoleculeEnv( cyclohexane_elem);

    // iterates through all the fragment complete objects in the fragment ensemble
    for (chemistry::FragmentEnsemble::const_iterator itr = molecules.Begin(), itr_end = molecules.End(); itr != itr_end; ++itr)
    {
      chemistry::FragmentComplete molecule( *itr);
      chemistry::MoleculeEnvironment environment_elem( chemistry::AtomEnvironment2::e_ElemRC, molecule);
      CheckMoleculeEnv( environment_elem);
      // Checks if the Tanimoto index of two identical molecular atom environments is 1
      //std::cout<< environment_elem.TanimotoScore( environment_elem) << std::endl;
      // Checks if the Tanimoto index of two identical molecular atom environments is 1
      //std::cout<< cyclohexane_elem(molecule, cyclohexane) << std::endl;
      //std::cout<< environment_elem.TanimotoScore( cyclohexane_elem) << std::endl;
    }
    return 0;
  } // Run

  // Print out all the atom environments of a molecule
  static void CheckMoleculeEnv( const chemistry::MoleculeEnvironment &MOL_ENV)
  {
    chemistry::MoleculeEnvironment::t_MoleculeEnv mol_env(MOL_ENV.GetMoleculeEnvironment());
    for (chemistry::MoleculeEnvironment::t_MoleculeEnv::iterator iter = mol_env.Begin(), end_iter = mol_env.End(); iter != end_iter; ++iter)
    {
      /*chemistry::AtomEnvironment2::t_AtomEnvironment atom_env( iter->GetAtomEnvironment() );

      //check that the number and the identities of the atoms in the atom environment are correct
      std::cout << "The atom of interest is: " << iter->GetAtomOfInterestIndex() << std::endl;
      int i = 0;
      for (chemistry::AtomEnvironment2::t_AtomEnvironment::const_iterator itr = atom_env.Begin(), itr_end = atom_env.End(); itr != itr_end; ++itr, ++i)
      {
        std::cout << i << ": ";
        for (storage::Map< size_t, size_t>::const_iterator map_itr = itr->Begin(), map_itr_end = itr->End(); map_itr!= map_itr_end; ++map_itr)
        {
          std::cout << map_itr->first << "(" << map_itr->second << ") ";
        }
        std::cout <<std::endl;
      }
      // check the the string of unhashed atom environment */
      std::cout << iter->UnHash() << std::endl;
    }
  }

  static const ExampleClass::EnumType s_Instance;

}; //end ExampleChemistryMoleculeEnvironment

const ExampleClass::EnumType ExampleChemistryMoleculeEnvironment::s_Instance
(
  GetExamples().AddEnum( ExampleChemistryMoleculeEnvironment())
);

} // namespace bcl
