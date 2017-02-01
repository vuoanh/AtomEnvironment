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
#include "chemistry/bcl_chemistry_atom_environment_two.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_environment_two.cpp
  //!
  //! @author vuot2
  //! @date   08/19/2016
  //! @remarks
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomEnvironmentTwo :
public ExampleInterface
{


public:
  ExampleChemistryAtomEnvironmentTwo *Clone() const
  {
    return new ExampleChemistryAtomEnvironmentTwo( *this);
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
    //const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
    const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "azetin3ol.sdf"));
    BCL_ExampleMustOpenInputFile( input_sdf, filename);

    // load information into small_mol_conformation
    chemistry::FragmentEnsemble molecules( input_sdf);
    // close the input stream
    io::File::CloseClearFStream( input_sdf);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    // iterates through all the fragment complete objects in the fragment ensemble
    for (chemistry::FragmentEnsemble::iterator itr = molecules.Begin(), itr_end = molecules.End(); itr != itr_end; ++itr)
    {
      chemistry::FragmentComplete molecule( *itr);
      chemistry::AtomEnvironment2 element_AE
        ( 2, chemistry::AtomEnvironment2::e_Element, molecule);
      chemistry::AtomEnvironment2 element_AE2( chemistry::AtomEnvironment2::e_Element, element_AE.UnHash());
      std::cout << element_AE.UnHash() << std::endl;
      std::cout << element_AE2.UnHash() << std::endl;
      BCL_ExampleCheck(element_AE,element_AE2);
      //CheckAtomEnv( element_AE);

      chemistry::AtomEnvironment2 elemRC_AE
      ( 2, chemistry::AtomEnvironment2::e_ElemRC, molecule);
      chemistry::AtomEnvironment2 elemRC_AE2( chemistry::AtomEnvironment2::e_ElemRC, elemRC_AE.UnHash() );
      std::cout << elemRC_AE.UnHash() << std::endl;
      std::cout << elemRC_AE2.UnHash() << std::endl;
      BCL_ExampleCheck(elemRC_AE,elemRC_AE2);

      chemistry::AtomEnvironment2 atom_AE
      ( 2, chemistry::AtomEnvironment2::e_Atom, molecule);
      chemistry::AtomEnvironment2 atom_AE2( chemistry::AtomEnvironment2::e_Atom, atom_AE.UnHash() );
      std::cout << atom_AE.UnHash() << std::endl;
      std::cout << atom_AE2.UnHash() << std::endl;
      BCL_ExampleCheck (atom_AE,atom_AE2);

      chemistry::AtomEnvironment2 atomRC_AE
      ( 2, chemistry::AtomEnvironment2::e_AtomRC, molecule);
      chemistry::AtomEnvironment2 atomRC_AE2( chemistry::AtomEnvironment2::e_AtomRC, atomRC_AE.UnHash() );
      std::cout << atomRC_AE.UnHash() << std::endl;
      std::cout << atomRC_AE2.UnHash() << std::endl;
      BCL_ExampleCheck( atomRC_AE, atomRC_AE2);
    }
    return 0;
  } // Run

  //! @brief check if the StringToAE works correctly

  // Print out all the atom environments of a molecule
  static void CheckAtomEnv( const chemistry::AtomEnvironment2 &MOLECULE)
  {
    chemistry::AtomEnvironment2::t_AtomEnvironment atom_env( MOLECULE.GetAtomEnvironment() );



    //check that the number and the identities of the atoms in the atom environment are correct
    std::cout << "The atom of interest is: " << MOLECULE.GetAtomOfInterestIndex() << std::endl;
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
    // check the the string of unhashed atom environment
    std::cout << "the uncoded string elem is: " << MOLECULE.UnHash() << std::endl;
  }

  static const ExampleClass::EnumType s_Instance;

}; //end ExampleChemistryAtomEnvironmentTwo

const ExampleClass::EnumType ExampleChemistryAtomEnvironmentTwo::s_Instance
(
  GetExamples().AddEnum( ExampleChemistryAtomEnvironmentTwo())
);

} // namespace bcl
