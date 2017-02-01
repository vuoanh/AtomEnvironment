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
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)/*

#ifndef BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_
#define BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_
// include the namespace headers
#include "bcl_chemistry.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment_two.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"
// external includes - sorted a;phabetically

namespace bcl
{
  namespace chemistry
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomEnvironment2
    //! @brief stores the atoms whose distances from a atom of interest are no more than a certain number.
    //!
    //! @see @link example_chemistry_atom_environment_oanh.cpp @endlink
    //! @author vuot2
    //! @date 06/29/2016
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MoleculeEnvironment:
      public ConformationComparisonInterface
    {

    public:
      typedef storage::Vector< AtomEnvironment2> t_MoleculeEnv;
      typedef  AtomEnvironment2::AtomTypeEnum t_AtomTypeEnum;
      typedef AtomEnvironment2::Atom_type t_AtomType;
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Elem;
      static const util::SiPtr< const util::ObjectInterface> s_ElemRC;
      static const util::SiPtr< const util::ObjectInterface> s_Atom;
      static const util::SiPtr< const util::ObjectInterface> s_AtomRC;

    private:
      //////////
      // data //
      //////////

      //! @brief atom type as string
      //! @param ATOM_TYPE the name of the atom type
      //! @return the string for the atom type
      static const std::string &GetAtomTypeName( const t_AtomType &ATOM_TYPE);
      t_MoleculeEnv m_MoleculeEnv;
      t_AtomTypeEnum m_AtomType;

    public:

      //////////////////
      // Constructions//
      //////////////////

      //! default constructor
      MoleculeEnvironment();

      //! constructor with Atom_type specification
      MoleculeEnvironment( const t_AtomType &ATOM_TYPE);

      //! @brief constructor for building an Molecule environment from bond distance limit, atom type, and fragment complete
      MoleculeEnvironment ( const t_AtomType &ATOM_TYPE, const FragmentComplete &FRAGMENT );

      MoleculeEnvironment( const MoleculeEnvironment& OTHER);

      //! virtual copy constructor
      MoleculeEnvironment *Clone() const;

      //! @brief compares two atom environments
      bool operator ==(const MoleculeEnvironment &ATOM) const;

      //! @brief compares two atom environments
      bool operator !=(const MoleculeEnvironment &MOLECULE) const;

      //! @brief compute the tanimoto score between this and other molecule environment
      double operator()
      (
        const ConformationInterface& OTHER,
        const ConformationInterface& THIS
      ) const ;


      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instances;

      //////////////////
      // data access ///
      //////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the atom type that this atom enviroment calculator calculates
      const t_AtomType &GetAtomType() const;

      //! @brief returns the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief creates atom environments from every atom of the molecule
      const t_MoleculeEnv &GetMoleculeEnvironment() const;

      // @brief Compute the TanimotoScore between this and OTHER
      double TanimotoScore (const MoleculeEnvironment &OTHER) const;

      protected:

      //////////////////////
      // input and output //
      //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      private:

      // @brief converts vector of sorted objects into a map of objects and their counts
      storage::Map<AtomEnvironment2, size_t> ConvertVectorToMap() const;

      //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
      template< typename K>
      static void AddToMap( const K &KEY, storage::Map<K, size_t> &MAP);
    };
  }
}
#endif /* BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_ */
