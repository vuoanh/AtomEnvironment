// (c) Copyright BCL @ Vanderbilt University 2015
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

// include the namespace headers
#include "chemistry/bcl_chemistry.h"
// include the header of this class
#include "chemistry/bcl_chemistry_molecule_environment.h"
#include "util/bcl_util_serializable_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment_two.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_fragment_complete.h" //!> input is a fragment complete
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_enumerated.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_function_wrapper.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //////////
    // data //
    //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_Elem
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironment2::e_Element))
    );
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_ElemRC
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironment2::e_ElemRC))
    );
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_Atom
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironment2::e_Atom))
    );
    const util::SiPtr< const util::ObjectInterface> MoleculeEnvironment::s_AtomRC
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MoleculeEnvironment( AtomEnvironment2::e_AtomRC))
    );

    //! @brief atom type name as string
    //! @param ATOM_TYPE the atom type
    //! @return the string for the measure
    const std::string &MoleculeEnvironment::GetAtomTypeName( const t_AtomType &ATOM_TYPE)
    {
      return AtomEnvironment2::GetAtomTypeInfo( ATOM_TYPE).First();
    }

    //////////////////
    // Constructions//
    //////////////////

    //! constructor with Atom_type specification
    MoleculeEnvironment::MoleculeEnvironment( const t_AtomType &ATOM_TYPE): m_AtomType( ATOM_TYPE)
    {
    }

    //! @brief constructor for building an atom environment from bond distance limit, index of atom, atom type, and fragment complete
    MoleculeEnvironment::MoleculeEnvironment( const t_AtomType &ATOM_TYPE, const FragmentComplete & FRAGMENT):
                                             m_AtomType( ATOM_TYPE)
    {
      size_t atom_num = FRAGMENT.GetSize();
      for( size_t i = 0; i < atom_num; ++i)
      {
        if( FRAGMENT.GetAtomVector()( i).GetAtomType()->GetElementType() != GetElementTypes().e_Hydrogen)
        {
          AtomEnvironment2 AtomEnv( i, m_AtomType.GetEnum(), FRAGMENT);
          m_MoleculeEnv.PushBack( AtomEnv);
        }
      }
      m_MoleculeEnv.Sort( std::less< AtomEnvironment2>() );
    }

    MoleculeEnvironment::MoleculeEnvironment( const MoleculeEnvironment& OTHER):
        m_MoleculeEnv(OTHER.m_MoleculeEnv), m_AtomType( OTHER.m_AtomType)
    {
    }

    //! copy constructor
    MoleculeEnvironment *MoleculeEnvironment::Clone() const
    {
      return new MoleculeEnvironment( *this);
    }

    //! @brief compares two atom environments
    bool MoleculeEnvironment::operator ==( const MoleculeEnvironment &MOLECULE) const
    {
      return ( m_MoleculeEnv == MOLECULE.m_MoleculeEnv);
    }

    //! @brief compares two atom environments
    bool MoleculeEnvironment::operator !=(const MoleculeEnvironment &MOLECULE) const
    {
      return (m_MoleculeEnv != MOLECULE.m_MoleculeEnv);
    }

    //! @brief compute the tanimoto score between molecule environments of THIS and OTHER
    double MoleculeEnvironment::operator()
    (
      const ConformationInterface& THIS,
      const ConformationInterface& OTHER
    ) const
    {
      MoleculeEnvironment this_molecule( m_AtomType, FragmentComplete(THIS));
      MoleculeEnvironment other_molecule( m_AtomType, FragmentComplete(OTHER));
      return this_molecule.TanimotoScore( other_molecule);
    }
    //! @brief compute the tanimoto score between this and other molecule environment
    //double operator()

    //////////////////////
    // input and output //
    //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeEnvironment::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( AtomEnvironment2::GetAtomTypeInfo( m_AtomType).Second());
      return serializer;
    }
    //////////////////
    // data access ///
    //////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeEnvironment::GetClassIdentifier() const
    {
      return GetStaticClassName(*this);
    }

    //! @brief get the atom type that this atom enviroment calculator calculates
    const MoleculeEnvironment::t_AtomType &MoleculeEnvironment::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns the data label
    //! @return data label as string
    const std::string &MoleculeEnvironment::GetAlias() const
    {
      return ( AtomEnvironment2::GetAtomTypeInfo( m_AtomType).First() );
    }

    //! @brief returns the atom environment
    const MoleculeEnvironment::t_MoleculeEnv &MoleculeEnvironment::GetMoleculeEnvironment() const
    {
      return m_MoleculeEnv;
    }

    // @brief converts vector of sorted objects into a map of objects and their counts
    storage::Map<AtomEnvironment2, size_t> MoleculeEnvironment::ConvertVectorToMap() const
    {
      storage::Map<AtomEnvironment2, size_t> map;
      for (t_MoleculeEnv::const_iterator iter = m_MoleculeEnv.Begin(), end_iter = m_MoleculeEnv.End(); iter != end_iter; ++iter)
      {
        AddToMap< AtomEnvironment2>( *iter, map);
      }
      return map;
    }

    //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
    template< typename K>
    void MoleculeEnvironment::AddToMap( const K &KEY, storage::Map< K, size_t> &MAP)
    {
      std::pair< typename storage::Map< K, size_t>::iterator, bool> insert_pair
      (
        MAP.Insert( storage::Pair< K, size_t> ( KEY, size_t(1) ) )
      );
      if( ! ( insert_pair.second ) )
      {
        insert_pair.first->second += 1;
      }
    }

    // @brief Compute the TanimotoScore between this and OTHER
    double MoleculeEnvironment::TanimotoScore (const MoleculeEnvironment &OTHER) const
    {
      double sumA ( GetMoleculeEnvironment().GetSize());
      double sumB ( OTHER.GetMoleculeEnvironment().GetSize());
      storage::Map< AtomEnvironment2, size_t> lhs( ConvertVectorToMap() );
      t_MoleculeEnv key_lhs( lhs.GetKeysAsVector() );
      storage::Map< AtomEnvironment2, size_t> rhs( OTHER.ConvertVectorToMap() );
      t_MoleculeEnv key_rhs( rhs.GetKeysAsVector() );
      std::vector< AtomEnvironment2> shared_atom;
      std::set_intersection
      (
        key_lhs.Begin(), key_lhs.End(), key_rhs.Begin(), key_rhs.End(), std::back_inserter( shared_atom)
      );
      double sumAB( 0);
      for ( std::vector< AtomEnvironment2>::const_iterator iter = shared_atom.begin(), end_iter = shared_atom.end(); iter != end_iter; ++iter)
      {
        size_t A( lhs.Find( *iter)->second);
        size_t B( rhs.Find( *iter)->second);
        sumAB += std::min(A, B);
      }
      //std::cout << sumA << " " << sumB << " " << sumAB << std::endl;
      if (! sumAB) return 0;                                          //!> If sumAB != 0
      return double( sumAB / (sumA + sumB - sumAB));
    }
  }
}





