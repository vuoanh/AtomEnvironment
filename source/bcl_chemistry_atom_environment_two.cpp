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
#include "chemistry/bcl_chemistry_atom_environment_two.h"
#include "util/bcl_util_serializable_interface.h"

// includes from bcl - sorted alphabetically
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
#include "function/bcl_function_member_unary_const.h"
#include "util/bcl_util_string_functions.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "sdf/bcl_sdf_bond_info.h"

using bcl::sdf::BondInfo;
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //////////
    // data //
    //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t atom_type( 0); atom_type < AtomEnvironment2::s_NumberAtomTypes; ++atom_type)
        {
          last_instance =
              util::Enumerated< util::SerializableInterface>::AddInstance
              (
                new AtomEnvironment2( static_cast< AtomEnvironment2::Atom_type>( atom_type))
              );
        }
        return last_instance;
      }
    }

    //single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomEnvironment2::s_Instances( AddInstances());

    //! @brief Get the name, description, and function of the given atom type
    //! @param ATOM_TYPE the atom type of interest
    //! @return the short name or abbreviation for the class
    const storage::Triplet< std::string, std::string, AtomEnvironment2::t_functions >
    & AtomEnvironment2::GetAtomTypeInfo( const Atom_type &ATOM_TYPE)
    {
      typedef storage::Triplet< std::string, std::string, AtomEnvironment2::t_functions > t_Info;
      AtomEnvironment2::t_functions  elem_function( &chemistry::AtomEnvironment2::ElementHash, &chemistry::AtomEnvironment2::UnElemHash );
      AtomEnvironment2::t_functions  elem_functionRC( &chemistry::AtomEnvironment2::ElemRCHash, &chemistry::AtomEnvironment2::UnElemRCHash );
      AtomEnvironment2::t_functions  atom_function( &chemistry::AtomEnvironment2::AtomHash, &chemistry::AtomEnvironment2::UnAtomHash );
      AtomEnvironment2::t_functions  atom_functionRC( &chemistry::AtomEnvironment2::AtomRCHash, &chemistry::AtomEnvironment2::UnAtomRCHash );
      static const t_Info s_info[ s_NumberAtomTypes] =
      {
          t_Info( "Element", "Atomic Number and bond orders", elem_function),
          t_Info( "ElemRC", "Atomic number, bond orders, and whether the atom is in a ring or/and aromatic", elem_functionRC),
          t_Info( "Atom", "Atomic Number and bond orders", atom_function),
          t_Info( "AtomRC", "Atomic Number and bond orders", atom_functionRC),
      };
      return s_info[ ATOM_TYPE];
    };

    //! @brief atom type name as string
    //! @param ATOM_TYPE the atom type
    //! @return the string for the measure
    const std::string &AtomEnvironment2::GetAtomTypeName( const Atom_type &ATOM_TYPE)
    {
      return GetAtomTypeInfo( ATOM_TYPE).First();
    }

    //////////////////
    // Constructions//
    //////////////////

    //! default constructor
    AtomEnvironment2::AtomEnvironment2():m_BondNumLimit(2)
    {
    }

    //! constructor with Atom_type specification
    AtomEnvironment2::AtomEnvironment2( const Atom_type &ATOM_TYPE): m_BondNumLimit(2), m_AtomType( ATOM_TYPE)
    {
    }

    //! constructor from AE string representation
    AtomEnvironment2::AtomEnvironment2( const Atom_type &ATOM_TYPE, const std::string &STRING):
        m_BondNumLimit(2),
        m_AtomType( ATOM_TYPE),
        m_AtomEnvironment(StringToAE(STRING, ATOM_TYPE))
    {
    }

    //! @brief constructor for building an atom environment from bond distance limit, index of atom, atom type, and fragment complete
    AtomEnvironment2::AtomEnvironment2( size_t ATOM_INDEX, const Atom_type &ATOM_TYPE, const ConformationInterface & FRAGMENT):
                                m_AtomIndex( ATOM_INDEX), m_BondNumLimit(2), m_AtomType( ATOM_TYPE),
                                m_AtomEnvironment(m_BondNumLimit + 1, storage::Map< size_t, size_t>())
    {
      t_output s_output(AtomEnvironment2::GetEncodedAtomEnvironment( FRAGMENT));
      storage::Vector< std::string > encoded_atom_environment = s_output.First();
      storage::Vector< size_t> bond_orders = s_output.Second();
      DecodeAtomEnvironment
      (
        encoded_atom_environment,
        HashMolecule( FRAGMENT, bond_orders),
        m_BondNumLimit
      );
    }

    //! copy constructor
    AtomEnvironment2 *AtomEnvironment2::Clone() const
    {
      return new AtomEnvironment2( *this);
    }

    //! @brief compares two atom environments
    bool AtomEnvironment2::operator ==( const AtomEnvironment2 &ATOM) const
    {
      return ( m_AtomEnvironment == ATOM.m_AtomEnvironment);
    }

    //! @brief compares two atom environments
    bool AtomEnvironment2::operator <( const AtomEnvironment2 &ATOM) const
    {
      return ( m_AtomEnvironment < ATOM.m_AtomEnvironment);
    }

    //! @brief compares two atom environments
    bool AtomEnvironment2::operator !=(const AtomEnvironment2 &ATOM) const
    {
      return (m_AtomEnvironment != ATOM.m_AtomEnvironment);
    }

    //////////////////////
    // input and output //
    //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomEnvironment2::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetAtomTypeInfo( m_AtomType).Second());
      return serializer;
    }
    //////////////////
    // data access ///
    //////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomEnvironment2::GetClassIdentifier() const
    {
      return GetStaticClassName(*this);
    }

    //! @brief access the atom of interest
    std::size_t AtomEnvironment2::GetAtomOfInterestIndex() const
    {
      return m_AtomIndex;
    }

    //! @brief get the atom type that this atom enviroment calculator calculates
    const AtomEnvironment2::Atom_type &AtomEnvironment2::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief returns the data label
    //! @return data label as string
    const std::string &AtomEnvironment2::GetAlias() const
    {
      return ( GetAtomTypeInfo( m_AtomType).First() );
    }

    //! @brief returns the atom environment
    const AtomEnvironment2::t_AtomEnvironment &AtomEnvironment2::GetAtomEnvironment() const
    {
      return m_AtomEnvironment;
    }

    //! @brief read the fragment complete of the molecule and and the index of the atom of interest,
    //!        and then provide encoded version of the atom environment.
    //! @param FRAGMENT the fragment complete representing the molecule
    //! @param BOND_LIMIT the limit number of bond away from the atom of interest
    AtomEnvironment2::t_output AtomEnvironment2::GetEncodedAtomEnvironment
    (
      const ConformationInterface &FRAGMENT
    )
    {
      size_t moleculeSize = FRAGMENT.GetSize(); //!> get the number of atoms in the molecule
      //!> initializes encoded environment with all '0' values
      storage::Vector< std::string> encoded_atom_environment(m_BondNumLimit + 1, std::string(moleculeSize, '0'));
      storage::Vector< size_t> bond_orders( moleculeSize, size_t(0));
      encoded_atom_environment(0)[ m_AtomIndex ] = '1';  // the atom itself would be the only '1' in the first string
      // this function will recursively do the rest of work
      GetEncodedAtomEnvironmentHelper(FRAGMENT, m_AtomIndex, m_BondNumLimit, 1, encoded_atom_environment, bond_orders);
      AtomEnvironment2::t_output s_output( encoded_atom_environment, bond_orders);
      return s_output;
    }

    //! @brief recursively find the neighbor atoms in the atom environment of the atom of interest and then update
    //! @param ENCODED_ATOM_ENVIRONMENT the encoded atom enviroment, each string represents a phere,
    //          each character represent absence/presece of the corresponding atom
    void AtomEnvironment2::GetEncodedAtomEnvironmentHelper
    (
      const ConformationInterface &FRAGMENT, size_t ATOM_INDEX, size_t BOND_LIMIT,
      size_t BOND_NUM, storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      storage::Vector< size_t> &BOND_ORDERS
    )
    {
      if (BOND_NUM <= BOND_LIMIT)
      {
        const storage::Vector< BondConformational> &bond_vector (FRAGMENT.GetAtomsIterator() + ATOM_INDEX).GetBonds());
        for
        (
            storage::Vector<BondConformational>::const_iterator itr = bond_vector.Begin(),
            itr_end = FRAGMENT.GetBondInfo().End(); itr != itr_end;
            ++itr
        )
        {
          if ( itr->IsValid() && ( itr->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen) )
          {
            size_t bond_order(itr->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrder) );
            size_t target_atom_index = FRAGMENT.GetAtomIndex( itr->GetTargetAtom());
            UpdateEncodedAtomEnvironment( BOND_NUM, target_atom_index, ENCODED_ATOM_ENVIRONMENT, bond_order, BOND_ORDERS);
            GetEncodedAtomEnvironmentHelper
            (
              FRAGMENT, target_atom_index, BOND_LIMIT, BOND_NUM + 1, ENCODED_ATOM_ENVIRONMENT, BOND_ORDERS
            );
          }
        }
      }
    }
    //! @brief returns the corresponding character representing the bond order
    //! @param BOND_ORDER a integer representing the bond order
    std::string AtomEnvironment2::GetBondOrder(size_t BOND_ORDER)
    {
      switch( BOND_ORDER)
      {
        case 1: return "-";
        case 2: return "=";
        case 3: return "#";
        default: return "_";
      }
    }

    //! @brief updates the encoded atom environment
    //! @param TARGET_ATOM_INDEX the index of the target atom, which bond to the atom of interest
    //! @param ENCODED_ATOM_ENVIRONMENT a vector of string that represents the neighbor atoms in different spheres
    void AtomEnvironment2::UpdateEncodedAtomEnvironment
    (
      size_t BOND_NUM, size_t TARGET_ATOM_INDEX,
      storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      size_t BOND_ORDER, storage::Vector< size_t> &BOND_ORDERS
    )
    {
      if ( ( ENCODED_ATOM_ENVIRONMENT( BOND_NUM)[ TARGET_ATOM_INDEX] != '0') && ( BOND_ORDERS( TARGET_ATOM_INDEX) < BOND_ORDER))
      {
        BOND_ORDERS( TARGET_ATOM_INDEX) = BOND_ORDER;
        return;
      }
      bool update = true;
      for ( size_t i = 0; i <= BOND_NUM; i++)
      {
        if ( ENCODED_ATOM_ENVIRONMENT( i)[ TARGET_ATOM_INDEX] != '0')
        {
          update = false;
          break;
        }
      }
      if ( update)
      {
        ENCODED_ATOM_ENVIRONMENT( BOND_NUM)[ TARGET_ATOM_INDEX] = '1';
        BOND_ORDERS( TARGET_ATOM_INDEX) = BOND_ORDER;
      }
    }

    //! @brief converts from the encoded atom environment to the normal representation of the atom environment
    //! @param ENCODED_ATOM_ENVIRONMENT encoded atom environment in term of a vector of string of 0 and 1
    //! @param HASHED_MOLECULE is vector of hashed atoms of the molecule
    void AtomEnvironment2::DecodeAtomEnvironment
    (
      const storage::Vector< std::string >&ENCODED_ATOM_ENVIRONMENT,
      const storage::Vector< size_t >&HASHED_MOLECULE,
      size_t BOND_LIMIT
    )
    {
      size_t size( HASHED_MOLECULE.GetSize());
      for( size_t bond_num( 0); bond_num <= BOND_LIMIT; bond_num++)
      {
        for( size_t atom_index = 0; atom_index <= size;  atom_index++)
        {
          if(ENCODED_ATOM_ENVIRONMENT( bond_num)[ atom_index] == '1')
          {
            this->AddAtom( HASHED_MOLECULE( atom_index), bond_num);   //!> adds the hashed atom into m_AtomEnvironment
          }
        }
      }
    }

    //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
    template< typename K>
    void AtomEnvironment2::AddToMap( const K &KEY, storage::Map<K, size_t> &MAP)
    {
      std::pair< typename storage::Map< K, size_t>::iterator, bool> insert_pair
      (
        MAP.Insert( storage::Pair<K, size_t> (KEY, size_t(1) ) )
      );
      if( ! ( insert_pair.second ) )
      {
        insert_pair.first->second += 1;
      }
    }

    //! @brief adds the hashed atom into m_AtomEnvironment
    //! @param HASHED_ATOM
    //! @param BOND_DISTANCE
    void AtomEnvironment2::AddAtom( size_t HASHED_ATOM, size_t BOND_DISTANCE)
    {
      AddToMap< size_t>( HASHED_ATOM, m_AtomEnvironment( BOND_DISTANCE));
    }

    ////////////////////////////////////////
    // helper functions for hash functions /
    ////////////////////////////////////////

    //! @brief returns the chemical symbol of an element from its atomic number ATOM_NUMBER
    const std::string &AtomEnvironment2::GetAtomicSymbol( size_t ATOM_NUMBER, bool ATOM_FLAG)
    {
      if (ATOM_FLAG) // if this fingerprint uses atom type
      {
        AtomType atom( ATOM_NUMBER);
        return atom.GetName();
      } else // if this fingerprint uses element type
      {
        ElementType element( ATOM_NUMBER);
        return element->GetChemicalSymbol();
      }
    }

    //! @brief return the enum index of the atom type
    const size_t AtomEnvironment2::GetAtomTypeIndex( const AtomConformationalInterface &ATOM)
    {
      return ATOM.GetAtomType()->GetIndex();
    }

    //! @brief return the enum index of the element type
    //! @param ATOM: the index of the atom of interest
    const size_t AtomEnvironment2::GetAtomicNumber( const AtomConformationalInterface &ATOM)
    {
      return ATOM.GetElementType().GetIndex();
    }

    //! @brief Checks if the atom is in a ring or not
    //! @return 1 if the atom in a ring, 0 otherwise
    const size_t AtomEnvironment2::RingOrAromatic( const AtomConformationalInterface &ATOM)
    {
      if( ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1) > 0)
      {
        return 2;
      }
      if( ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1) > 0)
      {
        return 1;
      }
      return 0;
    }

    //! @brief return the string symbol from the "ring or aromatic" information
    const std::string AtomEnvironment2::GetRCSymbol( size_t RC)
    {
      if (RC == 2)
      {
        return ".r";
      }
      if (RC == 1)
      {
        return ".c";
      }
      return "";
    }

    ///////////////////////////////
    // Hash and unhash functions /
    //////////////////////////////

    //! @brief choose the hash function based on the choice of atom type, then hash every atoms of a fragment complete
    //! @return a number which is a compressed version of the atom type
    storage::Vector<size_t> AtomEnvironment2::HashMolecule
    (
      const ConformationInterface &FRAGMENT,
      const storage::Vector< size_t> &BOND_ORDERS
    ) const
    {
      size_t size = FRAGMENT.GetNumberAtoms();
      storage::Vector< size_t> hashedAtoms( size, size_t(0));
      size_t index(0);
      for
      (
          iterate::Generic< const AtomConformationalInterface> itr( FRAGMENT.GetAtomsIterator());
          itr.NotAtEnd();
          ++itr, ++index
      )
      {
        function::MemberUnaryConst<AtomEnvironment2, const  AtomEnvironment2::t_input&, size_t>
        function( GetAtomTypeInfo( m_AtomType).Third().First());
        t_input s_input( *itr, BOND_ORDERS(index) );
        hashedAtoms( index) = function( *this, s_input); //!> add the hashed atom into the corresponding array element
      }
      return hashedAtoms;
    }

    //! @brief unhashes the atom environment into its string representation
    std::string AtomEnvironment2::UnHash() const
    {
      function::MemberUnaryConst<AtomEnvironment2, const size_t&, std::string> function( GetAtomTypeInfo( m_AtomType).Third().Second());
      std::ostringstream oss;
      for ( size_t i = 0; i <= m_BondNumLimit; ++i)
      {
        oss << "[";
        for ( storage::Map< size_t, size_t>::const_iterator itr = m_AtomEnvironment( i).Begin(), itr_end = m_AtomEnvironment( i).End(); itr != itr_end; ++itr)
        {
          for ( size_t j = 0; j < itr->second; ++j)
          {
            oss << function( *this, itr->first) + " ";
          }
        }
        oss << "]";
      }
      return oss.str();
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironment2::ElementHash( const t_input &INPUT) const
    {
      return ( INPUT.Second() << 7) | GetAtomicNumber( INPUT.First() );
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironment2::ElemRCHash( const t_input &INPUT) const
    {
      return ( AtomEnvironment2::RingOrAromatic( INPUT.First()) << 9) | ElementHash( INPUT) ;
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironment2::UnElemHash( const size_t &HASHEDATOM) const
    { // need to check the number of bits?
      std::string symbol( GetAtomicSymbol( HASHEDATOM & 127, false) );
      std::string bond_order( GetBondOrder( ( HASHEDATOM >> 7) & 3) );
      return ( bond_order + symbol);
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironment2::UnElemRCHash(const size_t &HASHEDATOM) const
    {
      std::ostringstream oss;
      oss << UnElemHash( HASHEDATOM) << GetRCSymbol( ( HASHEDATOM >> 9) & 3);
      return oss.str();
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironment2::AtomHash( const t_input &INPUT) const
    {
      return ( INPUT.Second() << 8) | GetAtomTypeIndex( INPUT.First() );
    }

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomEnvironment2::AtomRCHash( const t_input &INPUT) const
    {
      return ( AtomEnvironment2::RingOrAromatic( INPUT.First()) << 10) | AtomHash( INPUT) ;
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironment2::UnAtomHash( const size_t &HASHEDATOM) const
    { // need to check the number of bits?
      std::string symbol( GetAtomicSymbol( HASHEDATOM & 255, true) );
      std::string bond_order( GetBondOrder( ( HASHEDATOM >> 8) & 3) );
      return ( bond_order + symbol);
    }

    //! @brief converts the hashed atom into its string representation
    std::string AtomEnvironment2::UnAtomRCHash(const size_t &HASHEDATOM) const
    {
      std::ostringstream oss;
      oss << UnAtomHash( HASHEDATOM) << GetRCSymbol( ( HASHEDATOM >> 10) & 3);
      return oss.str();
    }

    //! @brief converts the hashed string back to the AE
    AtomEnvironment2::t_AtomEnvironment AtomEnvironment2::StringToAE(const std::string &STRING, const Atom_type ATOM_TYPE)
    {
      t_AtomEnvironment AE;
      storage::Vector< std::string> tokens(util::SplitString(STRING, "[]"));
      size_t count(0);
      for ( storage::Vector< std::string>::const_iterator iter = tokens.Begin(); iter != tokens.End(); ++iter, count++)
      {
        storage::Vector< std::string> atoms(util::SplitString(*iter, " "));
        storage::Map< size_t, size_t> layer;
        for ( storage::Vector< std::string>::const_iterator itr = atoms.Begin(); itr != atoms.End(); ++itr)
        {
          size_t hashed_atom;
          if (ATOM_TYPE == e_Element)
          {
            hashed_atom = ElementStringToAE(*itr);
          }
          else if (ATOM_TYPE == e_ElemRC)
          {
            hashed_atom = ElemRCStringToAE(*itr);
          }
          else if (ATOM_TYPE == e_Atom)
          {
            hashed_atom = AtomStringToAE(*itr);
          }
          else
          {
            hashed_atom = AtomRCStringToAE(*itr);
          }
          AddToMap<size_t>( hashed_atom, layer);
        }
        AE.PushBack(layer);
      }
      return AE;
    }

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironment2::ElementStringToAE(const std::string &STRING)
    {
      return CharToBond(STRING[0]) << 7 | StringToElement( STRING.substr( 1, std::string::npos));
    }

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironment2::ElemRCStringToAE(const std::string &STRING)
    {
      storage::Vector< std::string> tokens(util::SplitString(STRING, "."));
      return StringToRC(*(tokens[1])) << 9 | ElementStringToAE(*(tokens[0]));
    }

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironment2::AtomStringToAE(const std::string &STRING)
    {
      return CharToBond(STRING[0]) << 8 | StringToAtom( STRING.substr(1, std::string::npos));
    }

    //! @brief converts the hashed string back to the Element AE
    size_t AtomEnvironment2::AtomRCStringToAE(const std::string &STRING)
    {
      storage::Vector< std::string> tokens(util::SplitString(STRING, "."));
      return StringToRC( *( tokens[1])) << 10 | AtomStringToAE( *( tokens[0]));
    }

    //! @brief converts the string back to bond info
    size_t AtomEnvironment2::CharToBond(const char &CHAR)
    {
      if (CHAR == '-') return 1;
      if (CHAR == '=') return 2;
      if (CHAR == '#') return 3;
      return 0;
    }

    //! @brief converts the string back to aromatic/cyclic info
    size_t AtomEnvironment2::StringToRC(const std::string &STRING)
    {
      if (STRING == "c") return 1;
      if (STRING == "r") return 2;
      return 0;
    }

    //! @brief converts the string back to the Element info
    size_t AtomEnvironment2::StringToElement(const std::string &STRING)
    {
      return chemistry::GetElementTypes().ElementTypeLookup(STRING).GetIndex();
    }

    //! @brief converts the string back to the Atom type info
    size_t AtomEnvironment2::StringToAtom(const std::string &STRING)
    {
      return AtomType(STRING).GetIndex();
    }
  }
}
