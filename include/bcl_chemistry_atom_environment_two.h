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

#ifndef BCL_CHEMISTRY_ATOM_ENVIRONMENT_TWO_H_
#define BCL_CHEMISTRY_ATOM_ENVIRONMENT_TWO_H_
// include the namespace headers
#include "bcl_chemistry.h"
#include "util/bcl_util_serializable_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h" //!> input is a fragment complete
#include "storage/bcl_storage_vector.h" //!> for map between atom type and count of the atom types
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"
#include "util/bcl_util_function_wrapper.h"
#include "function/bcl_function_member_unary_const.h"
#include "storage/bcl_storage_triplet.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
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
    class BCL_API AtomEnvironment2 :
  public util::SerializableInterface
  {

    friend class MoleculeEnvironment;

  public:
    typedef storage::Vector< storage::Map< size_t, size_t> > t_AtomEnvironment;
    typedef storage::Vector< AtomEnvironment2> t_MoleculeEnv;
    // Different choices of atom types
    enum Atom_type
    {
      e_Element,  //!< Element atom type (atomic number, bond orders)
      e_ElemRC,   //!< Element RC atom type (atomic number, bond orders, is in ring?, is Aromatic?)
      e_Atom,     //!< Atom type (BCL Atom Type, bond orders, is in ring?, is Aromatic?)
      e_AtomRC,   //!< Atom type RC (BCL Atom Type, bond orders, is in ring?, is Aromatic?)
      s_NumberAtomTypes
    };

    //! @brief atom type as string
    //! @param ATOM_TYPE the name of the atom type
    //! @return the string for the atom type
    static const std::string &GetAtomTypeName( const Atom_type &ATOM_TYPE);

    //! @brief Initialization enum I/O helper
    typedef util::WrapperEnum< Atom_type, &GetAtomTypeName, s_NumberAtomTypes> AtomTypeEnum;

  private:

    typedef storage::Pair<const AtomConformationalInterface&, size_t> t_input;               //!> the wrapper parameter
    typedef storage::Pair< storage::Vector< std::string >, storage::Vector< size_t> > t_output;
    typedef storage::Pair
        <
        function::MemberUnaryConst<AtomEnvironment2, const  AtomEnvironment2::t_input&, size_t>,
        function::MemberUnaryConst<AtomEnvironment2, const size_t&, std::string>
    > t_functions;
    size_t m_AtomIndex;                                                       //!> the index of center of the atom environment
    size_t m_BondNumLimit;                                       //!> maximum number of bond distance
    AtomTypeEnum m_AtomType;                                                       //!> The type of atom type hash
    // the actual atom environment, size_t is the counts of the AtomType
    // Map of compressed atom type (size_t) and their counts (size_t)
    t_AtomEnvironment m_AtomEnvironment;

  public:

    //! @brief Get the name, description and function of the given atom type
    //! @param ATOM_TYPE the atom type used to build the fingerprints
    //! @return the short name or abbreviation of the class
    static const storage::Triplet
    <
    std::string,
    std::string,
    t_functions
    > &GetAtomTypeInfo( const Atom_type &ATOM_TYPE);

    //////////////////
    // Constructions//
    //////////////////

    //! default constructor
    AtomEnvironment2();

    //! constructor with Atom_type specification
    AtomEnvironment2( const Atom_type &ATOM_TYPE);

    //! constructor from AE string representation
    AtomEnvironment2( const Atom_type &ATOM_TYPE, const std::string &STRING);

    //! @brief constructor for building an atom environment from bond distance limit, index of atom, atom type, and fragment complete
    AtomEnvironment2 ( size_t ATOM_INDEX, const Atom_type &ATOM_TYPE, const ConformationInterface &FRAGMENT );

    //! virtual copy constructor
    AtomEnvironment2 *Clone() const;

    //! @brief compares two atom environments
    bool operator <( const AtomEnvironment2 &ATOM) const;

    //! @brief compares two atom environments
    bool operator ==(const AtomEnvironment2 &ATOM) const;

    //! @brief compares two atom environments
    bool operator !=(const AtomEnvironment2 &ATOM) const;

    //! instances of the class
    static const util::SiPtr< const util::ObjectInterface> s_Instances;

    //////////////////
    // data access ///
    //////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const;

    //! @brief access the atom of interest
    std::size_t GetAtomOfInterestIndex() const;

    //! @brief get the atom type that this atom enviroment calculator calculates
    const Atom_type &GetAtomType() const;

    //! @brief returns the data label
    //! @return data label as string
    const std::string &GetAlias() const;

    //! @brief unhashes the atom environment into its string representation
    std::string UnHash() const;

    //! @brief returns the atom environment
    const t_AtomEnvironment &GetAtomEnvironment() const;

  protected:

    //////////////////////
    // input and output //
    //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GetSerializer() const;

  private:

    //! Helper functions of the last constructor

    //! @brief read the fragment complete of the molecule and and the index of the atom of interest,
    //!        and then provide encoded version of the atom environment.
    AtomEnvironment2::t_output GetEncodedAtomEnvironment( const ConformationInterface &FRAGMENT);

    //! @brief recursively find the neighbor atoms in the atom environment of the atom of interest
    //!        and then update in corresponding character in the encoded atom environment
    //! @param BOND_NUM the number of bond away from the atom of interest
    //! @param ENCODED_ATOM_ENVIRONMENT the encoded atom enviroment, each string represents a phere,
    //!        each character represent absence/presece of the corresponding atom
    static void GetEncodedAtomEnvironmentHelper
    (
      const ConformationInterface &FRAGMENT, size_t ATOM_INDEX, size_t BOND_LIMIT,
      size_t BOND_NUM, storage::Vector< std::string> &ENCODED_ATOM_ENVIRONMENT,
      storage::Vector< size_t> &BOND_ORDERS
    );

    //! @brief returns the corresponding character representing the bond order
    static std::string GetBondOrder(size_t BOND_ORDER);

    //! @brief updates the encoded atom environment
    //! @param TARGET_ATOM_INDEX the index of the target atom, which bond to the atom of interest
    //! @param ENCODED_ATOM_ENVIRONMENT a vector of string that represents} the neighbor atoms in different spheres
    static void UpdateEncodedAtomEnvironment
    (
      size_t BOND_NUM, size_t TARGET_ATOM_INDEX,
      storage::Vector<std::string> &ENCODED_ATOM_ENVIRONMENT,
      size_t BOND_ORDER, storage::Vector<size_t> &BOND_ORDERS
    );

    //! @brief converts from the encoded atom environment to the normal representation of the atom environment
    //! @param ENCODED_ATOM_ENVIRONMENT encoded atom environment in term of a vector of string of 0 and 1
    void DecodeAtomEnvironment
    (
      const storage::Vector< std::string >&ENCODED_ATOM_ENVIRONMENT,
      const storage::Vector< size_t >&HASHED_MOLECULE,
      size_t BOND_LIMIT
    );

    //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
    template< typename K>
    static void AddToMap( const K &KEY, storage::Map<K, size_t> &MAP);

    // @brief converts vector of sorted objects into a map of objects and their counts
    static storage::Map<AtomEnvironment2, size_t> ConvertVetorToMap(const storage::Vector<AtomEnvironment2> &VEC);

    //! @brief adds the hashed atom into m_AtomEnvironment
    //! @param HASHED_ATOM
    //! @param BOND_DISTANCE
    void AddAtom( size_t HASHED_ATOM, size_t BOND_DISTANCE);

    //! @brief compute an atom environment for a atom of interest in a vector of AtomComplete
    //! @param FRAGMENT the fragment complete representing the molecule
    //! @param ATOM: the index of the atom of interest
    static const std::string &GetAtomicSymbol( size_t ATOM_NUMBER, bool ATOM_FLAG);

    //! @brief returns the atomic number for a atom of interest
    //! @param ATOM: the index of the atom of interest
    static const size_t GetAtomicNumber( const AtomComplete &ATOM);

    //! @brief returns the hashed atomtype index for a atom of interest
    //! @param ATOM: the index of the atom of interest
    static const size_t GetAtomTypeIndex( const AtomComplete &ATOM);

    //! @brief Checks if the atom is in a aromatic ring or not
    //! @return 1 if the atom in aromatic ring, 0 otherwise
    static const size_t RingOrAromatic( const AtomComplete &ATOM);

    //! @brief return the string symbol from the "ring or aromatic" information
    static const std::string GetRCSymbol( size_t RC);

    ///////////////////////////////
    // Hash and unhash functions //
    ///////////////////////////////

    //! @brief choose the hash function based on the choice of atom type, then hash every atoms of a fragment complete
    //! @param ATOM_TYPE the choice of the atom type
    //! @param Fragment the atom that is converted into the atom type
    //! @return a number which is a compressed version of the atom type
    storage::Vector< size_t > HashMolecule
    (
      const ConformationInterface &FRAGMENT,
      const storage::Vector< size_t> &BOND_ORDERS
    ) const;

    //! @brief hashes the Element atom type info of the ATOM
    //! @param ATOM the atom of interest
    //! @return the size_t hashed value of the ATOM
    size_t ElementHash( const t_input &INPUT) const;

    //! @brief hashes the Element atom type info of the ATOM
    //! @param ATOM the atom of interest
    //! @return the size_t hashed value of the ATOM
    size_t ElemRHash( const t_input &INPUT) const;

    //! @brief hashes the Element atom type info of the ATOM
    //! @param ATOM the atom of interest
    //! @return the size_t hashed value of the ATOM
    size_t ElemRCHash( const t_input &INPUT) const;

    //! @brief converts the hashed atom into its string representation
    std::string UnElemHash(const size_t &HASHEDATOM) const;

    //! @brief converts the hashed atom into its string representation
    std::string UnElemRCHash(const size_t &HASHEDATOM) const;

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomHash( const t_input &INPUT) const;

    //! @brief hashes the Element atom type info of the INPUT
    size_t AtomRCHash( const t_input &INPUT) const;

    //! @brief converts the hashed atom into its string representation
    std::string UnAtomHash( const size_t &HASHEDATOM) const;

    //! @brief converts the hashed atom into its string representation
    std::string UnAtomRCHash(const size_t &HASHEDATOM) const;

    //! @brief converts the hashed string back to the AE
    t_AtomEnvironment StringToAE(const std::string &STRING, const Atom_type ATOM_TYPE);

    //! @brief converts the hashed string back to the Element AE
    size_t ElementStringToAE(const std::string &STRING);

    //! @brief converts the hashed string back to the Element AE
    size_t ElemRCStringToAE(const std::string &STRING);

    //! @brief converts the hashed string back to the Element AE
    size_t AtomStringToAE(const std::string &STRING);

    //! @brief converts the hashed string back to the Element AE
    size_t AtomRCStringToAE(const std::string &STRING);

    //! @brief converts the char back to bond info
    size_t CharToBond(const char &CHAR);

    //! @brief converts the char back to aromatic/cyclic info
    size_t StringToRC(const std::string &STRING);

    //! @brief converts the string back to the Element info
    size_t StringToElement(const std::string &STRING);

    //! @brief converts the string back to the Atom type info
    size_t StringToAtom(const std::string &STRING);

  };
  }
}
#endif /* BCL_CHEMISTRY_ATOM_ENVIRONMENT_TWO_H_ */
